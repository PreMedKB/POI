#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
POI Analysis Modules.
Personal multi-omics data, including somatic variants, germline variants, gene fusion, 
expression status, copy number variants, msi and tmb are combined for drug recommendations.
"""


import pymysql, re, json, toml, time, os
import pandas as pd
import numpy as np
from src.functions import connect_database, pymysql_cursor, load_logging_cfg
from src.direct import drug_direct
from src.repurpose import drug_repurpose
from panno.genotype_resolution import resolution
from panno.clinical_annotation import annotation
from panno.pgx_report import report
from src.summary import summary_report


# Rotation of inquiries: If find the submit case, poi system will run.
def loop_func(func, second):
  while True:
    try:
      func()
      time.sleep(second)
    except Exception:
      continue


### Main
def main():
  # Connect database
  cfile = './conf/config.toml'
  cf = toml.load(cfile, _dict=dict)
  DEFAULT = cf['DEFAULT']
  db, cursor = connect_database(DEFAULT)
  logger = load_logging_cfg(DEFAULT)
  
  # Get task which status is waiting(0)
  info = cursor.execute('SELECT * FROM POI.Case WHERE Status = 0 LIMIT 1;')
  # info = cursor.execute('SELECT * FROM POI.Case WHERE Token = "example1";')
  info = cursor.fetchall()
  if info:
    ## Update status: 0 waiting; 1 running; 2 success; 3 failed
    status = 1
    headers = ['ID', 'Title', 'Token', 'Status', 'FilePath', 'DateTime', 'Tissue', 'Disease', 'Population', 'RefGenome', 'SomaticVCF', 'GermlineVCF', 'CNV', 'GeneFusion', 'TumorExp', 'NormalExp', 'TMB', 'MSI']
    case = dict(zip(headers, info[0]))
    token = case['Token']
    
    print('Find a case in the waiting status.\nToken is %s. Start analysing...' % token)
    pymysql_cursor('UPDATE POI.Case SET Status = "%s" WHERE Token = "%s";' % (status, token))
    
    FinalResult = {}
    t0=time.time()
    
    # Path exist?
    fp = case['FilePath']
    if not os.path.exists(fp):
      status = 3
    
    ############### BASIC INFORMATION #################
    dic_front2race = {"African/African American": "afr", "Latino/Admixed American": "amr", "South Asian": "sas", "East Asian": "eas", "Non-Finnish European": "nfe", "Ashkenazi Jewish": "asj", "Finnish": "fin", "Other": "oth"}
    race = dic_front2race[case['Population']]
    assembly = case['RefGenome']
    
    # Narrowing Therapy IDs based on the input disease
    tissue = case['Tissue']
    disease = case['Disease']
    FinalResult = {'Sample': token, 'Race': case['Population'], 'Tissue': tissue, 'Disease': disease}
    
    # Find all child diseases: child_id & parent_id
    child_id = pymysql_cursor('SELECT MetaID FROM OT.DiseaseDisplay2Meta WHERE DisplayID IN (SELECT ChildID FROM OT.DiseaseDisplayChild WHERE DisplayID IN (SELECT ID FROM OT.DiseaseDisplay WHERE NormTissue = "%s" AND NormDisease = "%s"));' % (tissue, disease))
    if type(child_id) != list:
      child_id = [child_id]
    
    parent_id = pymysql_cursor('SELECT MetaID FROM OT.DiseaseDisplay2Meta WHERE DisplayID IN (SELECT DisplayID FROM OT.DiseaseDisplayChild WHERE ChildID IN (SELECT ID FROM OT.DiseaseDisplay WHERE NormTissue = "%s" AND NormDisease = "%s"));' % (tissue, disease))
    if type(parent_id) != list:
      parent_id = [parent_id]
    
    disease_id = list(set(child_id + parent_id))
    
    ############### LOAD FILES #################
    try:
      files = os.listdir(fp)
      Input = []; alterations_num = 0
      ### Somatic VCF
      if case['SomaticVCF'] == 1:
        Input.append("somatic variant")
        for f in files:
          if re.findall('.*SomaticVCF.*', f):
            somatic_vcf = os.path.join(fp, re.findall('.*SomaticVCF.*', f)[0])
            input_somatic = os.path.join(fp, 'Somatic.%s_multianno.txt' % assembly)
            if os.path.exists(input_somatic) is False:
              ## Annotation by ANNOVAR
              command = 'bcftools norm -m -both %s -o %s' % (somatic_vcf, os.path.join(fp, 'SomaticVCF.norm.vcf'))
              os.system(command)
              command = '/mnt/annovar/table_annovar.pl %s /mnt/annovar/humandb -buildver %s -out %s -remove -protocol refGene,clinvar_20220320,gnomad211_exome,dbnsfp42c -operation g,f,f,f -nastring . -vcfinput -thread 2' % (os.path.join(fp, 'SomaticVCF.norm.vcf'), assembly, os.path.join(fp, 'Somatic'))
              os.system(command)
            # Read the ANNOVAR multianno txt
            svcf = pd.read_csv(input_somatic, sep="\t")
            somatic_anno = svcf#[svcf.Otherinfo10 == 'PASS']
        
      elif case['SomaticVCF'] == 0 or status == 3:
        somatic_vcf = None
        somatic_anno = pd.DataFrame()
      
      ### Germline VCF
      if case['GermlineVCF'] == 1:
        Input.append("germline variant")
        for f in files:
          if re.findall('.*GermlineVCF.*', f):
            germline_vcf = os.path.join(fp, re.findall('.*GermlineVCF.*', f)[0])
        input_germline = os.path.join(fp, 'Germline.%s_multianno.txt' % assembly)
        if os.path.exists(input_germline) is False:
          ## Annotation by ANNOVAR
          command = 'bcftools norm -m -both %s -o %s' % (germline_vcf, os.path.join(fp, 'GermlineVCF.norm.vcf'))
          os.system(command)
          command = '/mnt/annovar/table_annovar.pl %s /mnt/annovar/humandb -buildver %s -out %s -remove -protocol refGene,clinvar_20220320,gnomad211_exome,dbnsfp42c -operation g,f,f,f -nastring . -vcfinput -thread 2' % (os.path.join(fp, 'GermlineVCF.norm.vcf'), assembly, os.path.join(fp, 'Germline'))
          os.system(command)
        # Read the ANNOVAR multianno txt
        gvcf = pd.read_csv(input_germline, sep="\t")
        germline_anno = gvcf#[(gvcf.Otherinfo10 == 'PASS')]
      elif case['GermlineVCF'] == 0 or status == 3:
        germline_vcf = None
        germline_anno = pd.DataFrame()
      
      ### CNV txt/tsv
      if case['CNV'] == 1:
        Input.append("CNV")
        for f in files:
          if re.findall('.*CNV.*', f):
            input_cnv = os.path.join(fp, re.findall('.*CNV.*', f)[0])
            # Load file
            cnv_df = pd.read_csv(input_cnv, sep="\t|,", engine='python')#, usecols=[0,1,2])
            if cnv_df.shape[1] == 6:
              cnv_df = cnv_df.iloc[:,[0,1,3]]
              cnv_df.columns = ["symbol", "estimation", "copy_number"]
            else:
              cnv_df.columns = ["symbol", "estimation"]
            # Replace Chinese
            cnv_df = cnv_df.replace('正常', 'neutral').replace('缺失', 'loss').replace('扩增', 'gain')
            cnv_df.symbol = cnv_df.symbol.str.split('|', expand = True)[0]
            # Filtering the genes, which is not include in our database
            # cnv_genes = open("./assets/cnv_genes.txt", "r").read().strip().split('\n')
            symbol = pymysql_cursor('SELECT DISTINCT Symbol FROM Gene WHERE ID IN (SELECT ID from Gene WHERE ID IN (SELECT GeneID FROM Variant WHERE Name LIKE "%amplification%" OR Name LIKE "%loss%"));')
            alias = pymysql_cursor('SELECT DISTINCT Alias FROM GeneAlias WHERE length(Alias) < 10 AND GeneID IN (SELECT ID from Gene WHERE ID IN (SELECT GeneID FROM Variant WHERE Name LIKE "%amplification%" OR Name LIKE "%loss%"));')
            cnv_genes = list(set(symbol + alias))
            cnv = cnv_df[cnv_df.symbol.isin(cnv_genes)] 
      elif case['CNV'] == 0 or status == 3:
        cnv = pd.DataFrame()
      
      ### Fusion txt/tsv
      if case['GeneFusion'] == 1:
        Input.append("gene fusion")
        for f in files:
          if re.findall('.*GeneFusion.*', f):
            input_fusion = os.path.join(fp, re.findall('.*GeneFusion.*', f)[0])
            fusion_df = pd.read_csv(input_fusion, sep="\t|,", engine='python', usecols=[0,1])
            fusion_df.columns = ["gene1", "gene2"]
            # Filtering the genes, which is not include in our database (OT)
            fusion_genes = open("./assets/fusion_genes.txt", "r").read().strip().split('\n')
            fusion = fusion_df[(fusion_df.gene1.isin(fusion_genes)) | (fusion_df.gene2.isin(fusion_genes))]
      elif case['GeneFusion'] == 0 or status == 3:
        fusion = pd.DataFrame()
      
      ### TMB
      if case['TMB'] == 'yes':
        Input.append("TMB")
        tmb = 'TMB-H'
      else:
        tmb = None
      
      ### MSI
      if case['MSI'] == 'yes':
        Input.append("MSI")
        msi = 'MSI-H'
      else:
        msi = None
      
      ### RNA gene expression
      if case['TumorExp'] == 1:
        Input.append("gene expression")
        tumor_exp = None; normal_exp = None
        for f in files:
          if re.findall('.*TumorExp.*', f):
            tumor_exp = os.path.join(fp, re.findall('.*TumorExp.*', f)[0])
          if re.findall('.*NormalExp.*', f):
            normal_exp = os.path.join(fp, re.findall('.*NormalExp.*', f)[0])
        # Load the input data and generate the rna_df.
        tt_df = pd.read_csv(tumor_exp, sep="\t|,", engine='python', usecols=[0,1])
        if tt_df.iloc[0, 0].startswith("ENSG"):
          if tt_df.shape[1] == 1:
            tt_df = pd.DataFrame({'gene_id': tt_df.index.to_list(), 'TT': tt_df.iloc[:,0].to_list()})
          else:
            tt_df.columns = ['gene_id', 'TT']
          # Double quotes and version number
          tt_df.gene_id = tt_df.gene_id.str.strip('"|\'').str.split('.', expand=True)[0]
        else:
          tt_df.columns = ['gene_symbol', 'TT']
        if normal_exp:
          tp_df = pd.read_csv(normal_exp, sep="\t|,", engine='python', usecols=[0,1])
          if tp_df.iloc[0, 0].startswith("ENSG"):
            if tp_df.shape[1] == 1:
              tp_df = pd.DataFrame({'gene_id': tp_df.index.to_list(), 'TP': tp_df.iloc[:,0].to_list()})
            else:
              tp_df.columns = ['gene_id', 'TP']
            # Double quotes and version number
            tp_df.gene_id = tp_df.gene_id.str.strip('"|\'').str.split('.', expand=True)[0]
          else:
            tp_df.columns = ['gene_symbol', 'TP']
          # Merge
          rna_exp = pd.merge(tt_df, tp_df)
        else:
          rna_exp = tt_df
      else:
        rna_exp = pd.DataFrame()
    except Exception:
      status = 3
    
    t1 = time.time(); preprocess_time = round(t1-t0, 3)
    
    
    ############### START THE ANALYSIS #################
    ##### POI: drug recommendation (direct)
    if status != 3:
      print('DIRECT drug recommendation module...')
      try:
        assembly_dic = {'hg38': 'GRCh38', 'hg19': 'GRCh37'}
        if assembly.startswith('hg'):
          assembly = assembly_dic[assembly]
        (merged_therapy, alt_num, snpindel_biomarker, cnv_biomarker, fusion_biomarker, rna_biomarker) = drug_direct(tissue, disease, assembly, race, somatic_anno, germline_anno, cnv, fusion, tmb, msi, rna_exp, db, cursor)
        alterations_num = alterations_num + alt_num
        status = 1
      except Exception:
        status = 3
    t2 = time.time(); direct_time = round(t2-t1,3)
    
    ##### POI: drug repurposing (indirect)
    if status != 3:
      print('INDIRECT drug recommendation (repurposing) module...')
      try:
        if not somatic_anno.empty or not germline_anno.empty:
          (repurposing_therapy, alt_num, indirect_biomarker) = drug_repurpose(somatic_anno, germline_anno, race, disease_id, db, cursor)
          alterations_num = alterations_num + alt_num
        else:
          indirect_biomarker = pd.DataFrame()
          repurposing_therapy = pd.DataFrame()
        status = 1
        
      except Exception:
        status = 3
    t3 = time.time(); indirect_time = round(t3-t2,3)
    
    ##### POI: drug response (PAnno)
    if status != 3:
      print('PANNO Drug response (pharmacogenomics) module...')
      try:
        if germline_vcf:
          panno_path = '/mnt/case/%s/panno.html' % token
          dic_diplotype, dic_rs2gt, hla_subtypes = resolution(case['Population'], germline_vcf)
          if not dic_diplotype.empty:
            chemo_summary, prescribing_info, multi_var, single_var, phenotype_predict, clinical_anno = annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
            report(race, chemo_summary, prescribing_info, multi_var, single_var, phenotype_predict, clinical_anno, panno_path, token)
            # chemo_summary
            for key in chemo_summary.keys():
                for i in range(0, len(chemo_summary[key])):
                  chemo_summary[key][i] = chemo_summary[key][i].title()
            # Merge guidelines with clinical annotations
            chemotherapy = prescribing_info.merge(clinical_anno, how='left')
            multi_var = multi_var.rename({'Effect on Protein':'Effect_on_Protein', 'Definition of Alleles':'Definition_of_Alleles', 'Variant Call':'Variant_Call'})
            single_var = single_var.rename({'Effect on Protein':'Effect_on_Protein', 'Definition of Alleles':'Definition_of_Alleles', 'Variant Call':'Variant_Call'})
          else:
             chemo_summary = {'Avoid': [], 'Caution': [], 'Routine': []}
             chemotherapy = pd.DataFrame()
             multi_var = pd.DataFrame()
             single_var = pd.DataFrame()
        else:
          chemo_summary = {'Avoid': [], 'Caution': [], 'Routine': []}
          chemotherapy = pd.DataFrame()
          multi_var = pd.DataFrame()
          single_var = pd.DataFrame()
        status = 1
      except Exception:
        status = 3
    t4 = time.time(); response_time = round(t4-t3,3)
    
    ##### Get the POI report
    if status != 3:
      print('Reporting ...')
      
      try:
        (therapy_summary, detail, therapy_num) = summary_report(merged_therapy, repurposing_therapy, chemotherapy, snpindel_biomarker, cnv_biomarker, fusion_biomarker, rna_biomarker, indirect_biomarker, multi_var, single_var, disease_id, db, cursor)
        FinalResult["Input"] = ", ".join(Input)
        FinalResult["Alterations"] = alterations_num
        FinalResult["Therapies"] = therapy_num
        FinalResult['Summary'] = {"Target": therapy_summary, "Chemo": chemo_summary}
        FinalResult['Detail'] = detail
        status = 2
      except Exception:
        status = 3
    
    t5 = time.time(); summarize_time = round(t5-t4,3)
    FinalResult['Time'] = {'Preprocess': preprocess_time, 'Direct': direct_time, 'Indirect': indirect_time, 'Response': response_time, 'Summarize': summarize_time}
    FinalResult['Status'] = status
    
    ##### Export the POI report
    path = '/mnt/case/%s/summary.json' % token
    with open(path, "w") as f:
      json.dump(FinalResult, f, indent=2, ensure_ascii=False)
    
    # Compress the results
    #os.system('zip -r -j /mnt/case/%s/premedkb-poi-report.zip /mnt/case/%s/summary.json /mnt/case/%s/panno.html' % (token, token, token))
    os.system('zip -r -j /mnt/case/%s/premedkb-poi-report.zip /mnt/case/%s/summary.json' % (token, token))
    ## Update status: 0 waiting; 1 running; 2 success; 3 failed
    pymysql_cursor('UPDATE POI.Case SET Status = "%s" WHERE Token = "%s";' % (status, token))
    print("Finish This Task: %s. Final Status: %s.\n\n\n\n\n" % (token, status))
  
  # Close connect
  cursor.close()
  db.close()


if __name__ == "__main__":
  loop_func(main, 30)