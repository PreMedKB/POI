#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
POI Analysis Modules for test.
Personal multi-omics data, including somatic variants, germline variants, gene fusion, 
expression status, copy number variants, msi and tmb are combined for drug recommendations.
"""


import pymysql, re, json, toml, time, os, sys, getopt
import pandas as pd
import numpy as np
from src.functions import connect_database, pymysql_cursor, load_logging_cfg
from src.direct import drug_direct
from src.repurpose import drug_repurpose
from panno.genotype_resolution import resolution
from panno.clinical_annotation import annotation
from panno.pgx_report import report
from src.summary import summary_report

# Connect database
cfile = './conf/config.toml'
cf = toml.load(cfile, _dict=dict)
DEFAULT = cf['DEFAULT']
db, cursor = connect_database(DEFAULT)
logger = load_logging_cfg(DEFAULT)

# define mkdir
def mkdir(path):
  path=path.rstrip("/")
  isExists=os.path.exists(path)
  if not isExists:
    os.makedirs(path)
    #print (path+'succeed')
    return True
  else:
    #print (path+'failed')
    return False

# get input

somatic_vcf=''
somatic_anno_txt=''
germline_vcf=''
germline_anno_txt=''
input_cnv=''
input_fusion=''
input_tmb=''
input_msi=''
tumor_exp=''
normal_exp=''
fp='' # fp, the directory of the data
version='v20221220'

# test (/mnt/poi-testdata)
Population = "East Asian"; RefGenome='hg38'
tissue = 'Breast'; disease = 'Metastatic Triple-Negative Breast Carcinoma'
sample_name='FUSCCTNBC008'
somatic_vcf='/mnt/poi-testdata/poi-tnbc/FUSCCTNBC008/FUSCCTNBC008.TNseq.filter.PASS.vcf'
somatic_anno_txt='/mnt/poi-testdata/poi-tnbc/FUSCCTNBC008/FUSCCTNBC008.TNseq.filter.PASS.hg38_multianno.txt'
germline_vcf='/mnt/poi-testdata/poi-tnbc/FUSCCTNBC008/FUSCCTNBC008.Haplotyper.PASS.vcf'
germline_anno_txt='/mnt/poi-testdata/poi-tnbc/FUSCCTNBC008/FUSCCTNBC008.Haplotyper.PASS.hg38_multianno.txt'
tumor_exp='/mnt/poi-testdata/poi-tnbc/FUSCCTNBC008/FUSCCTNBC008_exp.csv'
normal_exp='/mnt/poi-testdata/poi-tnbc/FUSCCTNBC008/FUSCCTNBC008_PT_exp.csv'
fp='/mnt/poi-testdata/poi-tnbc/FUSCCTNBC008/' 


input_cnv='/mnt/poi-testdata/poi-mskimpack/P-0000649-T01-IM3/P-0000649-T01-IM3_cnv.csv'
input_fusion='/mnt/poi-testdata/poi-mskimpack/P-0000649-T01-IM3/P-0000649-T01-IM3_fusion.tsv'
tmb='TMB-H'
msi='MSI-H'


############### BASIC INFORMATION #################
dic_front2race = {"African/African American": "afr", "Latino/Admixed American": "amr", "South Asian": "sas", "East Asian": "eas", "Non-Finnish European": "nfe", "Ashkenazi Jewish": "asj", "Finnish": "fin", "Other": "oth"}
race = dic_front2race[Population]
assembly = RefGenome

# Narrowing Therapy IDs based on the input disease
FinalResult = {'Sample': sample_name, 'Race': Population, 'Tissue': tissue, 'Disease': disease}

# Find all child diseases: child_id & parent_id
child_id = pymysql_cursor('SELECT MetaID FROM OT.DiseaseDisplay2Meta WHERE DisplayID IN (SELECT ChildID FROM OT.DiseaseDisplayChild WHERE DisplayID IN (SELECT ID FROM OT.DiseaseDisplay WHERE NormTissue = "%s" AND NormDisease = "%s"));' % (tissue, disease))
if type(child_id) != list:
  child_id = [child_id]

parent_id = pymysql_cursor('SELECT MetaID FROM OT.DiseaseDisplay2Meta WHERE DisplayID IN (SELECT DisplayID FROM OT.DiseaseDisplayChild WHERE ChildID IN (SELECT ID FROM OT.DiseaseDisplay WHERE NormTissue = "%s" AND NormDisease = "%s"));' % (tissue, disease))
if type(parent_id) != list:
  parent_id = [parent_id]

disease_id = list(set(child_id + parent_id))


############### LOAD FILES #################
Input = []; alterations_num = 0
status = 1; Input = []; alterations_num = 0
print(id)

# Path exist?
if not os.path.exists(fp):
  status = 3

t1 = time.time()
### Somatic VCF
if somatic_vcf:
  Input.append("Somatic variant")
  if os.path.exists(somatic_anno_txt) is False:
    ## Annotation by ANNOVAR
    command = 'bcftools norm -m -both %s -o %s' % (somatic_vcf, os.path.join(fp, 'SomaticVCF.norm.vcf'))
    os.system(command)
    command = '/mnt/annovar/table_annovar.pl %s /mnt/annovar/humandb -buildver %s -out %s -remove -protocol refGene,clinvar_20220320,gnomad211_exome,dbnsfp42c -operation g,f,f,f -nastring . -vcfinput -thread 2' % (os.path.join(fp, 'SomaticVCF.norm.vcf'), assembly, os.path.join(fp, 'Somatic'))
    os.system(command)
  # Read the ANNOVAR multianno txt
  svcf = pd.read_csv(somatic_anno_txt, sep="\t")
  somatic_anno = svcf#[svcf.Otherinfo10 == 'PASS']
elif somatic_vcf == '' or status == 3:
  somatic_vcf = None
  somatic_anno = pd.DataFrame()

### Germline VCF
if germline_vcf:
  Input.append("Germline variant")
  if os.path.exists(germline_anno_txt) is False:
    ## Annotation by ANNOVAR
    command = 'bcftools norm -m -both %s -o %s' % (germline_vcf, os.path.join(fp, 'GermlineVCF.norm.vcf'))
    os.system(command)
    command = '/mnt/annovar/table_annovar.pl %s /mnt/annovar/humandb -buildver %s -out %s -remove -protocol refGene,clinvar_20220320,gnomad211_exome,dbnsfp42c -operation g,f,f,f -nastring . -vcfinput -thread 2' % (os.path.join(fp, 'GermlineVCF.norm.vcf'), assembly, os.path.join(fp, 'Germline'))
    os.system(command)
  # Read the ANNOVAR multianno txt
  gvcf = pd.read_csv(germline_anno_txt, sep="\t")
  germline_anno = gvcf#[(gvcf.Otherinfo10 == 'PASS')]
elif germline_vcf == '' or status == 3:
  germline_vcf = None
  germline_anno = pd.DataFrame()

### CNV txt/tsv
if input_cnv:
  Input.append("CNV")
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
elif input_cnv == '' or status == 3:
  cnv = pd.DataFrame()

### Fusion txt/tsv
if input_fusion:
  Input.append("Gene fusion")
  fusion_df = pd.read_csv(input_fusion, sep="\t|,", engine='python', usecols=[0,1])
  fusion_df.columns = ["gene1", "gene2"]
  # Filtering the genes, which is not include in our database (OT)
  fusion_genes = open("./assets/fusion_genes.txt", "r").read().strip().split('\n')
  fusion = fusion_df[(fusion_df.gene1.isin(fusion_genes)) | (fusion_df.gene2.isin(fusion_genes))]
elif input_fusion == '' or status == 3:
  fusion = pd.DataFrame()

# ### TMB
# tmb = None
# if input_tmb:
#   tmb_df = pd.read_csv(input_tmb, sep="\t|,", engine='python', usecols=[0,1,2])
#   if float(tmb_df['tmb']) > 10:
#     Input.append("TMB")
#     tmb = 'TMB-H'
#   else:
#     tmb = None

# ### MSI
# msi = None
# if input_msi:
#   msi_df = pd.read_csv(input_msi, sep="\t|,", engine='python')
#   if float(msi_df['MSI.MANTIS.Score']) > 0.4:
#     Input.append("MSI")
#     msi = 'MSI-H'
#   else:
#     msi = None

### RNA gene expression
if tumor_exp:
  Input.append("Expression profile")
  # Load the input data and generate the rna_df.
  tt_df = pd.read_csv(tumor_exp, sep="\t|,", engine='python', usecols=[0,1])
  if tt_df.shape[1] == 1:
    tt_df = pd.DataFrame({'gene_id': tt_df.index.to_list(), 'TT': tt_df.iloc[:,0].to_list()})
  else:
    tt_df.columns = ['gene_id', 'TT']
  # Double quotes and version number
  tt_df.gene_id = tt_df.gene_id.str.strip('"|\'').str.split('.', expand=True)[0]
  if normal_exp:
    tp_df = pd.read_csv(normal_exp, sep="\t|,", engine='python', usecols=[0,1])
    if tp_df.shape[1] == 1:
      tp_df = pd.DataFrame({'gene_id': tp_df.index.to_list(), 'TP': tp_df.iloc[:,0].to_list()})
    else:
      tp_df.columns = ['gene_id', 'TP']
    # Double quotes and version number
    tp_df.gene_id = tp_df.gene_id.str.strip('"|\'').str.split('.', expand=True)[0]
    # Merge
    rna_exp = pd.merge(tt_df, tp_df, on='gene_id')
  else:
    rna_exp = tt_df
else:
  rna_exp = pd.DataFrame()

t2 = time.time()
preprocess_time = round(t2-t1, 3)


############### START THE ANALYSIS #################
##### POI: drug recommendation (direct)
if status != 3:
  print('DIRECT drug recommendation module...')
  try:
    assembly_dic = {'hg38': 'GRCh38', 'hg19': 'GRCh37'}
    if assembly.startswith('hg'):
      assembly = assembly_dic[assembly]
    (merged_therapy, alt_num, snpindel_biomarker, cnv_biomarker, fusion_biomarker, rna_biomarker) = drug_direct(tissue, disease, assembly, somatic_anno, germline_anno, cnv, fusion, tmb, msi, rna_exp, db, cursor)
    alterations_num = alterations_num + alt_num
    status = 1
  except Exception:
    status = 3

t3 = time.time(); direct_time = round(t3-t2,3)

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

t4 = time.time(); indirect_time = round(t4-t3,3)

##### POI: drug response (PAnno)
if status != 3:
  print('PANNO Drug response (pharmacogenomics) module...')
  try:
    if germline_vcf:
      panno_path = '/mnt/case/%s/panno.html' % sample_name
      dic_diplotype, dic_rs2gt, hla_subtypes = resolution(Population, germline_vcf)
      # Summary, PrescribingInfo, MultiLocus, SingleLocus, PhenotypePrediction, ClinicalAnnotation = annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
      # report(race, Summary, PrescribingInfo, MultiLocus, SingleLocus, PhenotypePrediction, ClinicalAnnotation, panno_path, sample_name)
      chemo_summary, prescribing_info, multi_var, single_var, phenotype_predict, clinical_anno = annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
      report(race, chemo_summary, prescribing_info, multi_var, single_var, phenotype_predict, clinical_anno, panno_path, sample_name)
      # Merge guidelines with clinical annotations
      chemotherapy = prescribing_info.merge(clinical_anno, how='left')
      multi_var = multi_var.rename({'Detected Alleles':'Variant Call'})
      single_var = single_var.rename({'Detected Alleles':'Variant Call'})
      # alterations_num = alterations_num + alt_num
    else:
      Summary = {'Avoid': [], 'Caution': [], 'Routine': []}
      chemotherapy = pd.DataFrame()
      multi_var = pd.DataFrame()
      single_var = pd.DataFrame()
    status = 1
  except Exception:
    status = 3

t5 = time.time(); response_time = round(t5-t4,3)

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
mkdir('/mnt/poi-2022test/FUSCCTNBC/%s/' % sample_name)
path = '/mnt/poi-2022test/FUSCCTNBC/%s/%s.%s.summary.json' % (sample_name, sample_name, version)
with open(path, "w") as f:
  json.dump(FinalResult, f, indent=2, ensure_ascii=False)

# Close connect
cursor.close()
db.close()