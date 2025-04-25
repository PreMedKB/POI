#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pymysql, re, getopt, sys, vcf, json, toml, time, os
import pandas as pd
import numpy as np
from src.functions import pymysql_cursor, search_set, search_snpindel, search_others, parse_aachange


def drug_direct(tissue, disease, assembly, race, somatic_anno, germline_anno, cnv, fusion, tmb, msi, rna_exp, db, cursor):
  # Decipher the input disease, we also need to find the tissue info through input
  global detected_var
  detected_var = []
  # Find therapies separately
  if not somatic_anno.empty:
    print('Processing Somatic Annotated VCF or TXT...')
    somatic_res = parse_txt(somatic_anno, assembly, 'Somatic', race)
    somatic_therapy = somatic_res[0]
    detected_var.extend(somatic_res[1])
    somatic_biomarker = somatic_res[2]
  else:
    somatic_therapy = pd.DataFrame()
    somatic_biomarker = pd.DataFrame()
  
  if not germline_anno.empty:
    print('Processing Germline Annotated VCF or TXT...')
    # Filter by pathogenic
    germline_anno_p = germline_anno[germline_anno.CLNSIG.str.contains(r'Pathogenic|Likely_pathogenic', regex=True)]
    genelist = pymysql_cursor('SELECT DISTINCT Symbol FROM Variant WHERE ID IN (SELECT VariantID FROM Therapy WHERE ID IN (SELECT TherapyID FROM TherapyHasDetail WHERE DetailID IN (SELECT ID FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name REGEXP "Germline"))));')
    germline_res = parse_txt(germline_anno_p, assembly, 'Germline', race, genelist)
    germline_therapy = germline_res[0]
    detected_var.extend(germline_res[1])
    germline_biomarker = germline_res[2]
  else:
    germline_therapy = pd.DataFrame()
    germline_biomarker = pd.DataFrame()
  
  snpindel_biomarker = pd.concat([somatic_biomarker, germline_biomarker], axis=0)

  if not cnv.empty:
    print('Processing CNV...')
    cnv_res = parse_cnv(cnv, assembly)
    cnv_therapy = cnv_res[0]
    detected_var.extend(cnv_res[1])
    cnv_biomarker = cnv_res[2]
  else:
    cnv_therapy = pd.DataFrame()
    cnv_biomarker = pd.DataFrame()
  
  if not fusion.empty:
    print('Processing Fusion...')
    fusion_res = parse_fusion(fusion, assembly)
    fusion_therapy = fusion_res[0]
    detected_var.extend(fusion_res[1])
    fusion_biomarker = fusion_res[2]
  else:
    fusion_therapy = pd.DataFrame()
    fusion_biomarker = pd.DataFrame()
  
  if tmb:
    print('Processing TMB...')
    tmb_res = parse_tmb_msi(tmb, assembly)
    tmb_therapy = tmb_res[0]
    detected_var.extend(tmb_res[1])
  else:
    tmb_therapy = pd.DataFrame()
  
  if msi:
    print('Processing MSI...')
    msi_res = parse_tmb_msi(msi, assembly)
    msi_therapy = msi_res[0]
    detected_var.extend(msi_res[1])
  else:
    msi_therapy = pd.DataFrame()
  
  if not rna_exp.empty:
    print('Processing RNA...')
    rna_res = parse_rna(rna_exp, tissue, disease, assembly)
    rna_therapy = rna_res[0]
    detected_var.extend(rna_res[1])
    rna_biomarker = rna_res[2]
  else:
    rna_therapy = pd.DataFrame()
    rna_biomarker = pd.DataFrame()
  
  # Find therapies through `Set`
  # Process the detected_var to confirm every item is unique id
  print('Deciphering multiple interrelated biomarkers...')
  detected_var_df = pd.DataFrame(detected_var, columns=['POID', 'VariantID', 'Gene', 'Alteration', 'Source']).drop_duplicates().dropna().reset_index(drop=True)
  alt_num = detected_var_df[detected_var_df.POID != ''].POID.drop_duplicates().shape[0]
  if detected_var_df.empty is False:
    detected_var_df = detected_var_df.drop(['VariantID'], axis=1).join(detected_var_df['VariantID'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('VariantID')).reset_index(drop=True)
    set_therapy = search_set(detected_var_df, assembly)
  else:
    set_therapy = pd.DataFrame()
  
  ### Merge the result
  merged_therapy = pd.concat([somatic_therapy, germline_therapy, cnv_therapy, fusion_therapy, tmb_therapy, msi_therapy, rna_therapy, set_therapy], axis = 0)
#  print(merged_therapy)
  return (merged_therapy, alt_num, snpindel_biomarker, cnv_biomarker, fusion_biomarker, rna_biomarker)



###  Parse Input Files & Expert Therapies
def parse_txt(anno, assembly, source, race, genelist=None):
  vcf_therapy = pd.DataFrame()
  metavariants = []
  if genelist:
    anno = anno[anno['Gene.refGene'].isin(genelist)].reset_index(drop=True)
  else:
    anno = anno
  
  # Add columns: POID, Function
  anno = anno.reset_index(drop=True)
  anno.insert(0, 'POID', ['%s.%s' % (source, str(i+1)) for i in anno.index.to_list()])
  # ExonicFunc.refGene for splicing, intronic, UTR, etc., are all '.'.
  anno.insert(1, 'Function', anno[['ExonicFunc.refGene']])
  splicing_index = anno[anno['Func.refGene'].str.contains('splicing')].index
  anno.loc[splicing_index, 'Function'] = 'splicing'
  # For wildtype and whole gene deletions
  wild_type_genes = pymysql_cursor('SELECT DISTINCT Symbol FROM Variant WHERE Name LIKE "%Wildtype%";')
  wt_del = []
  for gene in wild_type_genes:
    if anno[anno['Gene.refGene'] == gene].empty:
      wt_del.append([gene, 'Wildtype', ''])
  
  # Parse every single variant and find the therapies
  for index, row in anno.iterrows():
    gene = row['Gene.refGene']
    # chr, pos, ref, alt = row['#CHROM'], row['POS'], row['REF'], row['ALT']
    chr, pos, ref, alt = row['Chr'], row['Start'], row['Ref'], row['Alt']
    # Whole gene deletion, add it into wt_del
    if re.search('wholegene', row['AAChange.refGene']) and re.search('deletion', row['ExonicFunc.refGene']):
      wt_del.append([gene, row['ExonicFunc.refGene'], row.POID])
    # Small variants
    else:
      nm_exon_cc_pc = parse_aachange(row['Gene.refGene'], row['AAChange.refGene'], row['GeneDetail.refGene']) # return [] for splicing
      fuc = row.Function
      # ClinVar
      if '|' in row['CLNSIG']:
        clinsig_match = re.compile(r'Pathogenic|Likely_pathogenic').findall(row['CLNSIG'])
        if clinsig_match:
          clinsig = clinsig_match[0]
        else:
          clinsig = row['CLNSIG'].split('|')[0]
      else:
        clinsig = row['CLNSIG']
      
      # Split the multiple protein changes
      var_multi = []
      for i in range(0, len(nm_exon_cc_pc)):
        var_multi.append([chr, pos, ref, alt, nm_exon_cc_pc[i][0], nm_exon_cc_pc[i][1], nm_exon_cc_pc[i][2], nm_exon_cc_pc[i][3], fuc, clinsig])
      if fuc == 'splicing' and nm_exon_cc_pc != []:
        var_multi.append([chr, pos, ref, alt, nm_exon_cc_pc[i][0], nm_exon_cc_pc[i][1], nm_exon_cc_pc[i][2], nm_exon_cc_pc[i][3], fuc, clinsig])
      ### Search database to get drugs
      for var in var_multi:
        res = search_snpindel(gene, var, source, assembly)
        therapy = res[0]
        if therapy.empty is False:
          # therapy.insert(0, 'FUSCC', "|".join([str(i) for i in var]))
          therapy.insert(0, 'POID', row.POID)
          vcf_therapy = pd.concat([vcf_therapy, therapy], axis=0)
#          path = '/mnt/case/%s/vcf_therapy.xlsx' % token
#          vcf_therapy.to_excel(path, index=False)
        if res[1]:
          metavariants.append((row.POID,) + res[1])
  
  # Parse wildtypes and wholegene deletions
  for gene, var, POID in wt_del:
    res = search_snpindel(gene, var, source, assembly)
    therapy = res[0]
    if therapy.empty is False:
      # therapy.insert(0, 'FUSCC', var)
      therapy.insert(0, 'POID', POID)
      vcf_therapy = pd.concat([vcf_therapy, therapy], axis=0)
#      path = '/mnt/case/%s/vcf_therapy1.xlsx' % token
#      vcf_therapy.to_excel(path, index=False)
    if res[1]:
      metavariants.append((POID,) + res[1])
  
  # Biomarker details
  snpindel_var = pd.DataFrame(metavariants, columns=['POID', 'VariantID', 'Gene', 'Alteration', 'Source'])
  snpindel_biomarker = anno[anno.POID.isin(snpindel_var.POID.to_list())]
  snpindel_biomarker = pd.merge(snpindel_biomarker, snpindel_var)
  cols = ['POID', 'Gene', 'Origin', 'Location', 'Ref>Alt', 'Transcript', 'Exon', 'Protein_Change', 'Function', 'Clinical_Significance', 'Pathogenic_Prediction', 'Population_AF']
  if snpindel_biomarker.empty:
    snpindel_biomarker = pd.DataFrame(columns=cols)
  else:
    snpindel_biomarker[['Location', 'Ref>Alt', 'Transcript', 'Exon', 'Protein_Change']] = snpindel_biomarker.Alteration.str.split(',', expand=True)
    ref_pop = 'AF_%s' % race
    snpindel_biomarker = snpindel_biomarker.copy().rename(columns={'Source':'Origin', 'CLNSIG':'Clinical_Significance', ref_pop: 'Population_AF'})
    snpindel_biomarker.insert(0, 'Pathogenic_Prediction', 'Tolerated')
#    d_index = snpindel_biomarker[(snpindel_biomarker.SIFT_pred=="D")|(snpindel_biomarker.LRT_pred=="D")|(snpindel_biomarker.MutationTaster_pred=="A")|(snpindel_biomarker.MutationTaster_pred=="D")|(snpindel_biomarker.MutationAssessor_pred=="H")|(snpindel_biomarker.MutationAssessor_pred=="M")|(snpindel_biomarker.FATHMM_pred=="D")|(snpindel_biomarker.PROVEAN_pred=="D")|(snpindel_biomarker.MetaSVM_pred=="D")|(snpindel_biomarker.MetaLR_pred=="D")].index
#    snpindel_biomarker.loc[d_index, 'Pathogenic_Prediction'] = 'Deleterious'
    snpindel_biomarker = snpindel_biomarker[cols]
    
  return(vcf_therapy, metavariants, snpindel_biomarker)


def parse_cnv(cnv, assembly):
  cnv_therapy = pd.DataFrame()
  metavariants = []
  
  cnv = cnv.reset_index(drop=True)
  cnv.insert(0, 'POID', ['%s.%s' % ("CNV", str(i+1)) for i in cnv.index.to_list()])
  
  # Row-by-row to find therapy
  for index, row in cnv.iterrows():
    # Search database get drug
    gene = row["symbol"].split('|')[0]
    var = row["estimation"]
    res = search_others(gene, var, "CNV", assembly)
    therapy = res[0]
    if therapy.empty is False:
      # if len(row) == 3:
      #   therapy.insert(0, 'FUSCC', row.copy_number)
      therapy.insert(0, 'POID', row.POID)
      cnv_therapy = pd.concat([cnv_therapy, therapy], axis=0)
    if res[1]:
      metavariants.append((row.POID,) + res[1])
  
  # Biomarker details
  cnv_var = pd.DataFrame(metavariants, columns=['POID', 'VariantID', 'Gene', 'Alteration', 'Source'])
  cnv = cnv.rename(columns = {'symbol':'Gene', 'estimation':'CNV_State'})
  cnv_biomarker = pd.concat([cnv, pd.DataFrame(columns=['Copy_Change', 'LOH_State', 'Cytoband', 'Location'])], axis=1)
  cnv_biomarker = cnv_biomarker[cnv_biomarker.POID.isin(cnv_var.POID.to_list())]
  cnv_biomarker = cnv_biomarker[['POID', 'Gene', 'Copy_Change', 'CNV_State', 'LOH_State', 'Cytoband', 'Location']]
  
  return(cnv_therapy, metavariants, cnv_biomarker)


def parse_fusion(fusion, assembly):
  fusion_therapy = pd.DataFrame()
  metavariants = []
  
  fusion = fusion.reset_index(drop=True)
  fusion.insert(0, 'POID', ['%s.%s' % ("Fusion", str(i+1)) for i in fusion.index.to_list()])
  
  # Row-by-row to find therapy
  for index, row in fusion.iterrows():
    # Search database get drug
    gene = sorted([row.gene1, row.gene2])
    var = ["%s-%s Fusion" % (gene[0], gene[1]), "%s Fusions" % gene[0], "%s Fusions" % gene[1]]
    res = search_others(gene, var, "Fusion", assembly)
    therapy = res[0]
    if therapy.empty is False:
      therapy.insert(0, 'POID', row.POID)
      fusion_therapy = pd.concat([fusion_therapy, therapy], axis=0)
    if res[1]:
      metavariants.append((row.POID,) + res[1])
  
  # Biomarker details
  fusion_var = pd.DataFrame(metavariants, columns=['POID', 'VariantID', 'Gene', 'Alteration', 'Source'])
  fusion_biomarker = pd.concat([fusion, pd.DataFrame(columns=['Breakpoint', 'Event_Type', 'Sample', 'Cytogenic_Description'])], axis=1)
  fusion_biomarker = fusion_biomarker[fusion_biomarker.POID.isin(fusion_var.POID.to_list())]
  fusion_biomarker = fusion_biomarker[['POID', 'gene1', 'gene2', 'Breakpoint', 'Event_Type', 'Sample', 'Cytogenic_Description']]
  
  return(fusion_therapy, metavariants, fusion_biomarker)


def parse_tmb_msi(tmb_or_msi, assembly):
  metavariants = []
  gene = "Other Biomarkers"
  var_type = tmb_or_msi.split('-')[0]
  res = search_others(gene, tmb_or_msi, var_type, assembly)
  therapy = res[0]
  therapy.insert(0, 'POID', tmb_or_msi)
  metavariants.append((tmb_or_msi,) + res[1])
  
  return(therapy, metavariants)


def parse_rna(rna_exp, tissue, disease, assembly):
  rna_exp = rna_exp.reset_index(drop=True)
  rna_exp.insert(0, 'POID', ['%s.%s' % ("RNA", str(i+1)) for i in rna_exp.index.to_list()])
  
  # Preprocessing.
  # Import gene length information
  # Transform the format of expression: count2CPM
  gene_e = rna_exp
  gene_e['CPM_TT']=gene_e['TT']*1000000/rna_exp['TT'].sum()
  if 'TP' in gene_e.columns:
    gene_e['CPM_TP']=gene_e['TP']*1000000/rna_exp['TP'].sum()
  
  # Transform the ID into Symbol.
  geneid2name = pd.read_csv("./assets/geneid2name.csv", sep=",")
  if 'gene_id' in gene_e.columns:
    rna_df_symbol = pd.merge(geneid2name, gene_e, on='gene_id')
    # Summation of duplicate values by gene symbol, keep the max value
    rna_df_symbol = rna_df_symbol.sort_values(['gene_symbol', 'TT'], ascending=[False, False])
    rna_symbol = rna_df_symbol.drop_duplicates(subset='gene_symbol', keep='first')
  else:
    rna_df_symbol=gene_e.sort_values(['gene_symbol', 'TT'], ascending=[False, False])
    rna_symbol = rna_df_symbol.drop_duplicates(subset='gene_symbol', keep='first')
  
  # Filtering. raw read count >= 2 & tmp >= 0.5 are retained.
  # reference: https://f1000research.com/articles/5-1438 #Filtering to remove low counts
  if 'TP' in gene_e.columns:
    rna_anno = rna_symbol[(rna_symbol[["TT", "TP"]].T >= 10).all()]
    converted_df = rna_anno[(rna_anno[["CPM_TT", "CPM_TP"]].T >= 0.5).all()]
  
  # Determine if the excluded genes are drug-related.
  rna_targets = pymysql_cursor('SELECT DISTINCT Symbol FROM MetaLite.Gene WHERE ID IN (SELECT GeneID FROM MetaLite.Variant WHERE Name REGEXP "expression");')
  # exp_genes = open("./assets/exp_genes.txt", "r").read().strip().split('\n')
  excluded_genes = list(set(rna_targets) - set(converted_df['gene_symbol'].to_list()))
  # under_exp = pd.DataFrame(excluded_genes, columns=['gene_symbol'])
  # under_exp["Alteration"] = "Underexpression"
  print('Total number of drug-related genes: %s.\nNumber of filtered genes in the previous step: %s.' % (len(rna_targets), len(excluded_genes)))
  # Calculate WINTER score: TP exists.
  if "TP" in converted_df.columns:
    # Get log2(fold-change)
    converted_df.insert(converted_df.shape[1], 'log2fc', np.log2(converted_df.CPM_TT+0.01) - np.log2(converted_df.CPM_TP+0.01))
    # Judge the expression level
    low_exp = converted_df[converted_df['log2fc'] <= -3.5] # which should use the threshold from the database
    low_exp.insert(low_exp.shape[1], 'Alteration', "Underexpression")
    over_exp = converted_df[converted_df['log2fc'] >= 3.5]
    over_exp.insert(over_exp.shape[1], 'Alteration', "Overexpression")
    # Merge the above three dataframe
    # biomarkers = pd.concat([low_exp[['gene_symbol', 'Alteration']], over_exp[['gene_symbol', 'Alteration']]], axis=0)
    biomarkers = pd.concat([low_exp, over_exp], axis=0)
  # regard those overespression genes as expression genes
  exp = over_exp.replace({"Alteration":{"Overexpression":"Expression"}})
  biomarkers_final = pd.concat([biomarkers, exp], axis=0)
  
  # get CPM from RNA database; calculate the percentile of outlier gene in the patient
  outlier = set(biomarkers_final['gene_symbol'][biomarkers_final.gene_symbol.isin(rna_targets)].to_list())
  outlier_info = pd.DataFrame(columns=['gene_symbol', 'pecentile', 'rank_c'])
  for i in outlier:
    cpm_list = pymysql_cursor('SELECT CPMTumor FROM RNA.Expression WHERE Symbol = "%s" AND ID IN (SELECT ExpressionID FROM RNA.Expression2TCGA WHERE TCGAID IN (SELECT ID FROM RNA.TCGA WHERE NormTissue = "%s"));' % (i, tissue))
    cpm_p = biomarkers_final['CPM_TT'][biomarkers_final.gene_symbol==i].to_list()
    cpm_a = (cpm_list + cpm_p); cpm_a.sort()
    rank = cpm_a.index(cpm_p[0])+1; pecentile = rank/len(cpm_a)*100
    rank_c = str(rank) + "/" + str(len(cpm_a))
    temp = pd.DataFrame({'gene_symbol': [i],'pecentile': [pecentile],'rank_c': [rank_c]})
    outlier_info = pd.concat([outlier_info,temp],ignore_index=True)

  # the dataframe for Biomarker Detail
  biomarkers_show = pd.merge(biomarkers_final, outlier_info,on = 'gene_symbol')
  # 2 digits
  biomarkers_show['CPM_TT'] = round(biomarkers_show['CPM_TT'],2)
  biomarkers_show['CPM_TP'] = round(biomarkers_show['CPM_TP'],2)
  biomarkers_show['pecentile'] = round(biomarkers_show['pecentile'],2)
  
  # Find therapies
  rna_therapy = pd.DataFrame()
  metavariants = []
  
  biomarkers_final = biomarkers_final[biomarkers_final.gene_symbol.isin(rna_targets)].reset_index(drop=True)
  for index, row in biomarkers_final.iterrows():
    gene = row.gene_symbol
    var = row.Alteration
    res = search_others(gene, var, "RNA", assembly)
    therapy = res[0]
    if therapy.empty is False:
      therapy.insert(0, 'POID', row.POID)
      rna_therapy = pd.concat([rna_therapy, therapy], axis=0)
    if res[1]:
      metavariants.append((row.POID,) + res[1])
  
  # Biomarker details
  rna_var = pd.DataFrame(metavariants, columns=['POID', 'VariantID', 'Gene', 'Alteration', 'Source'])
  biomarkers_show = biomarkers_show.rename(columns={'gene_symbol':'Gene', 'CPM_TT':'Tumor_CPM', 'CPM_TP':'Normal_CPM', 'log2fc':'log2FC', 'pecentile':'Percentile', 'rank_c':'Rank' })
  rna_biomarker = pd.merge(rna_var, biomarkers_show)
  rna_biomarker = rna_biomarker[['POID', 'Gene', 'Alteration', 'Tumor_CPM', 'Normal_CPM', 'log2FC', 'Rank', 'Percentile']]
  
  return(rna_therapy, metavariants, rna_biomarker)
