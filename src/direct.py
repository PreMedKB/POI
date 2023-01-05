#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pymysql, re, getopt, sys, vcf, json, toml, time, os
import pandas as pd
import numpy as np
from src.functions import pymysql_cursor, search_set, search_snpindel, search_others, parse_aachange


def drug_direct(tissue, disease, assembly, somatic_anno, germline_anno, cnv, fusion, tmb, msi, rna_exp, db, cursor):
  # Decipher the input disease, we also need to find the tissue info through input
  global detected_var
  detected_var = []
  
  # Find therapies separately
  if not somatic_anno.empty:
    print('Processing Somatic Annotated VCF or TXT...')
    somatic_res = parse_txt(somatic_anno, assembly, 'Somatic')
    somatic_therapy = somatic_res[0]
    detected_var.extend(somatic_res[1])
  else:
    somatic_therapy = pd.DataFrame()
  
  if not germline_anno.empty:
    print('Processing Germline Annotated VCF or TXT...')
    # Filter by pathogenic
    germline_anno_p = germline_anno[germline_anno.CLNSIG.str.contains(r'Pathogenic|Likely_pathogenic', regex=True)]
    genelist = pymysql_cursor('SELECT DISTINCT Symbol FROM Variant WHERE ID IN (SELECT VariantID FROM Therapy WHERE ID IN (SELECT TherapyID FROM TherapyHasDetail WHERE DetailID IN (SELECT ID FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name REGEXP "Germline"))));')
    germline_res = parse_txt(germline_anno_p, assembly, 'Germline', genelist)
    germline_therapy = germline_res[0]
    detected_var.extend(germline_res[1])
  else:
    germline_therapy = pd.DataFrame()
  
  if not cnv.empty:
    print('Processing CNV...')
    cnv_res = parse_cnv(cnv, assembly)
    cnv_therapy = cnv_res[0]
    detected_var.extend(cnv_res[1])
  else:
    cnv_therapy = pd.DataFrame()
  
  if not fusion.empty:
    print('Processing Fusion...')
    fusion_res = parse_fusion(fusion, assembly)
    fusion_therapy = fusion_res[0]
    detected_var.extend(fusion_res[1])
  else:
    fusion_therapy = pd.DataFrame()
  
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
  else:
    rna_therapy = pd.DataFrame()
  
  # Find therapies through `Set`
  # Process the detected_var to confirm every item is unique id
  print('Deciphering multiple interrelated biomarkers...')
  detected_var_df = pd.DataFrame(detected_var, columns=['VariantID', 'Gene', 'Alteration', 'Source']).dropna()
  alt_num = len(detected_var)
  if detected_var_df.empty is False:
    detected_var_df = detected_var_df.drop(['VariantID'], axis=1).join(detected_var_df['VariantID'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('VariantID'))
    set_therapy = search_set(detected_var_df, assembly)
  else:
    set_therapy = pd.DataFrame()
  
  
  ### Merge the result
  merged_therapy = pd.concat([somatic_therapy, germline_therapy, cnv_therapy, fusion_therapy, tmb_therapy, msi_therapy, rna_therapy, set_therapy], axis = 0)
  return (merged_therapy, alt_num)



###  Parse Input Files & Expert Therapies
def parse_txt(anno, assembly, source, genelist=None):
  vcf_therapy = pd.DataFrame()
  metavariants = []
  if genelist:
    anno = anno[anno['Gene.refGene'].isin(genelist)].reset_index(drop=True)
  else:
    anno = anno
  
  # For wildtype and whole gene deletions
  wild_type_genes = pymysql_cursor('SELECT DISTINCT Symbol FROM Variant WHERE Name LIKE "%Wildtype%";')
  wt_del = []
  for gene in wild_type_genes:
    if anno[anno['Gene.refGene'] == gene].empty:
      wt_del.append([gene, 'Wildtype'])
  
  # Parse every single variant and find the therapies
  for index, row in anno.iterrows():
    gene = row['Gene.refGene']
    # chr, pos, ref, alt = row['#CHROM'], row['POS'], row['REF'], row['ALT']
    chr, pos, ref, alt = row['Chr'], row['Start'], row['Ref'], row['Alt']
    # Whole gene deletion, add it into wt_del
    if re.search('wholegene', row['AAChange.refGene']) and re.search('deletion', row['ExonicFunc.refGene']):
      wt_del.append([gene, row['ExonicFunc.refGene']])
    # Small variants
    else:
      nm_exon_cc_pc = parse_aachange(row['AAChange.refGene'], row['GeneDetail.refGene']) # return [] for splicing
      # ClinVar
      if '|' in row['CLNSIG']:
        clinsig_match = re.compile(r'Pathogenic|Likely_pathogenic').findall(row['CLNSIG'])
        if clinsig_match:
          clinsig = clinsig_match[0]
        else:
          clinsig = row['CLNSIG'].split('|')[0]
      else:
        clinsig = row['CLNSIG']
      # ExonicFunc.refGene for splicing, intronic, UTR, etc., are all '.'.
      fuc = 'splicing' if 'splicing' in row['Func.refGene'] else row['ExonicFunc.refGene']
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
          therapy.insert(0, 'FUSCC', "|".join([str(i) for i in var]))
          vcf_therapy = pd.concat([vcf_therapy, therapy], axis=0)
        if res[1]:
          metavariants.append(res[1])
  
  # Parse wildtypes and wholegene deletions
  for gene, var in wt_del:
    res = search_snpindel(gene, var, source, assembly)
    therapy = res[0]
    if therapy.empty is False:
      therapy.insert(0, 'FUSCC', var)
      vcf_therapy = pd.concat([vcf_therapy, therapy], axis=0)
    if res[1]:
      metavariants.append(res[1])
  
  return(vcf_therapy, metavariants)


def parse_cnv(cnv, assembly):
  cnv_therapy = pd.DataFrame()
  metavariants = []
  # Row-by-row to find therapy
  for index, row in cnv.iterrows():
    # Search database get drug
    gene = row["symbol"].split('|')[0]
    var = row["estimation"]
    res = search_others(gene, var, "CNV", assembly)
    therapy = res[0]
    if therapy.empty is False:
      if len(row) == 3:
        therapy.insert(0, 'FUSCC', row.copy_number)
      cnv_therapy = pd.concat([cnv_therapy, therapy], axis=0)
    if res[1]:
      metavariants.append(res[1])
  
  return(cnv_therapy, metavariants)


def parse_fusion(fusion, assembly):
  fusion_therapy = pd.DataFrame()
  metavariants = []
  # Row-by-row to find therapy
  for index, row in fusion.iterrows():
    # Search database get drug
    gene = sorted([row.gene1, row.gene2])
    var = ["%s-%s Fusion" % (gene[0], gene[1]), "%s Fusions" % gene[0], "%s Fusions" % gene[1]]
    res = search_others(gene, var, "Fusion", assembly)
    therapy = res[0]
    if therapy.empty is False:
      fusion_therapy = pd.concat([fusion_therapy, therapy], axis=0)
    if res[1]:
      metavariants.append(res[1])
  
  return(fusion_therapy, metavariants)


def parse_tmb_msi(tmb_or_msi, assembly):
  metavariants = []
  gene = "Other Biomarkers"
  var_type = tmb_or_msi.split('-')[0]
  res = search_others(gene, tmb_or_msi, var_type, assembly)
  therapy = res[0]
  metavariants.append(res[1])
  
  return(therapy, metavariants)


def parse_rna(rna_exp, tissue, disease, assembly):
  rna_therapy = pd.DataFrame()
  metavariants = []
  # Preprocessing.
  # Import gene length information
  # Transform the format of expression: count2CPM
  gene_e = rna_exp
  gene_e['CPM_TT']=gene_e['TT']*1000000/rna_exp['TT'].sum()
  if 'TP' in gene_e.columns:
    gene_e['CPM_TP']=gene_e['TP']*1000000/rna_exp['TP'].sum()
  
  # Transform the ID into Symbol.
  geneid2name = pd.read_csv("./assets/geneid2name.csv", sep=",")
  rna_df_symbol = pd.merge(geneid2name, gene_e, on='gene_id')
  # Summation of duplicate values by gene symbol, keep the max value
  rna_df_symbol = rna_df_symbol.sort_values(['gene_symbol', 'TT'], ascending=[False, False])
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
    biomarkers = pd.concat([low_exp[['gene_symbol', 'Alteration']], over_exp[['gene_symbol', 'Alteration']]], axis=0)
  # regard those overespression genes as expression genes
  exp = over_exp
  biomarkers_final = pd.concat([biomarkers, exp], axis=0)
  # Find therapies
  biomarkers_final = biomarkers_final[biomarkers_final.gene_symbol.isin(rna_targets)].reset_index(drop=True)
  for index, row in biomarkers_final.iterrows():
    gene = row.gene_symbol
    var = row.Alteration
    res = search_others(gene, var, "RNA", assembly)
    therapy = res[0]
    if therapy.empty is False:
      rna_therapy = pd.concat([rna_therapy, therapy], axis=0)
    if res[1]:
      metavariants.append(res[1])
  
  return(rna_therapy, metavariants)