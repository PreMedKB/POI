#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Drug repurposing by genomic variants.
"""

import re, vcf, json, toml
import pandas as pd
import numpy as np
from src.functions import connect_database, load_logging_cfg, pymysql_cursor, get_metaid, parse_aachange

def drug_repurpose(somatic_anno, germline_anno, race, disease_id, db, cursor):
  
  repurposing_therapy = pd.DataFrame()
  
  ## Related gene symbols
  related_target = []
  target = pymysql_cursor('SELECT DISTINCT Symbol FROM Gene WHERE Symbol != "-" AND ID IN (SELECT GeneID FROM TherapyTarget WHERE TherapyID IN (SELECT ID FROM Therapy WHERE DiseaseID IN (%s)));' % ",".join([str(i) for i in disease_id]))
  target_alias = pymysql_cursor('SELECT DISTINCT Alias FROM GeneAlias WHERE LENGTH (Alias) < 10 AND ID IN (SELECT GeneID FROM TherapyTarget WHERE TherapyID IN (SELECT ID FROM Therapy WHERE DiseaseID IN (%s)));' % ",".join([str(i) for i in disease_id]))
  if target:
    related_target.extend(target) if type(target) is list else related_target.append(target)
  if target_alias:
    related_target.extend(target_alias) if type(target_alias) is list else related_target.append(target_alias)
  
  related_target = set(related_target)
  print("%s genes are associated with this tumor type." % len(related_target))
  
  # Merge input annotation data
  variant = pd.DataFrame()
  if not somatic_anno.empty:
    somatic_anno.insert(0, 'Source', 'Somatic')
    variant = pd.concat([variant, somatic_anno], axis = 0)
  if not germline_anno.empty:
    germline_anno.insert(0, 'Source', 'Germline')
    # Filter by pathogenic
    germline_anno_p = germline_anno[germline_anno.CLNSIG.str.contains(r'Pathogenic|Likely_pathogenic', regex=True)]
    variant = pd.concat([variant, germline_anno_p], axis = 0)
  
  # Filtering 1: only save this tumor type related genes
  if variant.empty is False:
    variant = variant[variant['Gene.refGene'].isin(related_target)]
    ref_pop = 'AF_%s' % race
    cols = ['Source', 'Chr', 'Start', 'Ref', 'Alt', 'Gene.refGene', 'AAChange.refGene', 'ExonicFunc.refGene', ref_pop, 'CLNSIG', 'SIFT_pred'] #, 'Polyphen2_HVAR_pred', 'cosmic95_coding', 'cosmic95_noncoding']
    varf = variant[cols].reset_index(drop=True)
    # Only retain the first protein change
    varf = varf.copy().rename(columns={ref_pop: 'Population_AF', 'Gene.refGene':'Gene', 'AAChange.refGene':'AAChange', 'ExonicFunc.refGene':'Consequence'})
    
    if varf.empty is False:
      # varf = varf.drop(['AAChange'], axis=1).join(varf['AAChange'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('AAChange')).drop_duplicates().reset_index(drop = True)
      # Only save the first one
      varf = varf.drop(['AAChange'], axis=1).join(varf['AAChange'].str.split(',', expand=True)[0].rename('AAChange')).drop_duplicates().reset_index(drop = True)
      # varf['AAChange'] = varf['AAChange'].copy().str.split(',', expand=True)[0].map(parse_aachange)
    
    # varf.CLNSIG = varf.CLNSIG.map(lambda x: x.replace('_', ' ').title())
    # varf.columns = ['Source', 'Chr', 'Start', 'Ref', 'Alt', 'Gene', 'AAChange', 'Consequence', 'Population_AF', 'CLNSIG', 'SIFT_pred']#, 'Polyphen2_HVAR_pred']
    
    # Filtering 2
    varf['Population_AF'] = varf['Population_AF'].replace('.', np.nan).astype(float)
    varf = varf[(varf['SIFT_pred'] == 'D') | (varf['CLNSIG'].str.contains(r'Pathogenic|Likely_pathogenic', regex=True)) | (varf['Population_AF'] < 0.01)] # | (varf['Polyphen2_HVAR_pred'] != 'B')] # dbnsfp33a有Polyphen的信息，但是注释需要将近一个小时，是42c的十倍
    varf.insert(1, 'Alteration', '')
    varf.insert(2, 'Clinical_Significance', '')
    
    for index, row in varf.iterrows():
      # Alteration
      # chrom = "%s:%s:%s:%s" % (row.Chr, str(row.Start), row.Ref, row.Alt)
      # nm_exon_cc_pc = parse_aachange(row['AAChange'], row['Gene'])
      # if nm_exon_cc_pc:
      #   nm, exon, cc, pc = nm_exon_cc_pc[0]
      #   Alteration = ",".join([chrom, nm, exon, cc, pc]).strip(",")
      # else:
      #   Alteration = ",".join([chrom])
      Alteration = row['AAChange'].strip(row['Gene']).strip(':')
      varf.loc[index, 'Alteration'] = Alteration
      # Clinical Significance: Polyphen-2, SIFT, ClinVar
      clin_sig = []
      if row.CLNSIG != '.':
        clin_sig.append("%s (ClinVar)" % row.CLNSIG)
      if row.SIFT_pred != '.':
        clin_sig.append("%s (SIFT)" % {'D': 'Deleterious', 'T': 'Tolerated'}[row.SIFT_pred])
      # if row.Polyphen2_HVAR_pred != '.':
      #   clin_sig.append("%s (Polyphen2)" % {'B': 'Deleterious', 'T': 'Tolerated'}[row.Polyphen2_HVAR_pred])
      varf.loc[index, 'Clinical_Significance'] = "; ".join(clin_sig)
    
    biomarker_summary = varf[['Gene', 'Alteration', 'Source', 'Consequence', 'Population_AF', 'Clinical_Significance']]
  else:
    biomarker_summary = pd.DataFrame()
  
  
  #######################################################################################
  ### Based on the biomarker summary result, we will try to find the repurposing drugs
  alt_gene = biomarker_summary.Gene.drop_duplicates().to_list()
  print('%s genes remained after the filtering of ClinVar, SIFT and Population AF.' % len(alt_gene))
  ## Load pathway datasets
  hallmark = json.loads(open('./assets/h.all.v7.4.symbols.json').read())
  
  ## Find insection
  vargene_path_ins = []
  n = 0
  for key in hallmark.keys():
    # Whether there are drug targets within the pathway
    path = key
    path_gene = hallmark[key]
    path_has_target = list(set(path_gene)&set(related_target))
    if path_has_target:
      path_has_altg = list(set(set(path_gene))&set(alt_gene))
      if path_has_altg != []:
        print(path)
        for dt in path_has_target:
          for im in path_has_altg:
            joint_symbol = " ".join(sorted([dt, im]))
            # Find the combined score
            score = pymysql_cursor('SELECT DISTINCT Score FROM PPI WHERE JointSymbol = "%s";' % joint_symbol)
            if score:
              if type(score) == list:
                score = max(score)
              vargene_path_ins.append([im, dt, path, score])
  
  print("%s drug targeted abnormal genes are enriched in the same pathway." % len(vargene_path_ins))
  reuse = pd.DataFrame(vargene_path_ins, columns=['Gene', 'Associated_Gene', 'Pathway', 'PPI_Score'])
  reuse.insert(1, 'User_Alteration', "-")
  reuse.insert(2, 'Source', "Somatic")
  reuse.insert(3, 'Population_AF', "-")
  reuse.insert(4, 'GeneID', "")
  
  # Filtering based on cut-off
  reuse = reuse[reuse.PPI_Score > 0.98]
  print("%s of the gene pairs have a confidencial score > 0.98." % reuse.shape[0])
  reuse.sort_values(by = "PPI_Score", inplace=True, ascending=False)
  reuse = reuse.reset_index(drop=True)
  
  alt_num = reuse.shape[0]
  
  ### !!! IMPORTANT: Find the therapy range
  detail_ids = pymysql_cursor('SELECT DISTINCT DetailID FROM TherapyHasDetail WHERE TherapyID IN (SELECT ID FROM Therapy WHERE DiseaseID IN (%s));' % ",".join([str(i) for i in disease_id]))
  associated_therapy = []
  if reuse.empty == False:
    column_names = pymysql_cursor("SELECT column_name FROM information_schema.columns WHERE table_schema='OT' AND table_name='TherapyDetail';")
    reuse_names = reuse.columns.to_list()
    column_names.extend(reuse_names)
    for index, row in reuse.iterrows():
      # Find all variants corresponding to each gene
      reuse.loc[index, 'User_Alteration'] = str(biomarker_summary[biomarker_summary.Gene == row.Gene].Alteration.to_list())
      reuse.loc[index, 'Population_AF'] = str(biomarker_summary[biomarker_summary.Gene == row.Gene].Population_AF.to_list())
      ass_mg_id = get_metaid(row.Associated_Gene, "gene")[0]
      reuse.loc[index, 'GeneID'] = ass_mg_id
    
      ### Searching by row: Drug and therapy data
      if row.Source == 'Somatic':
        therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name REGEXP "Somatic" OR Name =  "Unknown") AND ID IN (%s) AND ID IN \
          (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN \
            (SELECT TherapyID FROM TherapyTarget WHERE GeneID = "%s"));' % (",".join([str(i) for i in detail_ids]), ass_mg_id))
      elif row.Source == 'Germline':
        therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name REGEXP "Germline" OR Name =  "Unknown") AND ID IN (%s) AND ID IN \
          (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN \
            (SELECT TherapyID FROM TherapyTarget WHERE GeneID = "%s"));' % (",".join([str(i) for i in detail_ids]), ass_mg_id))
      therapy_detail = cursor.fetchall()
      for item in therapy_detail:
        item = list(item)
        item.extend(row.to_list())
        associated_therapy.append(item)
    
    # Transfer to dataframe
    repurposing_therapy = pd.DataFrame(associated_therapy, columns=column_names).drop_duplicates()
    print("There are %s repurposing therapies." % repurposing_therapy.Therapy.drop_duplicates().shape[0])
  else:
    repurposing_therapy = pd.DataFrame()
    print("There are 0 repurposing therapies.")
  
  return((biomarker_summary, repurposing_therapy, alt_num))

