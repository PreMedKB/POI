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
  

  def anno_preprocess(anno, race, source, related_target):
    anno = anno[anno['Gene.refGene'].isin(related_target)].reset_index(drop=True)
    anno.insert(0, 'POID', ['%s.%s.indirect' % (source, str(i+1)) for i in anno.index.to_list()])
    anno.insert(1, 'Source', source)
    anno.insert(2, 'Location', str(anno['Chr'])+':'+str(anno['Start']))
    anno.insert(3, 'Ref>Alt', anno.Ref+'>'+anno.Alt)
    
    # ExonicFunc.refGene for splicing, intronic, UTR, etc., are all '.'.
    anno.insert(4, 'Function', anno[['ExonicFunc.refGene']])
    splicing_index = anno[anno['Func.refGene'].str.contains('splicing')].index
    anno.loc[splicing_index, 'Function'] = 'splicing'
    
    # Only keep the first AAChange
    # anno['AAChange.refGene'] = anno['AAChange.refGene'].str.split(',', expand=True)[0]
    # anno['GeneDetail.refGene'] = anno['GeneDetail.refGene'].str.split(';', expand=True)[0]
    anno.insert(5, 'Transcript', '.')
    anno.insert(6, 'Exon', '.')
    anno.insert(7, 'Protein_Change', '.')
    for index, row in anno.iterrows():
      nm_exon_cc_pc = parse_aachange(row['Gene.refGene'], row['AAChange.refGene'], row['GeneDetail.refGene'])
      if nm_exon_cc_pc != []:
        anno.loc[index, 'Transcript'] = nm_exon_cc_pc[0][0]
        anno.loc[index, 'Exon'] = nm_exon_cc_pc[0][1]
        anno.loc[index, 'Protein_Change'] = nm_exon_cc_pc[0][3]
    
    # Rename
    ref_pop = 'AF_%s' % race
    anno = anno.copy().rename(columns={'Gene.refGene':'Gene', 'ExonicFunc.refGene':'Exonic Function', 'CLNSIG':'Clinical_Significance', ref_pop: 'Population_AF'})
    
    # Add 'Pathogenic_Prediction'
    anno['Population_AF'] = anno['Population_AF'].replace('.', np.nan).astype(float)
    anno.insert(0, 'Pathogenic_Prediction', 'Tolerated')
    d_index = anno[(anno.SIFT_pred=="D")|(anno.LRT_pred=="D")|(anno.MutationTaster_pred=="A")|(anno.MutationTaster_pred=="D")|(anno.MutationAssessor_pred=="H")|(anno.MutationAssessor_pred=="M")|(anno.FATHMM_pred=="D")|(anno.PROVEAN_pred=="D")|(anno.MetaSVM_pred=="D")|(anno.MetaLR_pred=="D")].index # (Polyphen2_HDIV_pred=="D")|(Polyphen2_HDIV_pred=="P")|(Polyphen2_HVAR_pred=="D")|(Polyphen2_HVAR_pred=="P")
    anno.loc[d_index, 'Pathogenic_Prediction'] = 'Deleterious'
    
    # Select columns due to the potential difference between germline and somatic txt
    cols = ['POID', 'Gene', 'Source', 'Location', 'Ref>Alt', 'Transcript', 'Exon', 'Protein_Change', 'Function', 'Clinical_Significance', 'Pathogenic_Prediction', 'Population_AF']
    anno_res = anno[cols]
    
    return(anno_res)
  
  
  # Merge input annotation data
  variant = pd.DataFrame()
  if not somatic_anno.empty:
    variant = pd.concat([variant, anno_preprocess(somatic_anno, race, 'Somatic', related_target)], axis = 0)
  if not germline_anno.empty:
    # # Filter by pathogenic
    # germline_anno_p = germline_anno[germline_anno.CLNSIG.str.contains(r'Pathogenic|Likely_pathogenic', regex=True)]
    # variant = pd.concat([variant, germline_anno_p], axis = 0)
    variant = pd.concat([variant, anno_preprocess(germline_anno, race, 'Germline', related_target)], axis = 0)
    
  
  # Filter by Population_AF, Pathogenic Prediction, Clinical Significance
  index1 = variant[variant['Clinical_Significance'].str.contains(r'Pathogenic|Likely_pathogenic', regex=True)].index.to_list()
  index2 = variant[variant['Population_AF'] < 0.01].index.to_list()
  index3 = variant[variant['Pathogenic_Prediction'] == 'Deleterious'].index.to_list()
  index_123 = index1 + index2 + index3
  statistic = pd.DataFrame(pd.value_counts(index_123), columns=['Count'])
  index_final = statistic[statistic.Count > 1].index.to_list()
  varf = variant[variant.index.isin(index_final)]
  

  ###################
  # Try to find the repurposing drugs
  alt_gene = varf.Gene.drop_duplicates().to_list()
  print('%s genes remained after the filtering of ClinVar, Population_AF, SIFT, LRT, MutationTaster, MutationAssessor, FATHMM, PROVEAN, MetaSVM, MetaLR.' % len(alt_gene))
  # Load pathway datasets
  hallmark = json.loads(open('./assets/h.all.v7.4.symbols.json').read())
  # Find insection
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
              vargene_path_ins.append([im, dt, path.title(), score])
  
  print("%s drug targeted abnormal genes are enriched in the same pathway." % len(vargene_path_ins))
  reuse = pd.DataFrame(vargene_path_ins, columns=['Gene', 'Associated_Gene', 'Pathway', 'PPI_Score'])
  
  # Filtering based on cut-off
  reuse.sort_values(by = "PPI_Score", inplace=True, ascending=False)
  cutoff = 0.99
  reuse = reuse[reuse.PPI_Score >= cutoff]
  print("%s of the gene pairs have a confidencial score > %s." % (reuse.shape[0], cutoff))
  reuse = reuse.reset_index(drop=True)
  
  # Merge reuse with varf
  reuse = reuse.merge(varf)
  alt_num = reuse.shape[0]
  
  ### !!! IMPORTANT: Find the therapy range
  detail_ids = pymysql_cursor('SELECT DISTINCT DetailID FROM TherapyHasDetail WHERE TherapyID IN (SELECT ID FROM Therapy WHERE DiseaseID IN (%s));' % ",".join([str(i) for i in disease_id]))
  associated_therapy = []
  
  column_names = pymysql_cursor("SELECT column_name FROM information_schema.columns WHERE table_schema='OT' AND table_name='TherapyDetail';")
  column_names = ['POID'] + column_names + reuse.columns.to_list()
  
  if reuse.empty == False:
    for index, row in reuse.iterrows():
      # Find all variants corresponding to each gene
      ass_mg_id = get_metaid(row.Associated_Gene, "gene")[0]
      reuse.loc[index, 'GeneID'] = ass_mg_id
    
      ### Searching by row: Drug and therapy data
      if row.Source == 'Somatic':
        therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name REGEXP "Somatic") AND ID IN (%s) AND ID IN \
          (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN \
            (SELECT TherapyID FROM TherapyTarget WHERE GeneID = "%s"));' % (",".join([str(i) for i in detail_ids]), ass_mg_id))
      elif row.Source == 'Germline':
        therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name REGEXP "Germline") AND ID IN (%s) AND ID IN \
          (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN \
            (SELECT TherapyID FROM TherapyTarget WHERE GeneID = "%s"));' % (",".join([str(i) for i in detail_ids]), ass_mg_id))
      therapy_detail = cursor.fetchall()
      for item in therapy_detail:
        item = [row.POID] + list(item) + row.to_list()
        associated_therapy.append(item)
    
    # Transfer to dataframe
    repurposing_therapy = pd.DataFrame(associated_therapy, columns=column_names).drop_duplicates()
    print("There are %s repurposing therapies." % repurposing_therapy.Therapy.drop_duplicates().shape[0])
  else:
    repurposing_therapy = pd.DataFrame(columns=column_names)
    print("There are 0 repurposing therapies.")
  
  # Filter indirect biomarker
  indirect_biomarker = varf[varf.POID.isin(repurposing_therapy.POID.drop_duplicates())]
  
  return((repurposing_therapy, alt_num, indirect_biomarker))

