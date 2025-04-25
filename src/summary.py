##################################################
###   Merge and Report Multi-Sources Results   ###
##################################################
import pandas as pd
import numpy as np
import json, re, toml, os, pymysql
from src.functions import connect_database, load_logging_cfg, pymysql_cursor


## Aggregation and splicing
def concat_func(x):
  return pd.Series({'Level_Details':'|'.join(x['Level_Details'].unique()), 'Guidelines': '|'.join(x['Guidelines'].unique())})


def summary_report(merged_therapy, repurposing_therapy, chemotherapy, snpindel_biomarker, cnv_biomarker, fusion_biomarker, rna_biomarker, indirect_biomarker, multi_var, single_var, disease_id, db, cursor):
  
  drug_blacklist = open('./assets/drug_blacklist.txt').read().split('\n')
  
  if not merged_therapy.empty:
    merged_therapy = merged_therapy[merged_therapy.NormLevelID != 31]
  if not repurposing_therapy.empty:
    repurposing_therapy = repurposing_therapy[repurposing_therapy.NormLevelID != 31]
  
  snpindel_biomarker.fillna('', inplace=True)
  cnv_biomarker.fillna('', inplace=True)
  fusion_biomarker.fillna('', inplace=True)
  rna_biomarker.fillna('', inplace=True)
  indirect_biomarker.fillna('', inplace=True)
  multi_var.fillna('', inplace=True)
  single_var.fillna('', inplace=True)

  ################################## Part 1: directed drug recommendation
  print("Part 1: directed drug recommendation")
  direct_evidence = []
  for index, row in merged_therapy.iterrows():
    # Origin
    Origin = pymysql_cursor('SELECT Name FROM OriginDic WHERE ID = "%s";' % row.OriginID)
    # Find MetaDrugID AND MetaDiseaseID
    therapy = cursor.execute('SELECT * FROM Therapy WHERE ID IN (SELECT TherapyID FROM TherapyHasDetail WHERE DetailID = "%s");' % row.ID)
    therapy = cursor.fetchall()
    column_names2 = pymysql_cursor("SELECT column_name FROM information_schema.columns WHERE table_schema='OT' AND table_name='Therapy';")
    therapy = pd.DataFrame(therapy, columns=column_names2)
    therapy_disease_ids = therapy.DiseaseID.drop_duplicates().to_list()
    dis_range = ",".join([str(i) for i in therapy_disease_ids])
    # DrugID
    drugid = therapy.DrugID.drop_duplicates().to_list()
    # DrugSetID
    DrugSetID = therapy.DrugSetID.drop_duplicates().to_list()
    if DrugSetID != [None]:
      drugid = pymysql_cursor('SELECT DISTINCT MetaID FROM SubsetHasElement WHERE Class = "Drug" AND SubsetID IN (SELECT SubsetID FROM SetHasSubset WHERE SetID IN (%s));' % (",".join([str(i) for i in DrugSetID])))
      if type(drugid) != list:
        drugid = [drugid]
    
    # Drug display
    Drugs = pymysql_cursor('SELECT Therapy FROM MetaLite.Therapy WHERE ID = "%s";' % row.NormTherapyID)
    if Drugs not in drug_blacklist:
      drug_range = ",".join([str(i) for i in drugid])
      ### Guidelines
      Guidelines = {}
      NCCNDrug = cursor.execute('SELECT DISTINCT NCCNRecommendedUse, Category, Route FROM NCCNDrug.Main WHERE \
        DrugID IN (SELECT DISTINCT DomainID FROM NCCNDrug.Domain2Meta WHERE Class = "Drug" AND MetaID IN (%s)) AND \
          DiseaseID IN (SELECT DISTINCT DomainID FROM NCCNDrug.Domain2Meta WHERE Class = "Disease" AND MetaID IN (%s));' % (drug_range, dis_range))
      NCCNDrug = cursor.fetchall()
      multi_nccn = []
      if NCCNDrug:
        for item in NCCNDrug:
          multi_nccn.append("Recommended Use: %s. Category: %s. Route: %s." % (item[0], item[1], item[2]))
        Guidelines["NCCN"] = list(set(multi_nccn))
      
      Response = row.Response
      if Response in ['No Sensitivity', 'May Decrease Sensitivity', 'Reduced Sensitivity', 'Overcomes acquired resistance']:
        continue
      elif 'Sensitivity' in Response:
        Response = 'Sensitivity'
      elif 'Resistance' in Response:
        Response = 'Resistance'
      
      LevelDetails = {}
      ### Raw Level
      RawLevel = pymysql_cursor('SELECT Name FROM EvidenceLevelDic WHERE ID = "%s";' % row.RawLevelID)
      SourceDB = pymysql_cursor('SELECT Name FROM SourceDic WHERE ID = "%s";' % row.SourceID)
      if SourceDB in ['FUSCC_BZ', 'FUSCC_CZ']:
        continue
      ### Displayed Level
      map_rule = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E'}#, '31':'E#'}
      tmp_Level = map_rule[str(row.NormLevelID)]
      # Judging Level
      # 1. 疾病是提交疾病本身、子类、直系父类，等级不变：类似于将非小细胞肺癌药用于肺癌，无需降级
      if set(therapy_disease_ids).intersection(disease_id):
        Level = tmp_Level
      # 2. 疾病是同类组织的其他分支疾病，或者直接是其他组织的适应症，需要降级
      else:
        # Level = 'C*' if tmp_Level == 'A' else 'E*'
        Level = 'C' if tmp_Level == 'A' else 'E'
      
      ### Level Detail
      First = "Source Database: %s." % SourceDB
      raw_level = '-' if RawLevel == '-' else RawLevel.split('_')[1]
      EvidenceLevel = "Evidence Level: %s" % raw_level
      TumorType = "Tumor Type: %s" % row.TumorType
      Alteration = "Alteration: %s" % row.Alteration
      Therapy = "Therapy: %s" % row.Therapy
      Details = [First, EvidenceLevel, TumorType, Alteration, Therapy]
      others = pymysql_cursor('SELECT Annotation FROM ClinicalAnnotation WHERE DetailID = "%s";' % row.ID)
      if others:
        if type(others) != list:
          others = [others]
        for oth in others:
          if oth.startswith('Indication: '):
            if 'FDA-approved' in oth:
              Guidelines["FDA"] = [oth]
          else:
            Details.append(oth)
      URL = row.URL
      Details.append(URL)
      LevelDetails[SourceDB] = Details
      
      ### Biomarker Detail
      Small_Variant = []; CNV = []; Fusion = []; Expression = []
      if snpindel_biomarker.empty == False:
        Small_Variant = snpindel_biomarker[snpindel_biomarker.Gene == row.Gene].to_dict('records')
        Small_Variant = [str(x) for x in Small_Variant]
      if cnv_biomarker.empty == False:
        CNV = cnv_biomarker[cnv_biomarker.Gene == row.Gene].to_dict('records')
        CNV = [str(x) for x in CNV]
      if fusion_biomarker.empty == False:
        fusion_f = fusion_biomarker[(fusion_biomarker.gene1 == row.Gene) | (fusion_biomarker.gene2 == row.Gene)]
        fusion_f.insert(0, 'Gene_Pair', fusion_f.gene1 + '::' + fusion_f.gene2)
        Fusion = fusion_f.drop(columns=['gene1', 'gene2']).to_dict('records')
        Fusion = [str(x) for x in Fusion]
      if rna_biomarker.empty == False:
        Expression = rna_biomarker[rna_biomarker.Gene == row.Gene].to_dict('records')
        Expression = [str(x) for x in Expression]
        
      ## Add result
      direct_evidence.append([row.Gene, row.Variant, row.Source, Drugs, Response, Level, str(LevelDetails), str(Guidelines), str(Small_Variant), str(CNV), str(Fusion), str(Expression)])
  
  DE = pd.DataFrame(direct_evidence, columns=['Gene', 'Alteration', 'Source', 'Drugs', 'Response', 'Level', 'Level_Details', 'Guidelines', 'Small_Variant', 'CNV', 'Fusion', 'Expression'], dtype='object').drop_duplicates()
  
  ### Transform and filter
  DE = DE.sort_values('Level').reset_index(drop=True)
  DE_F = DE.groupby(['Gene', 'Alteration', 'Source', 'Drugs', 'Response', 'Level', 'Small_Variant', 'CNV', 'Fusion', 'Expression']).apply(concat_func).reset_index()
  for index, row in DE_F.iterrows():
    DE_F.loc[index, 'Level_Details'] = str(row.Level_Details.split('|'))
    DE_F.loc[index, 'Guidelines'] = str(row.Guidelines.split('|'))
  
  if DE_F.empty is False:
    DE_F = DE_F.sort_values('Level').reset_index(drop=True)
  
  ### JSON
  direct_json = DE_F.to_dict(orient='records')
  for item in direct_json:
    Level_Details = {}
    for detail in eval(item['Level_Details']):
      Level_Details.update(eval(detail))
    item['Level_Details'] = Level_Details
    
    Guidelines = {}
    for guide in eval(item['Guidelines']):
      Guidelines.update(eval(guide))
    item['Guidelines'] = Guidelines
    
    Small_Variant = []
    for i in range(0,len(eval(item['Small_Variant']))):
      ele = eval(item['Small_Variant'])[i]
      Small_Variant.append(eval(ele))
    item['Small_Variant'] = Small_Variant
    
    CNV = []
    for i in range(0,len(eval(item['CNV']))):
      ele = eval(item['CNV'])[i]
      CNV.append(eval(ele))
    item['CNV'] = CNV
    
    Fusion = []
    for i in range(0,len(eval(item['Fusion']))):
      ele = eval(item['Fusion'])[i]
      Fusion.append(eval(ele))
    item['Fusion'] = Fusion

    Expression = []
    for i in range(0,len(eval(item['Expression']))):
      ele = eval(item['Expression'])[i]
      Expression.append(eval(ele))
    item['Expression'] = Expression
  
  
  ################################## Part 2: indirected drug recommendation
  print("Part 2: indirected drug recommendation")
  indirect_evidence = []
  for index, row in repurposing_therapy.iterrows():
    # Origin
    Origin = pymysql_cursor('SELECT Name FROM OriginDic WHERE ID = "%s";' % row.OriginID)
    # 找到MetaDrugID和MetaDiseaseID
    therapy = cursor.execute('SELECT * FROM Therapy WHERE ID IN (SELECT TherapyID FROM TherapyHasDetail WHERE DetailID = "%s");' % row.ID)
    therapy = cursor.fetchall()
    column_names2 = pymysql_cursor("SELECT column_name FROM information_schema.columns WHERE table_schema='OT' AND table_name='Therapy';")
    therapy = pd.DataFrame(therapy, columns=column_names2)
    therapy_disease_ids = therapy.DiseaseID.drop_duplicates().to_list()
    dis_range = ",".join([str(i) for i in therapy_disease_ids])
    # DrugID
    drugid = therapy.DrugID.drop_duplicates().to_list()
    # DrugSetID
    DrugSetID = therapy.DrugSetID.drop_duplicates().to_list()
    if DrugSetID != [None]:
      drugid = pymysql_cursor('SELECT DISTINCT MetaID FROM SubsetHasElement WHERE Class = "Drug" AND SubsetID IN \
        (SELECT SubsetID FROM SetHasSubset WHERE SetID IN (%s));' % (",".join([str(i) for i in DrugSetID])))
      if type(drugid) != list:
        drugid = [drugid]
    
    # Drugs = pymysql_cursor('SELECT DISTINCT NormDrug FROM DrugDisplay WHERE Drug = "%s";' % row.Therapy)
    # if Drugs == None and DrugSetID != [None]:
    #     Drugs = pymysql_cursor('SELECT DISTINCT NormDrug FROM DrugDisplay WHERE DrugSetID IN (%s);' % (",".join([str(i) for i in DrugSetID])))
    #     if type(Drugs) == list:
    #       Drugs = "%".join(Drugs)
    # if Drugs == None:
    #     Drugs = row.Therapy
    Drugs = pymysql_cursor('SELECT Therapy FROM MetaLite.Therapy WHERE ID = "%s";' % row.NormTherapyID)
    if Drugs not in drug_blacklist:
      drug_range = ",".join([str(i) for i in drugid])
      ### Guidelines
      Guidelines = {}
      NCCNDrug = cursor.execute('SELECT DISTINCT NCCNRecommendedUse, Category, Route FROM NCCNDrug.Main WHERE \
        DrugID IN (SELECT DISTINCT DomainID FROM NCCNDrug.Domain2Meta WHERE Class = "Drug" AND MetaID IN (%s)) AND \
          DiseaseID IN (SELECT DISTINCT DomainID FROM NCCNDrug.Domain2Meta WHERE Class = "Disease" AND MetaID IN (%s));' % (drug_range, dis_range))
      NCCNDrug = cursor.fetchall()
      multi_nccn = []
      if NCCNDrug:
        for item in NCCNDrug:
          multi_nccn.append("Recommended Use: %s. Category: %s. Route: %s." % (item[0], item[1], item[2]))
        Guidelines["NCCN"] = list(set(multi_nccn))
      
      Response = row.Response
      if Response in ['No Sensitivity', 'May Decrease Sensitivity', 'Reduced Sensitivity', 'Overcomes acquired resistance']:
        continue
      elif 'Sensitivity' in Response:
        Response = 'Sensitivity'
      elif 'Resistance' in Response:
        Response = 'Resistance'
      
      LevelDetails = {}
      ### Level
      Level = 'E'
      RawLevel = pymysql_cursor('SELECT Name FROM EvidenceLevelDic WHERE ID = "%s";' % row.RawLevelID)
      SourceDB = pymysql_cursor('SELECT Name FROM SourceDic WHERE ID = "%s";' % row.SourceID)
      if SourceDB in ['FUSCC_BZ', 'FUSCC_CZ']:
        continue
      ### Level Detail
      First = "Source Database: %s." % SourceDB
      raw_level = '-' if RawLevel == '-' else RawLevel.split('_')[1]
      EvidenceLevel = "Evidence Level of the Associated Gene: %s." % raw_level
      # Phenotype = "Phenotype: %s." % row.TumorType
      Alteration = "Alteration: %s." % row.Alteration
      # Therapy = "Therapy: %s." % row.Therapy
      Details = [First, EvidenceLevel, Alteration]
      others = pymysql_cursor('SELECT Annotation FROM ClinicalAnnotation WHERE DetailID = "%s";' % row.ID)
      if others:
        if type(others) != list:
          others = [others]
        for oth in others:
          if oth.startswith('Indication: '):
            if 'FDA-approved' in oth:
              Guidelines["FDA"] = oth
          else:
            Details.append(oth)
      URL = row.URL
      Details.append(URL)
      LevelDetails[SourceDB] = Details

      ### Biomarker Detail
      Small_Variant = []
      if indirect_biomarker.empty == False:
        Small_Variant = indirect_biomarker[indirect_biomarker.Gene == row.Gene].to_dict('records')
        Small_Variant = [str(x) for x in Small_Variant]
      
      ## Add result
      indirect_evidence.append([row.Gene, row.Source, row.Associated_Gene, row.Pathway, row.PPI_Score, Drugs, Response, Level, str(LevelDetails), str(Guidelines), str(Small_Variant)])
  
  ### Transform and filter
  INDE = pd.DataFrame(indirect_evidence, columns=['Gene', 'Source', 'Associated_Gene', 'Pathway', 'PPI_Score', 'Drugs', 'Response', 'Level', 'Level_Details', 'Guidelines', 'Small_Variant'], dtype='object').drop_duplicates()
  INDE = INDE.sort_values('Level').reset_index(drop=True)
  INDE_F = INDE.groupby(['Gene', 'Source', 'Associated_Gene', 'Pathway', 'PPI_Score', 'Drugs', 'Response', 'Level', 'Small_Variant']).apply(concat_func).reset_index()
  for index, row in INDE_F.iterrows():
    INDE_F.loc[index, 'Level_Details'] = str(row.Level_Details.split('|'))
    INDE_F.loc[index, 'Guidelines'] = str(row.Guidelines.split('|'))
  
  if INDE_F.empty is False:
    INDE_F = INDE_F.sort_values('Guidelines').reset_index(drop=True)
  
  ### JSON
  indirect_json = INDE_F.to_dict(orient='records')
  for item in indirect_json:
    Level_Details = {}
    for detail in eval(item['Level_Details']):
      Level_Details.update(eval(detail))
    item['Level_Details'] = Level_Details
    
    Guidelines = {}
    for guide in eval(item['Guidelines']):
      Guidelines.update(eval(guide))
    item['Guidelines'] = Guidelines
    
    Small_Variant = []
    for i in range(0,len(eval(item['Small_Variant']))):
      ele = eval(item['Small_Variant'])[i]
      Small_Variant.append(eval(ele))
    item['Small_Variant'] = Small_Variant
  

  ################################## Part 3: drug_response
  print("Part 3: drug_response")
  drug_response = []
  for index, row in chemotherapy.iterrows():
    SourceDB = 'PharmGKB'
    Origin = 'Germline'
    Gene = row.Gene
    
    # Drugs = pymysql_cursor('SELECT DISTINCT NormDrug FROM DrugDisplay WHERE Drug = "%s";' % row.Drug)
    # if Drugs == None:
    #     Drugs = row.Drug
    Drugs = row.Drug.title()#pymysql_cursor('SELECT Therapy FROM MetaLite.Therapy WHERE ID = "%s";' % row.NormTherapyID)
    
    ### Guidelines
    Guidelines = {}
    Guidelines[row.Source] = [row.Recommendation]
    Diplotype = row.Variant + ' ' + row.Diplotype; Diplotype = Diplotype.strip()
    Category = row.PhenotypeCategory
    Response = row.PAnnoPhenotype
    if type(Response) != str:
      Response = row.Phenotype
    LevelDetails = {}
    ### Raw Level
    if type(row.EvidenceLevel) is str:
      RawLevel = row.EvidenceLevel
      ### Displayed Level
      map_rule = {'1A':'A', '1B':'A', '2A':'B', '2B':'B'}
      Level = map_rule[RawLevel]
    else:
      RawLevel = '-'
      Level = 'A'
    
    ### Level Detail
    First = "Source Database: %s." % SourceDB
    raw_level = '-' if RawLevel == '-' else RawLevel
    EvidenceLevel = "Evidence Level: %s." % raw_level
    # Phenotype = "Phenotype: %s." % row.TumorType
    Alteration = "Alteration: %s." % row.Variant
    # Therapy = "Therapy: %s." % row.Therapy
    Details = [First, EvidenceLevel, Alteration]
    URL = "https://www.pharmgkb.org/clinicalAnnotation/%s" % row.CAID
    Details.append(URL)
    LevelDetails[SourceDB] = Details
    
    # Biomarker Detail
    Multi_Variant = multi_var[multi_var.Gene == Gene].to_dict('records')
    Multi_Variant = [str(x) for x in Multi_Variant]
    Single_Variant = single_var[single_var.Gene == Gene].to_dict('records')
    Single_Variant = [str(x) for x in Single_Variant]
    
    ## Add result
    drug_response.append([Gene, Diplotype, 'Germline', Drugs, Category, Response, Level, str(LevelDetails), str(Guidelines), str(Multi_Variant), str(Single_Variant)])
  
  DR = pd.DataFrame(drug_response, columns=['Gene', 'Diplotype', 'Source', 'Drugs', 'Category', 'Response', 'Level', 'Level_Details', 'Guidelines', 'Multi_Variant', 'Single_Variant'], dtype='object').drop_duplicates()
  DR = DR[DR.Gene != '-']
  ## 如果将Level Details和Guidelines合并
  # 1. 先排序
  DR = DR.sort_values('Level').reset_index(drop=True)
  # 2. 找到重复行，合并Level Details
  DR_F = DR.groupby(['Gene', 'Diplotype', 'Source', 'Drugs', 'Category', 'Response', 'Level', 'Multi_Variant', 'Single_Variant']).apply(concat_func).reset_index()
  for index, row in DR_F.iterrows():
    DR_F.loc[index, 'Level_Details'] = str(row.Level_Details.split('|'))
    DR_F.loc[index, 'Guidelines'] = str(row.Guidelines.split('|'))
  
  if DR_F.empty is False:
    DR_F = DR_F.sort_values('Level').reset_index(drop=True)
  
  ## JSON
  response_json = DR_F.to_dict(orient='records')
  for item in response_json:
    Level_Details = {}
    for detail in eval(item['Level_Details']):
      Level_Details.update(eval(detail))
    item['Level_Details'] = Level_Details
    
    Guidelines = {}
    for guide in eval(item['Guidelines']):
      Guidelines.update(eval(guide))
    item['Guidelines'] = Guidelines
    
    Multi_Variant = []
    for i in range(0,len(eval(item['Multi_Variant']))):
      ele = eval(item['Multi_Variant'])[i]
      Multi_Variant.append(eval(ele))
    item['Multi_Variant'] = Multi_Variant
    
    Single_Variant = []
    for i in range(0, len(eval(item['Single_Variant']))):
      ele = eval(item['Single_Variant'])[i]
      Single_Variant.append(eval(ele))
    item['Single_Variant'] = Single_Variant
  
  
  ################################## Part 4: Generating a single file
  print("Part 4: Detail Report")
  detail = {
    "Direct_Evidence": direct_json,
    "Indirect_Evidence": indirect_json,
    "Drug_Response": response_json,
  }
  
  ################################## Part 5: Summary section
  print("Part 5: Summary")
  DE.insert(0, 'Class', 'direct')
  INDE.insert(0, 'Class', 'indirect')
  DR.insert(0, 'Class', 'response')
  used_col = ['Drugs', 'Response', 'Level', 'Class']
  
  DE.Alteration.drop_duplicates().shape[0] + INDE.Gene.drop_duplicates().shape[0]
  # subDR = DR[DR.Response.str.contains('Toxicity, increased|Toxicity, moderate')]
  # merge_df = pd.concat([subDR[used_col], INDE[used_col], DE[used_col]], axis=0).drop_duplicates()
  merge_df = pd.concat([INDE[used_col], DE[used_col]], axis=0).drop_duplicates()
  
  # Keep the highest one
  merge_df = merge_df.sort_values(by = 'Level').drop_duplicates('Drugs', keep='first').reset_index(drop=True)
  # resistance / ad
  res = merge_df[(merge_df.Response.str.contains('.*Resistance.*|No Sensitivity', regex=True, case=False))].Drugs.drop_duplicates().to_list()
  # ad = merge_df[merge_df.Response.str.contains('Toxicity, increased|Toxicity, moderate')].Drugs.drop_duplicates().to_list()
  therapy_summary = {}; therapy_num = 0
  for level in ['A', 'B', 'C', 'D', 'E']:#merge_df.Level.drop_duplicates().to_list():
    sensitivity = []; resistance = []
    df = merge_df[merge_df.Level == level]
    for index, row in df.iterrows():
      # Cleaning: 临时解决方法
      # DRUG = re.sub('\-|\s', '_', row.Drugs.replace(" + ", "+"))
      DRUG = row.Drugs
      if row.Drugs in res:
        resistance.append(DRUG)
      else:
        sensitivity.append(DRUG)
    
    sensitivity.sort()
    resistance.sort()
    num = len(sensitivity)+len(resistance)
    therapy_num = therapy_num + num
    therapy_summary["Level_%s" % level] = {"sensitivity": sensitivity, "resistance": resistance, "num": num}
  
  return (therapy_summary, detail, therapy_num)

