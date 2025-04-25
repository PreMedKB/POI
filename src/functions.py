#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Basic functions.
"""

# import modules
import sys, toml, pymysql, re
import pandas as pd
import logging.config
from warnings import filterwarnings


### Config
def db_config(cfile):
  global DEFAULT, JSON, TSV, SCHEMA
  cf = toml.load(cfile, _dict=dict)
  DEFAULT = cf['DEFAULT']
  JSON = cf['JSON']
  TSV = cf['TSV']
  SCHEMA = cf['SCHEMA']


### Connect database
def connect_database(DEFAULT):
  global db, cursor 
  # Connect mysql
  try:
    db = pymysql.connect(host = DEFAULT['host'], user = DEFAULT['user'], port = DEFAULT['port'], passwd = DEFAULT['passwd'], charset = DEFAULT['charset'], local_infile = DEFAULT['local_infile'])
    # Create a cursor object that can execute SQL statements
    cursor = db.cursor()
    cursor.execute("USE %s;" % DEFAULT['database_name'])
  except pymysql.Error as e:
    db.rollback()
    logger.error(e)
    # Define "3" as the db connecting error
    sys.exit(3)
  except pymysql.Warning as w:
    db.rollback()
    logger.warning(w)
  
  return(db, cursor)


### Execute SQL command; Capture errors and warnings
def pymysql_cursor(sql):
  # try... except... can only capture errors
  filterwarnings("error", category = pymysql.Warning)
  try:
    # cursor.execute(sql) will automatically return the number of affected rows
    effected_rows = cursor.execute(sql)
    db.commit()
  except pymysql.Error as e:
    db.rollback()
    logger.error(e)
  except pymysql.Warning as w:
    db.rollback()
    logger.warning(w)
  else:
    if re.match('^SELECT', sql):
      # cursor.fetchone() return tuple; cursor.execute(sql) returns the number of rows affected 
      # one query result: ((1,)); multiple results: ((1,),(2,))
      result = cursor.fetchall()
      if len(result) == 1:
        single_result = result[0][0]
        return single_result
      elif len(result) > 1:
        multiple_results = []
        for ele in result:
          multiple_results.append(ele[0])
        return multiple_results
    elif re.match('^LOAD DATA LOCAL INFILE', sql):
      return cursor.rowcount



### Logger
def load_logging_cfg(DEFAULT):
  # Define three log output formats
  standard_format = '%(asctime)-15s  [%(threadName)s:%(thread)d, task_id:%(name)s, %(filename)s:%(lineno)d]' \
          '\n[%(levelname)s]  %(message)s\n' 
  simple_format = '[%(levelname)s]  %(asctime)-15s\n%(message)s\n'
  #id_simple_format = '[%(levelname)s][%(asctime)s]%(message)s'
  # Full path of log file
  logfile = DEFAULT['logfile']
  # Configuration of the logger
  LOGGING_DIC = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
      'standard': {
        'format': standard_format
      },
      'simple': {
        'format': simple_format
      },
    },
    'filters': {},
    'handlers': {
      # Print logs to terminal
      'console': {
        'level': 'DEBUG',
        'class': 'logging.StreamHandler',
        'formatter': 'simple',
      },
      # Print logs into file, collect logs above 'INFO' 
      'default': {
        'level': 'DEBUG',
        'class': 'logging.handlers.RotatingFileHandler', 
        'formatter': 'standard',
        'filename': logfile, 
        'maxBytes': 1024*1024*5, 
        'backupCount': 5,
        'encoding': 'utf-8',
      },
    },
    'loggers': {
      # Using logging.getLogger(__name__) to get the configuration
      '': {
        'handlers': ['default', 'console'],
        'level': 'DEBUG',
        'propagate': True,
      },
    },
  }  
  logging.config.dictConfig(LOGGING_DIC)
  global logger
  logger = logging.getLogger(__name__)
  return(logger)



## Function used to parse the AAVariant.
def parse_aachange(gene, aachange, genedetail):
  transcript = pd.read_csv('./assets/transcript.txt', sep = '\t')
  nm_exon_cc_pc = []
  if aachange != '.' and aachange != [None] and aachange != 'UNKNOWN':
    if type(aachange) is list: # VCF
      AAChange = aachange
    else: # TXT
      AAChange = aachange.split(',')
    for x in range(0, len(AAChange)):
      i = AAChange[x]
      nm = re.findall(r'\w+:(\w+):\w+:c.[\w\->]+:p.\w+', i)[0] if re.findall(r'\w+:(\w+):\w+:c.[\w\->]+:p.\w+', i) else ''
      exon = re.findall(r'\w+:\w+:(exon\d+):c.[\w\->]+:p.\w+', i)[0] if re.findall(r'\w+:\w+:(exon\d+):c.[\w\->]+:p.\w+', i) else ''
      cc = re.findall(r':(c.[\w\->]+)', i)[0] if re.findall(r':(c.[\w\->]+)', i) else ''
      pc = re.findall(r':(p.\w+)', i)[0] if re.findall(r':(p.\w+)', i) else ''
      # Filter by transcript
      if gene in transcript.gene_name.to_list():
        if nm == transcript[transcript.gene_name == gene]['nm_name'].to_list()[0]:
          nm_exon_cc_pc.append([nm, exon, cc, pc])
      elif x == 0:
        nm_exon_cc_pc.append([nm, exon, cc, pc])
  elif genedetail != '.' and aachange == '.':
    GeneDetail = genedetail.split(',')
    for x in range(0, len(GeneDetail)):
      i = GeneDetail[x]
      nm = re.findall(r'(\w+):\w+:c.[\w\->]+', i)[0] if re.findall(r'(\w+):\w+:c.[\w\->]+', i) else ''
      exon = re.findall(r'\w+:(exon\d+):c.[\w\->]+', i)[0] if re.findall(r'\w+:(exon\d+):c.[\w\->]+', i) else ''
      cc = re.findall(r':(c.[\w\->]+)', i)[0] if re.findall(r':(c.[\w\->]+)', i) else ''
      pc = ''
      # Filter by transcript
      if gene in transcript.gene_name.to_list():
        if nm == transcript[transcript.gene_name == gene]['nm_name'].to_list()[0]:
          nm_exon_cc_pc.append([nm, exon, cc, pc])
      elif x == 0:
        nm_exon_cc_pc.append([nm, exon, cc, pc])
  
  return (nm_exon_cc_pc)




############################################
###   Functions Used to Search Therapy   ###
############################################
def get_metaid(name, typ):
  if typ == "gene":
    meta_id = pymysql_cursor('SELECT ID FROM Gene WHERE Symbol = "%s";' % name)
    if meta_id is None:
      meta_id =  pymysql_cursor('SELECT GeneID FROM GeneAlias WHERE Alias = "%s";' % name)
    # Collect the meta_ids
    if meta_id is None:
      meta_id = pymysql_cursor('SELECT ID FROM Gene WHERE Symbol = "-";')
    if type(meta_id) != list:
      meta_id = [meta_id]
  
  return (meta_id)


def search_snpindel(gene, var, source, assembly):
  mgene_id = get_metaid(gene, "gene")
  mgid_range = ",".join([str(gid) for gid in mgene_id])
  mvar_ids = []
  assembly_ids = pymysql_cursor('SELECT ID FROM AssemblyDic WHERE Name IN ("%s", "na");' % assembly)
  assembly_range = ",".join([str(assembly_id) for assembly_id in assembly_ids])
  # Get MetaVariantID
  if type(var) is str:
    if 'deletion' in var:
      mvids = pymysql_cursor('SELECT ID FROM Variant WHERE Symbol = "{0}" AND Name LIKE "%deletion%" AND AssemblyID IN ({1});'.format(gene, assembly_range))
    else:
      mvids = pymysql_cursor('SELECT ID FROM Variant WHERE Symbol = "{0}" AND Name = "{1}" AND AssemblyID IN ({2});'.format(gene, var, assembly_range))
    if mvids:
      if type(mvids) != list:
        mvar_ids.append(mvids)
      else:
        mvar_ids.extend(mvids)
  else:
    chr, pos, ref, alt, nm, exon, cc, pc, fuc, clinsig = var[0], var[1], var[2], var[3], var[4], var[5], var[6], var[7], var[8], var[9]
    ### 1. Search by AAchange
    if pc != '':
      if pc.endswith('*fs') or pc.endswith('fs*'):
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND ID IN (SELECT VariantID FROM VariantHasTranscript WHERE TranscriptID IN (SELECT ID FROM Transcript WHERE HGVSpShort LIKE "{1}%")) AND AssemblyID IN ({2});'.format(mgid_range, pc, assembly_range))
      else:
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN (%s) AND ID IN (SELECT VariantID FROM VariantHasTranscript WHERE TranscriptID IN (SELECT ID FROM Transcript WHERE HGVSpShort = "%s")) AND AssemblyID IN (%s);' % (mgid_range, pc, assembly_range))
      if mvids:
        if type(mvids) != list:
          mvar_ids.append(mvids)
        else:
          mvar_ids.extend(mvids)
    ### 2. Search by Position
    chr = re.sub('chr|Chr|CHR', '', str(chr))
    mvids = pymysql_cursor('SELECT ID FROM Variant WHERE Chr = "%s" AND Pos = "%s" AND RefBase = "%s" AND AltBase IN (%s) AND AssemblyID IN (%s);' % (chr, pos, ref, ','.join(['"%s"' % a for a in alt]), assembly_range))
    if mvids:
      if type(mvids) != list:
        mvar_ids.append(mvids)
      else:
        mvar_ids.extend(mvids)
    ### 3. Truncating Mutations
    if fuc in ['stopgain', 'frameshift insertion', 'frameshift deletion', 'frameshift block substitution']:
      mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "{1}" AND AssemblyID IN ({2});'.format(mgid_range, '%trunc%', assembly_range))
      if mvids:
        if type(mvids) != list:
          mvar_ids.append(mvids)
        else:
          mvar_ids.extend(mvids)
    ### 4. Oncogenic Mutations
    if re.findall(r'Pathogenic|Likely_pathogenic', clinsig):
      mvids = pymysql_cursor('SELECT ID, Name FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "{1}" AND AssemblyID IN ({2});'.format(mgid_range, '%oncogenic%', assembly_range))
      if mvids:
        if type(mvids) != list:
          mvar_ids.append(mvids)
        else:
          mvar_ids.extend(mvids)
    ### 5. nonframeshift deletion (p.xxx_xxxdel) #& nonframeshift insertion (p.xxxdelinsxxx)
    if fuc in ['nonframeshift deletion', 'nonframeshift insertion']:
      if re.findall('p.\d+_\d+del', pc):
        start = re.findall('p.(\d+)_(\d+)del', pc)[0][0]
        end = re.findall('p.(\d+)_(\d+)del', pc)[0][1]
        mvids = []
        # \_ will be ineffective
        del_df = cursor.execute('SELECT ID, Name FROM Variant WHERE GeneID IN ({0}) AND AssemblyID IN ({1}) AND Name LIKE "{2}" and Name NOT LIKE "{3}";'.format(mgid_range, assembly_range, "%p.%_%del%", "%delins%"))
        del_df = cursor.fetchall()
        del_df = pd.DataFrame(del_df, columns = ['ID', 'Name'])#; print(del_df)
        for index, row in del_df.iterrows():
          # \_ will be ineffective. Recheck.
          match = re.findall('p.(?:[a-zA-Z]+)?(\d+)_(?:[a-zA-Z]+)?(\d+)del', row.Name)
          if match:
            start_def = match[0][0]
            end_def = match[0][1]
            # Determine the intersection
            check_list = [[int(start), int(end)], [int(start_def), int(end_def)]]
            check_list_sort = sorted(check_list, key=lambda l: l[0])
            for i in range(0, len(check_list_sort)-1):
              if check_list_sort[i+1][0] <= check_list_sort[i][1]:
                print('Nonframeshift INDELs: %s overlaps with %s!'% (str(check_list_sort[i]), str(check_list_sort[i+1])))
                mvids.append(row.ID)
      elif 'insertion' in fuc:
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN (%s) AND Name = "insertion" AND AssemblyID IN (%s);' % (mgid_range, assembly_range))
      elif 'deletion' in fuc:
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN (%s) AND Name = "deletion" AND AssemblyID IN (%s);' % (mgid_range, assembly_range))
      if mvids:
        if type(mvids) != list:
          mvar_ids.append(mvids)
        else:
          mvar_ids.extend(mvids)
    ### 6. Exon Variant: e.g., Exon 18 in-frame insertions, Exon 18 missense mutations
    exon = exon.replace('exon', '')
    mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "%{1}%" AND Name LIKE "%{2}%" AND Name LIKE "%{3}%" AND AssemblyID IN ({4});'.format(mgid_range, exon, fuc, 'exon', assembly_range))
    if mvids is None:
      if fuc == 'nonframeshift deletion': # NOT LIKE deletion
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "%{1}%" AND Name LIKE "%{2}%" AND Name NOT LIKE "%{3}%" AND Name LIKE "%{4}%" AND AssemblyID IN ({5});'.format(mgid_range, exon, 'deletion', 'frame', 'exon', assembly_range))
      elif fuc == 'nonframeshift insertion': # NOT LIKE insertion
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "%{1}%" AND Name LIKE "%{2}%" AND Name NOT LIKE "%{3}%" AND Name LIKE "%{4}%" AND AssemblyID IN ({5});'.format(mgid_range, exon, 'insertion', 'frame', 'exon', assembly_range))
      elif fuc == 'frameshift deletion':
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "%{1}%" AND Name LIKE "%{2}%" AND Name LIKE "%{3}%" AND Name LIKE "%{4}%" AND AssemblyID IN ({5});'.format(mgid_range, exon, 'deletion', 'frame', 'exon', assembly_range))
      elif fuc == 'frameshift insertion':
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "%{1}%" AND Name LIKE "%{2}%" AND Name LIKE "%{3}%" AND Name LIKE "%{4}%" AND AssemblyID IN ({5});'.format(mgid_range, exon, 'insertion', 'frame', 'exon', assembly_range))
      elif fuc == 'splicing':
        mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "%{1}%" AND Name LIKE "%{2}%" AND Name LIKE "%{3}%" AND AssemblyID IN ({4});'.format(mgid_range, exon, 'splic', 'exon', assembly_range))
    if mvids is None:
      mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN ({0}) AND Name LIKE "%{1}%" AND Name LIKE "%{2}%" AND Name LIKE "%{3}%" AND AssemblyID IN ({4});'.format(mgid_range, exon, 'mutation', 'exon', assembly_range))
    if mvids:
      if type(mvids) != list:
        mvar_ids.append(mvids)
      else:
        mvar_ids.extend(mvids)
    ### 7. Exonic function without considering specific exons
    mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN (%s) AND Name = "%s" AND AssemblyID IN (%s);' % (mgid_range, fuc, assembly_range))
    if mvids:
      if type(mvids) != list:
        mvar_ids.append(mvids)
      else:
        mvar_ids.extend(mvids)
    ### 8. Missense Mutation: https://annovar.openbioinformatics.org/en/latest/user-guide/gene/
    if fuc not in ['unknown', 'synonymous SNV']:
      mvids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN (%s) AND Name = "%s" AND AssemblyID IN (%s);' % (mgid_range, 'Missense Mutation', assembly_range))
    if mvids:
      if type(mvids) != list:
        mvar_ids.append(mvids)
      else:
        mvar_ids.extend(mvids)
    ### 9. Transcript and nucleotide change (nm, cc)
    # mvar_ids = pymysql_cursor('SELECT ID FROM Variant WHERE GeneID IN (%s) AND ID IN (SELECT VariantID FROM VariantHasTranscript WHERE TranscriptID IN (SELECT ID FROM Transcript WHERE NucleotideHGVS LIKE "%s%s%s" AND NucleotideHGVS LIKE "%s%s%s")) AND AssemblyID IN (%s);' % (mgid_range, '%', nm, '%', '%', cc, '%', assembly_range))
  
  ## Generate var_range for searching Therapy
  mvar_ids = list(set(mvar_ids))
  var_range = ",".join([str(var_id) for var_id in mvar_ids])
  if var_range != "":
    ## Somatic and Germline
    if re.search('Somatic', source, flags=re.IGNORECASE):
      therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name LIKE "%{0}%") AND ID IN (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN (SELECT ID FROM Therapy WHERE VariantID IN ({1})));'.format('Somatic', var_range))
    elif re.search('Germline', source, flags=re.IGNORECASE):
      therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name LIKE "%{0}%") AND ID IN (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN (SELECT ID FROM Therapy WHERE VariantID IN ({1})));'.format('Germline', var_range))
    therapy_detail = cursor.fetchall()
    # Change the format
    out = [list(item) for item in therapy_detail]
    column_names = ['ID', 'SourceID', 'EntryID', 'Support', 'Alteration', 'TumorType', 'Therapy', 'Allele', 'Response', 'URL', 'OriginID', 'RawLevelID', 'NormLevelID', 'NormTherapyID']
    out_df = pd.DataFrame(out, columns=column_names)
    Gene = gene
    if type(var) is str:
      Alteration = var
    else:
      chrom = "%s:%s" % (chr, pos)
      refalt = "%s>%s" % (ref, alt)
      if pc != []:
        Alteration = "chr%s,%s,%s,exon %s,%s" % (chrom, refalt, nm, exon, pc)
      elif nm != [] and exon != []:
        Alteration = "chr%s,%s,%s,exon %s" % (chrom, refalt, nm, exon)
      elif nm != []:
        Alteration = "chr%s,%s,%s" % (chrom, refalt, nm)
      elif exon != []:
        Alteration = "chr%s,%s,exon %s" % (chrom, refalt, exon)
      else:
        Alteration = "chr%s,%s" % (chrom, refalt)
    out_df.insert(0, 'Gene', Gene)
    out_df.insert(1, 'Variant', Alteration)
    out_df.insert(2, 'Source', source)
    # Judge whether AAchange (pc) of the input is excluded.
    out_df_index = out_df.index.to_list()
    for index, row in out_df.iterrows():
      exclude = re.findall(".*\(excluding (.*)\).*", row.Alteration)
      if exclude:
        exclude_vars = re.split(', | and ', exclude[0])
        for ex_var in exclude_vars:
          if re.match(r'.*[0-9]$', ex_var):
            AA_abbr = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
            for aa in AA_abbr:
              exclude_vars.append('%s%s' % (ex_var, aa))
        for ex_var in exclude_vars:
          ex_var = 'p.%s' % ex_var
          if ex_var == pc:
            out_df_index.remove(index)
            break
    out_df = out_df.iloc[out_df_index,]
    # Detected variants
    detected_var = (var_range, Gene, Alteration, source)
  else:
    out_df = pd.DataFrame()
    detected_var = ()
  
  return(out_df, detected_var)



def search_others(gene, var, var_type, assembly):
  assembly_ids = pymysql_cursor('SELECT ID FROM AssemblyDic WHERE Name IN ("%s", "na");' % assembly)
  assembly_range = ",".join([str(assembly_id) for assembly_id in assembly_ids])
  # Get MetaVariantID (Single variants)
  # Name in Variant may be TMB-H (Tumor Mutational Burden-High) or MSI-H (Microsatellite Instability-High)
  if var_type == "MSI" or var_type == "TMB":
    mgene_id = get_metaid(gene, "gene")
    mvar_ids = pymysql_cursor('SELECT ID FROM Variant WHERE Name LIKE "%s%s" AND AssemblyID IN (%s);' % (var, '%', assembly_range))
  
  # Already handled at the database construction stage: Expression, Overexpression, Underexpression, No Expression
  elif var_type == "RNA":
    mgene_id = get_metaid(gene, "gene")
    mgid_range = ",".join([str(gid) for gid in mgene_id])
    mvar_ids = pymysql_cursor('SELECT ID FROM Variant WHERE Name = "%s" AND GeneID = "%s" AND AssemblyID IN (%s);' % (var, mgid_range, assembly_range))
  
  # Copy Number Variation
  elif var_type == 'CNV':
    mgene_id = get_metaid(gene, "gene")
    mgid_range = ",".join([str(gid) for gid in mgene_id])
    if var == 'gain':
      var = "Amplification"
    elif var == 'loss':
      var = "Loss"
    else:
      var = "Non Amplification"
    mvar_ids = pymysql_cursor('SELECT ID FROM Variant WHERE Name = "%s" AND GeneID IN ("%s") AND AssemblyID IN (%s);' % (var, mgid_range, assembly_range))
    if var == "Non Amplification":
      var = "Neutral"

  # Gene Fusions
  elif var_type == "Fusion":
    mvar_ids = []
    # Gene 1
    mg_id1 = get_metaid(gene[0], "gene")
    ids1 = pymysql_cursor('SELECT ID FROM Variant WHERE (Name = "%s" OR Name = "%s") AND GeneID = "%s" AND AssemblyID IN (%s);' % (var[0], var[1], ",".join([str(gid) for gid in mg_id1]), assembly_range))
    if ids1:
      if type(ids1) == list:
        mvar_ids.extend(ids1)
      else:
        mvar_ids.append(ids1)
    # Gene 2
    mg_id2 = get_metaid(gene[1], "gene")
    ids2 = pymysql_cursor('SELECT ID FROM Variant WHERE (Name = "%s" OR Name = "%s") AND GeneID = "%s" AND AssemblyID IN (%s);' % (var[0], var[2], ",".join([str(gid) for gid in mg_id2]), assembly_range))
    if ids2:
      if type(ids2) == list:
        mvar_ids.extend(ids2)
      else:
        mvar_ids.append(ids2)
  
  # Collect the result
  if mvar_ids is None:
    mvar_ids = []
  else:
    if type(mvar_ids) != list:
      mvar_ids = [mvar_ids]
  # Generate var_range for searching Therapy
  var_range = ",".join([str(var_id) for var_id in mvar_ids])
  if var_range != "":
    therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name LIKE "%{0}%") AND ID IN (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN (SELECT ID FROM Therapy WHERE VariantID IN ({1})));'.format('Somatic', var_range))
    therapy_detail = cursor.fetchall()
    # Change the format
    out = [list(item) for item in therapy_detail]
    column_names = ['ID', 'SourceID', 'EntryID', 'Support', 'Alteration', 'TumorType', 'Therapy', 'Allele', 'Response', 'URL', 'OriginID', 'RawLevelID', 'NormLevelID', 'NormTherapyID']
    out_df = pd.DataFrame(out, columns=column_names)
    # Add Gene and Alteration
    if var_type == "Fusion":
      Gene = "-".join(gene)
      Alteration = "Fusion"
    else:
      Gene = gene
      Alteration = var
    out_df.insert(0, 'Gene', Gene)
    out_df.insert(1, 'Variant', Alteration)
    out_df.insert(2, 'Source', var_type)
    # Detected variants
    detected_var = (var_range, Gene, Alteration, var_type)
  else:
    out_df = pd.DataFrame()
    detected_var = ()
  
  return(out_df, detected_var)


def search_set(detected_var_df, assembly):
  assembly_ids = pymysql_cursor('SELECT ID FROM AssemblyDic WHERE Name IN ("%s", "na");' % assembly)
  assembly_range = ",".join([str(assembly_id) for assembly_id in assembly_ids])
  # Process the detected variants one by one
  varset_ids = []
  detected_var = detected_var_df.VariantID.astype('int').drop_duplicates().to_list()
  for index, row in detected_var_df.iterrows():
    var = row.VariantID
    # Get the related SetIDs
    set_ids = pymysql_cursor('SELECT DISTINCT SetID FROM SetHasSubset WHERE SubsetID IN (SELECT SubsetID FROM SubsetHasElement WHERE MetaID = "%s" AND Class = "Variant");' % var)
    if type(set_ids) != list:
      set_ids = [set_ids]
    # Determing wether this record is included.
    for SetID in set_ids: # set_ids = [None]时也可以查
      flag = 0
      # Dont need to get the relations of between subsets, because they are `all` or `combination`
      # Get the relation between the elements of subset
      subset_ids = pymysql_cursor('SELECT DISTINCT SubsetID FROM SetHasSubset WHERE SetID = "%s";' % SetID)
      ### One Set maps to multiple Subsets
      if type(subset_ids) is list:
        # 逐个Subset判断是否符合条件
        for subset_id in subset_ids:
          all_metaids = pymysql_cursor('SELECT DISTINCT ID FROM Variant WHERE ID IN (SELECT MetaID FROM SubsetHasElement WHERE SubsetID = "%s" AND Class = "Variant") AND AssemblyID IN (%s);' % (subset_id, assembly_range))
          if type(all_metaids) != list:
            all_metaids = [all_metaids]
          Relation = pymysql_cursor('SELECT DISTINCT Relation FROM SubsetHasElement WHERE SubsetID = "%s" AND Class = "Variant";' % subset_id)
          if type(Relation) == list:
            Relation = Relation[0]
          if Relation == "one or more":
            if list(set(all_metaids) & set(detected_var)):
              flag = 1
            else:
              flag = 0 # 因为是依次处理subset，所以可能出现这个subset可以而另一个不可以的情况，所以要让flag再等于0
          elif Relation == "not":
            if list(set(all_metaids) & set(detected_var)) == []:
              flag = 1
            else:
              flag = 0
          else:
            if (set(all_metaids) < set(detected_var)):
              flag = 1
            else:
              flag = 0
      ### One Set maps to only one Subset and this Subset maps to multiple Elements
      else:
        Relation = pymysql_cursor('SELECT DISTINCT Relation FROM SubsetHasElement WHERE SubsetID = "%s" AND Class = "Variant";' % subset_ids)
        if type(Relation) == list:
          Relation = Relation[0]
        if Relation == "one or more":
          flag = 1
        elif Relation == "not":
          flag = 0
        else: # contains all and null value
          all_metaids = pymysql_cursor('SELECT DISTINCT ID FROM Variant WHERE ID IN (SELECT MetaID FROM SubsetHasElement WHERE SubsetID = "%s" AND Class = "Variant") AND AssemblyID IN (%s);' % (subset_ids, assembly_range))
          if type(all_metaids) != list:
            all_metaids = [all_metaids]
          # Check if the variants contain all the required
          if (set(all_metaids) < set(detected_var)):
            flag = 1
      # Final judging
      if flag == 1:
        varset_ids.append([SetID, row.Gene, row.Alteration, row.Source, row.POID])
  
  ## Generate var_range for searching Therapy
  column_names = ['ID', 'SourceID', 'EntryID', 'Support', 'Alteration', 'TumorType', 'Therapy', 'Allele', 'Response', 'URL', 'OriginID', 'RawLevelID', 'NormLevelID', 'NormTherapyID']
  set_therapy = pd.DataFrame()
  for item in varset_ids:
    if item[3] == 'Germline':
      # Just one record is germline alteration in PreMedKB-POI v1.0
      therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name LIKE "%{0}%") AND ID IN (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN (SELECT ID FROM Therapy WHERE VariantSetID IN ({1})));'.format('Germline', item[0]))
    else:
      therapy_detail = cursor.execute('SELECT * FROM TherapyDetail WHERE SourceID != 4 AND OriginID IN (SELECT ID FROM OriginDic WHERE Name LIKE "%{0}%") AND ID IN (SELECT DetailID FROM TherapyHasDetail WHERE TherapyID IN (SELECT ID FROM Therapy WHERE VariantSetID IN ({1})));'.format('Somatic', item[0]))
    therapy_detail = cursor.fetchall()
    # Change the format
    out = [list(item) for item in therapy_detail]
    out_df = pd.DataFrame(out, columns=column_names)
    # Add Gene and Alteration
    out_df.insert(0, 'Gene', item[1])
    out_df.insert(1, 'Variant', item[2])
    out_df.insert(2, 'Source', item[3])
    out_df.insert(0, 'POID', item[4])
    set_therapy = pd.concat([set_therapy, out_df], axis=0)
  
  # Output
  return(set_therapy)




# def read_annova_txt(input_vcf):
#   # Read annotated txt and split by ALT & AF
#   vcf_reader = pd.read_csv(input_vcf, sep="\t")
#   vcf_single = vcf_reader[vcf_reader.ALT.str.contains(',') == False]
#   vcf_multi = vcf_reader[vcf_reader.ALT.str.contains(',')]
#   multi_split = []
#   for index, row in vcf_multi.iterrows():
#     alts = row.ALT.split(',')
#     if 'AF' in vcf_multi.columns.to_list():
#       afs = row.AF.split(',')
#     for i in range(0, len(alts)):
#       new_row = row
#       new_row.ALT = alts[i]
#       if 'AF' in vcf_multi.columns.to_list():
#         new_row.AF = afs[i]
#       multi_split.append(new_row.tolist())
# 
#   multi_split_df = pd.DataFrame(multi_split, columns=vcf_reader.columns)
#   merged_vcf = pd.concat([vcf_single, multi_split_df], axis=0)
#   return(merged_vcf)