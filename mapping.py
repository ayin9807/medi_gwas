import pandas as pd
from pandas import DataFrame
import sys

#--------
# Imports medi dataset with icd9 and rxcui descriptions to .csv file
# PARAMETERS:
# medi = medi spreadsheet
# icd9_desc = contains icd9 codes and their descriptions
# rxcui_desc = contains rxcui codes and their descriptions
def add_info_to_medi(medi, icd9_desc, rxcui_desc):
	# adding in icd9 descriptions
	df_icd9_desc = pd.read_table(icd9_desc, sep='	', header=None, usecols=[0, 1])
	df_icd9_desc.columns = ['ICD9', 'ICD9_DESC']

	# adding in rxcui descriptions into the medi spreadsheet
	df_rxcui_desc = pd.read_csv(rxcui_desc, encoding='latin-1').drop_duplicates().groupby('RXCUI_IN')['STR'].apply('; '.join)
	rxcui_desc = pd.DataFrame({'RXCUI_IN': df_rxcui_desc.index, 'STR': df_rxcui_desc.values})

	df_medi = pd.read_csv(medi)

	df_medi_desc = pd.merge(df_medi, rxcui_desc, how='left', on='RXCUI_IN')
	df_rxcui_icd9 = pd.merge(df_medi_desc, df_icd9_desc, how='left', on='ICD9')
	df_rxcui_icd9 = df_rxcui_icd9[['RXCUI_IN', 'STR', 'DRUG_DESC', 'ICD9', 'ICD9_DESC', 'INDICATION_DESCRIPTION', 'MENTIONEDBYRESOURCES', 
			'HIGHPRECISIONSUBSET', 'POSSIBLE_LABEL_USE']]

	df_rxcui_icd9.to_csv('medi_with_icd9_rxcui.csv', index=False)

#--------
# Imports medi_rxcui_icd9 dataset with icd9-phecode mappings to .csv file
# Maps drug (rxcui codes) with clinical phenotype (phecode) through icd9 codes
# PARAMETERS:
# medi_rxcui_icd9 = medi spreadsheet (created from add_info_to_medi function above) with rxcui + icd9 descriptions
# phecode_icd9_mapping = maps phecodes to icd9 codes
def drug_phenotype(phecode_icd9_mapping, medi_rxcui_icd9):
	df_rxcui_icd9 = pd.read_csv(medi_rxcui_icd9)
	df_phecode_icd9 = pd.read_csv(phecode_icd9_mapping, usecols=['ICD9', 'PheCode'])

	result = pd.merge(df_rxcui_icd9, df_phecode_icd9, how='left', on='ICD9').drop_duplicates().sort_values('RXCUI_IN')

	result.to_csv('drug_phenotype.csv', index=False)
	#print (result)

#--------
# Imports medi_rxcui_icd9 dataset with drug-targeted gene mappings to .csv file
# Maps drugs (rxcui codes) with corresponding targeted genes (HuGOIDs) through unii codes and DrugBank drug IDs
# PARAMETERS:
# unii_rxcui = contains mapping of unii codes to rxcui codes
# unii_drug = contains mapping of unii codes to HuGOIDs (DrugBank), needs to be .txt file
# medi_rxcui_icd9 = medi spreadsheet (created from add_info_to_medi function above) with rxcui + icd9 descriptions
# drug_gene = for each gene, contains list of drugs that target said gene
def drug_gene(unii_rxcui, unii_drug, drug_gene, medi_rxcui_icd9):
	df_unii_rxcui = pd.read_csv(unii_rxcui)
	df_unii_drug = pd.read_table(unii_drug, header=0, sep=':', usecols=['unii', 'drug_id'])
	df_rxcui_icd9 = pd.read_csv(medi_rxcui_icd9)
	
	# drugbank id and rxcui mapping
	data1 = pd.merge(df_unii_drug, df_unii_rxcui, how='left', 
		on='unii').drop('unii', axis=1).drop_duplicates()

	# splits drugs for each gene in individual cell
	data2 = pd.read_csv(drug_gene, usecols=['Drug IDs', 'Gene Name'])
	df_drugbank_gene = DataFrame(data2['Drug IDs'].str.split('; ').tolist(), 
		index=data2['Gene Name']).stack().reset_index()[[0, 'Gene Name']] # var1 variable is currently labeled 0
	df_drugbank_gene.columns = ['drug_id', 'Gene Name']
	df_drugbank_gene = df_drugbank_gene.dropna(how='any', axis=0)

	# for each drug combines all targeted genes into one cell
	data3 = df_drugbank_gene.drop_duplicates().groupby('drug_id')['Gene Name'].apply('; '.join)
	data4 = pd.DataFrame({'drug_id': data3.index, 'Gene Name': data3.values})

	drug_rxcui = pd.merge(data1, data4, how='left', on='drug_id').drop_duplicates()

	result = pd.merge(df_rxcui_icd9, drug_rxcui, how='left', on='RXCUI_IN')
	result.to_csv('drug_gene.csv', index=False)
	#print (result)

#--------
# Imports dataset with mappings between gwas phenotype and clinical phenotype through snp
# Merges gwas phenotype with phewas phenotype (phecode) by SNP
''' ** was not used when creating all-drugs-gwas-data.csv, mapping by SNP was deemed to be not
	as accurate as using gwas_catalog-phewas.csv, which already had gwas phenotype (disease column)
	mapped to phecode '''
# PARAMETERS:
# gwas = contains gwas/inrich associations between drug-phenotype pairs
# phewas = contains mapping between phecodes and gwas phenotypes (by name)
def gene_phenotype(gwas, phewas):
	df_gwas = pd.read_table(gwas, header=0, sep=' ')
	df_phewas = pd.read_csv(phewas, usecols=['snp', 'jd_code', 'jd_string'])

	result = pd.merge(df_gwas, df_phewas, how='left', on='snp').drop_duplicates()
	result.to_csv('gene_phenotype.csv', index=False)
	#print (result)



