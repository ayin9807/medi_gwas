import pandas as pd
import xlrd
import ast
import sys

#--------
# Taking phenotype with multiple jd_codes and condensing back down to number of drug-GWAS phenotype pairs
# Sums up number of matches for each drug-phenotype pair (can be > 1 since 1 phenotype can have multiple jd_codes) 
# Imports condensed dataset to .csv file
# PARAMETERS:
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
def condensing_medi(in_medi_dataset):
	df = pd.read_csv(in_medi_dataset)
	short = df[['phenotype','drug','IN_MEDI']].groupby(['phenotype','drug']).sum()
	short.to_csv('all-drugs-gwas-medi.csv')

#--------
# Counts number of jd_codes for each phenotype, sets count in 'count' column
# Returns resulting dataframe
# PARAMETERS:
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
def count_jdcode(in_medi_dataset):
	df = pd.read_csv(in_medi_dataset)
	data1 = df.filter(['phenotype', 'jd_code'], axis=1).drop_duplicates()

	# drops all phenotypes with no jd_nodes
	no_codes = data1[data1['jd_code'].isnull()].drop_duplicates()
	no_codes['count'] = 0

	data1 = data1.dropna(how='any', axis=0)
	data1['jd_code'] = data1['jd_code'].astype(str)
	#data1.to_csv('phenotypes_nona.txt', index=False, sep='	')

	# for each phenotype, group all phecodes and put into a list	
	temp = data1.groupby('phenotype')['jd_code'].apply(list)
	jd_list = pd.DataFrame({'phenotype':temp.index, 'jd_code':temp.values})
	jd_list['count'] = 0
	
	# count number of phecodes for each phenotype
	for index, row in jd_list.iterrows():
		jd_list.set_value(index, 'count', len(row['jd_code']))

	# reads phenotypes with no jd_codes
	result = jd_list.append(no_codes).sort_values('phenotype').reset_index(drop=True)
	return result

#--------
# Creates dataframe with 2 additional columns: 'drug' = drug name, 'drug_in_medi' = 1 if drug is in MEDI
# Returns resulting dataframe
# PARAMETERS:
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
def drug_medi(in_medi_dataset):
	df = pd.read_csv(in_medi_dataset)
	drugs = df['drug'].drop_duplicates().sort_values().tolist()
	data = pd.DataFrame({'drug': drugs})
	data['drug_in_medi'] = 0

	medi = pd.read_csv('MEDI_01212013.csv').filter(['DRUG_DESC'], axis=1).drop_duplicates()
	medi['DRUG_DESC'] = medi['DRUG_DESC'].str.lower()
	medi_list = medi['DRUG_DESC'].tolist()

	data.ix[data.drug.isin(medi_list), 'drug_in_medi'] = 1
	return data



	




