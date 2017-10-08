import pandas as pd
import xlrd
import ast
import sys

#--------
# Adds jd_string to entire in_medi dataset and imports resulting dataframe to .csv file
# PARAMETERS:
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
# gwas_phewas_dataset = contains mapping of GWAS phenotype, PheCode, and PheWAS phenotype
def add_jdstring(in_medi_dataset, gwas_phewas_dataset):
	df = pd.read_csv(in_medi_dataset)
	jd_string = pd.read_excel(gwas_phewas_dataset, 
		usecols=['jd_string', 'jd_code']).dropna(subset=['jd_code'], axis=0).drop_duplicates()
	jd_string['jd_code'] = jd_string['jd_code'].astype(str)

	# separates rows with phecode and rows without
	no_jd = df[df['jd_code'].isnull()]
	no_jd['jd_string'] = ""
	with_jd = df.dropna(subset=['jd_code'], axis=0)
	with_jd['jd_string'] = ""

	unique_phenotypes = with_jd.filter(['phenotype', 'jd_code', 'jd_string']).drop_duplicates()

	# iterates through jd_codes for each phenotype and adds corresponding jd_string
	for index1, row1 in unique_phenotypes.iterrows():
		jd_codes = ast.literal_eval(row1['jd_code'])
		temp = []
		for value in jd_codes:
			for index2, row2 in jd_string.iterrows():
				if value == row2['jd_code']:
					temp.append(row2['jd_string'])
					break

		phenotypes = ', '.join(temp)
		unique_phenotypes.set_value(index1, 'jd_string', phenotypes)

	unique_phenotypes = unique_phenotypes.drop('jd_code', axis=1)
	data = pd.merge(df, unique_phenotypes, how='left', on='phenotype')

	data.to_csv('all-drugs-gwas-data.csv', index=False)

#--------
# Adds num_overlap from GWAS findings for each drug-phenotype pair to database
# Imports resulting dataframe to .csv file
# PARAMETERS:
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
# drug_GWASphenotype_pairs = contains p values for each drug-GWAS phenotype pair
def add_overlapping(in_medi_dataset, drug_GWASphenotype_pairs):
	df = pd.read_csv(in_medi_dataset).sort_values(['phenotype', 'drug'])
	df2 = pd.read_table(drug_GWASphenotype_pairs, header=0, sep='	').sort_values(['phenotype', 'drug'])
	df['num_overlap'] = df2['num_overlap']

	df.to_csv('all-drugs-gwas-data.csv', index=False)
	