import pandas as pd
from string import ascii_lowercase
import xlrd
import sys

#--------
# Creates nested dictionary of drug names (key) to associated phecodes (values) in MEDI
# outer keys = first letter of drug name, inner keys = drug name, values = list of phecodes
# PARAMETERS:
# medi = dataframe object of the medi spreadsheet
def create_dic(medi):
	dic = {k: g["PheCode"].tolist() for k,g in medi.groupby("DRUG_DESC")}

	nested_dic = {}
	for c in ascii_lowercase:
		nested_dic[c] = {}
	for key in dic.keys():
		for c in ascii_lowercase:
			if key[0] == c:
				nested_dic[c][key] = dic[key]

	return nested_dic

#--------
# Iterates through each drug-phecode pair (that was analyzed for association in GWAS) and checks
# if the pair is in MEDI (set IN_MEDI column value to 1) or not (set to 0)
# Imports all pairs and their documentation in MEDI to .csv file
# PARAMETERS: 
# drug_GWASphenotype_pairs = contains p values for each drug-GWAS phenotype pair
# phecode_GWASphenotype_pairs = contains mapping of phecode to GWAS phenotype
# drug_clinicalphenotype_pairs = contains drug-clinical phenotype pairs 
def finding_matches(drug_GWASphenotype_pairs, phecode_GWASphenotype_pairs, 
					drug_clinicalphenotype_pairs):
	pairs = pd.read_table(drug_GWAS_pairs, header=0, sep=',', 
		usecols=['phenotype', 'drug', 'p'])
	pairs['phenotype'] = pairs['phenotype'].str.lower()
	
	# database mapping gwas phenotype ('disease' column) to phecode
	gwas_phecode = pd.read_excel(phecode_GWASphenotype_pairs, usecols=['disease', 'jd_string', 'jd_code']).drop_duplicates().dropna(how='any', axis=0)
	gwas_phecode['disease'] = gwas_phecode['disease'].str.lower()
	gwas_phecode.columns = ['phenotype', 'jd_string', 'jd_code']
	#gwas_phewas = pd.read_csv('gene_phenotype.csv',		# merged gwas phenotype with phewas phenotype (phecode) by SNP
		#usecols=['phenotype', 'jd_code']).dropna(how='any', axis=0) 	# not used since deemed less accurate

	# MEDI database including associated phecodes for each drug
	medi_phecode = pd.read_csv(drug_clinicalphenotype_pairs, 
		usecols=['DRUG_DESC', 'PheCode']).drop_duplicates().dropna(how='any', axis=0)
	medi_phecode['DRUG_DESC'] = medi_phecode['DRUG_DESC'].str.lower()

	# merge drug-phenotype pairs with phecodes corresponding by phenotype
	# also added in column 'IN_MEDI' to record whether pair is in medi_phecode or not
	pairs_phewas = pd.merge(pairs, gwas_phecode, how='left', on='phenotype').drop_duplicates()
	pairs_phewas['IN_MEDI'] = 0

	# calculates number of unique drugs in MEDI
	all_drugs_medi = medi_phecode['DRUG_DESC'].drop_duplicates()

	# calculates number of unique drugs in all-drugs-gwas
	all_drugs_pairs = pairs['drug'].drop_duplicates()

	# reduces MEDI database to only the drugs that were examined by GWAS
	medi_phecode_reduced = medi_phecode[medi_phecode['DRUG_DESC'].isin(all_drugs_pairs)].sort_values('DRUG_DESC')
	# pairs_phewas_reduced = pairs_phewas[pairs_phewas['drug'].isin(all_drugs_medi)].sort_values('drug').dropna(how='any', axis=0)

	# converts medi-phecode pairs into a dictionary
	medi_dic = create_dic(medi_phecode_reduced)

	# checks whether GWAS pairing is in MEDI
	for index, row in pairs_phewas.iterrows():
		for key in medi_dic.keys():
			if row['drug'][0] == key:
				for key2 in medi_dic[key].keys():
					if row['drug'] == key2:
						if row['jd_code'] in medi_dic[key][key2]:
							pairs_phewas.set_value(index, 'IN_MEDI', 1)

	print(pairs_phewas)
	pairs_phewas.to_csv('in_medi.csv', index=False)
