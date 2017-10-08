import pandas as pd
from string import ascii_lowercase
import xlrd
import numpy as np
import json
from scipy import stats
from math import isinf
import sys

#--------
# Returns list of drugs studied in gwas/inrich dataset
# PARAMETERS:
# gwas_dataset = gwas/inrich dataset containing all the drug-phenotype pairs
def get_drug_list(gwas_dataset):
	data = pd.read_table(gwas_dataset, header=0, sep='	')
	drugs = data['drug'].drop_duplicates().sort_values().tolist()

	return drugs

#--------
# Creates and returns dictionary with drug name (key) to hash value (value)
# Creates a 'value_to_name.json' file containing a dictionary with hash value (key) to drug name (value)
# PARAMETERS:
# gwas_dataset = gwas/inrich dataset containing all the drug-phenotype pairs
def create_match(gwas_dataset):
	drugs = get_drug_list(gwas_dataset)
	name_to_value = {}
	value_to_name = {}
	for i in range(987):
		name_to_value[drugs[i]] = i
		value_to_name[i] = drug[i]

	with open('value_to_name.json', 'w') as f:
		json.dump(value_to_name, f)

	return name_to_value

#--------
# Replaces drug names with hash values in the input dataframe, returns resulting dataframe
# PARAMETERS:
# gwas_dataset = gwas/inrich dataset containing all the drug-phenotype pairs
# medi_df = medi dataset in dataframe form
def replace_w_hash(gwas_dataset, medi_df):
	#df = pd.read_table('all-drugs-gwas.results.txt', header=0, sep='	')
	match = create_match(gwas_dataset)
	result = medi_df.replace({'drug': match})
	
	return result

#--------
# Finds significant pairs that correspond with an input p-value, returns a dataframe with only these pairs
# PARAMETERS:
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
# p = p-value cutoff
def significant_pairs(p, in_medi_dataset):
	df = pd.read_csv(in_medi_dataset).sort_values(['phenotype', 'drug'])

	# find # of significant pairs (763 for 0.05), then find how many of those are in MEDI (191 for 0.05)
	sig = df[df['p'] <= p].filter(['phenotype', 'jd_code', 'drug', 'p', 'pair_in_medi'])
	# sig_medi = sig[sig['pair_in_medi'] >= 1]

	return sig

#--------
# Creates and returns dictionary of drug hash value (key) to associated phecodes in MEDI
# PARAMETERS:
# drug_phenotype = contains drugs (in form of name) and their associated phecodes
# gwas_dataset = gwas/inrich dataset containing all the drug-phenotype pairs
def create_dict(drug_phenotype, gwas_dataset):
	medi = pd.read_csv(drug_phenotype)
	drugs = get_drug_list(gwas_dataset)

	medi['DRUG_DESC'] = medi['DRUG_DESC'].str.lower()

	# drops all rows where the drug is not in the list of drugs studied in gwas/inrich dataset
	medi_reduced = medi[medi['DRUG_DESC'].isin(drugs)].filter(['DRUG_DESC', 'PheCode']).dropna(how='any', axis=0)
	medi_reduced.columns = ['drug', 'jd_code']
	medi_reduced['jd_code'] = medi_reduced['jd_code'].astype(str)
	
	hash_medi = replace_w_hash(medi_reduced)
	dic = {k: g["jd_code"].tolist() for k,g in hash_medi.groupby("drug")}

	return dic

#--------
# Randomize drug column and finds how many of the permuted pairs are in MEDI
# PARAMETERS: 
# p = p-value cutoff
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
# drug_phenotype = contains drugs (in form of name) and their associated phecodes
# gwas_dataset = gwas/inrich dataset containing all the drug-phenotype pairs
def randomize(p, in_medi_dataset, gwas_dataset, drug_phenotype):
	# gets significant drug-phenotype pairs
	sig = significant_pairs(p, in_medi_dataset)
	sig = replace_w_hash(gwas_dataset, sig)
	
	# gets MEDI dictionary
	my_dict = create_dict(drug_phenotype, gwas_dataset)

	# create dataframe where drugs are randomly permuted
	random = sig.filter(['jd_code', 'drug', 'p'], axis=1)
	random['IN_MEDI'] = 0
	random['drug'] = np.random.permutation(random['drug'])

	# takes out pairs where phenotype has no jd_codes
	nan = random[random['jd_code'].isnull()]
	random = random.dropna(how='any', axis=0)

	# iterates through rest of rows to find whether pair is in MEDI
	for index, row in random.iterrows():
		for key in my_dict.keys():
			if row['drug'] == key:
				for i in range(len(row['jd_code'])):
					if row['jd_code'][i] in my_dict[key]:
						random.set_value(index, 'IN_MEDI', 1)
						break

	result = random.append(nan).sort_index()
	medi = result[result['IN_MEDI'] == 1]
	num_pairs = medi.shape[0]

	return num_pairs

#--------
# Runs randomized permutation ('randomize' function) n number of times, puts output in a list
# PARAMETERS:
# n = number of times to run permutation
# p = p-value cutoff
# in_medi_dataset = contains information about documentation of drug-GWAS phenotype pairs in MEDI
# drug_phenotype = contains drugs (in form of name) and their associated phecodes
# gwas_dataset = gwas/inrich dataset containing all the drug-phenotype pairs
def analysis(n, p, in_medi_dataset, gwas_dataset, drug_phenotype):
	data = []
	for i in range(n):
		x = randomize()
		data.append(x)

	with open('observed.json', 'w') as f:
		json.dump(data, f)




	




