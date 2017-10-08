The following scripts were used to find potential new disease targets for existing medications using GWAS/INRICH and MEDI data. Updated 10/07/2017.

- mapping.py was used to create the following three datasets:
1) drug_phenotype.csv - mapped RxCUI codes in MEDI to PheCodes by ICD9
2) drug_gene.csv - mapped RxCUI codes to targeted genes (HuGOIDs) through UNII codes and DrugBank IDs
3) gene_phenotype.csv - mapped GWAS phenotype to PheCode by SNP

- accessing_pair.py iterated through all the drug-PheCode pair analyzed by GWAS and checked if each one was documented in MEDI

- condense.py 
1) condensed drug-PheCode pairs back down to drug-phenotype pairs
2) added number of PheCodes for each phenotype
3) added whether any of the drugs studied were in MEDI

- add_info.py added jd_string and num_overlap from original GWAS dataset 

- random_permutation.py took all significant pairs by a certain threshold (generally 0.05) and randomly permuted drug and checked if newly assigned pairs were in MEDI
