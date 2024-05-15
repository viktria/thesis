#!/usr/bin/python3
import pandas as pd

# Read the data into a DataFrame
df = pd.read_csv('ipr_start_stop_filtered.txt', sep=',', header=None, names=['UniProtID', 'Domain_1','Start_1','Stop_1','Domain_2','Start_2','Stop_2'])
df1 = pd.read_csv('refiltered_pdb_chain_uniprot.tsv', sep='\t')
df1_selected = df1[['PDB', 'CHAIN', 'SP_PRIMARY']]

# Merge the two dataframes based on the UniProt ID
merged_df = pd.merge(df, df1_selected, left_on='UniProtID', right_on='SP_PRIMARY')

# Drop the redundant SP_PRIMARY column
merged_df.drop(columns=['UniProtID'], inplace=True)

# Reorder the columns
merged_df = merged_df[['SP_PRIMARY','PDB','CHAIN','Domain_1','Start_1','Stop_1','Domain_2','Start_2','Stop_2']]

# Save the mixed dataframe
merged_df.to_csv('/bigdisk/users/kelvi/unp_pdb_dom.txt', sep=' ', index=False)






