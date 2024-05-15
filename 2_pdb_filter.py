#!/usr/bin/python3

import pandas as pd

# Read domain data from txt file
domain_df = pd.read_csv('ipr_start_stop.txt', sep=',', header=None, 
                        names=['UniProt ID', 'Domain 1', 'Domain 2', 'Start1', 'Stop1', 'Start2', 'Stop2'])

# Read the first TSV file
tsvfile = pd.read_csv('pdb_chain_uniprot.tsv', sep='\t', skiprows=1, usecols=['PDB', 'CHAIN', 'SP_PRIMARY'] )

# Merge the two DataFrames based on 'SP_PRIMARY'
merged_df = pd.merge(tsvfile, domain_df, left_on='SP_PRIMARY', right_on='UniProt ID')

# Filter rows based on conditions
selected = merged_df[['PDB', 'CHAIN', 'SP_PRIMARY']].drop_duplicates()
selected2 = merged_df[['PDB', 'CHAIN', 'SP_PRIMARY','Domain 1', 'Start1', 'Stop1', 'Domain 2', 'Start2', 'Stop2']].drop_duplicates()

# Write the filtered rows to a new TSV file
selected.to_csv('filtered_pdb_chain_uniprot.tsv', sep='\t', index=False)
selected2.to_csv('filtered_pdb_chain_uniprot_domain.tsv', sep='\t', index=False)
