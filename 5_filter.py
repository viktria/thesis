#!/usr/bin/python3
import pandas as pd

# Read the data into a DataFrame
df = pd.read_csv('tmp', sep='\t', header=None, names=['UniProtID', 'domain1'])
df1 = pd.read_csv('ipr_start_stop.txt', sep=',', names=['UniProt_ID', 'domain_1', 'domain_2', 'start_1', 'stop_1', 'start_2', 'stop_2'])

# Merge the two dataframes based on the UniProt ID
merged_df = pd.merge(df, df1, left_on=['UniProtID','domain1'], right_on=['UniProt_ID','domain_1'])

# Reorder the columns
merged_df = merged_df[['UniProtID', 'domain_1','start_1', 'stop_1', 'domain_2','start_2', 'stop_2']]

# Save the mixed dataframe
merged_df.to_csv('ipr_start_stop_filtered.txt', sep=',', index=False,header=None)






