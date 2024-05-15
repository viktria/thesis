#!/usr/bin/python3

import pandas as pd
import numpy as np

# Read the first dataset
df1 = pd.read_csv("best_pdb_full.txt", sep=' ', header=None, index_col=False,
                        names=['SP_PRIMARY','PDB','CHAIN','Domain_1','Start_1','Stop_1','Domain_2','Start_2','Stop_2'])

# Read the second dataset
df2 = pd.read_csv("distance_pdb.tsv", sep="\t")

# Merge the datasets based on common columns
pdb_df = pd.merge(df1, df2, left_on=["PDB", "SP_PRIMARY",'CHAIN', "Domain_1", "Domain_2","Start_1"], right_on=["PDB ID", "UniProt ID",'CHAIN', "Domain_1", "Domain_2",'Merge_ID'], how="inner")

# Drop unnecessary columns
pdb_df.drop(columns=["PDB ID", "UniProt ID"], inplace=True)
pdb_df = pdb_df.loc[:, ~pdb_df.columns.duplicated()]
pdb_df["Distance"] = pd.to_numeric(pdb_df["Distance"], errors="coerce")

# Group the data by Domain_1 and Domain_2
grouped_df = pdb_df.groupby(["Domain_1", "Domain_2"])

# Calculate the mean and standard deviation for each group, excluding NaN values
grouped_stats = grouped_df.agg({"Distance": [lambda x: np.nanmean(x), lambda x: np.nanstd(x), lambda x: (~np.isnan(x)).sum()]})

# Rename the columns for clarity
grouped_stats.columns = ["Mean_Distance", "Std_Distance", "Count"]

# Reset index to make the group keys accessible as columns
grouped_stats.reset_index(inplace=True)

# Write the grouped statistics to a new TSV file
grouped_stats.to_csv("grouped_stats_pdb.tsv", sep="\t", index=False)

