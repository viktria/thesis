#!/usr/bin/python3

import pandas as pd
import numpy as np

data = []
with open("ipr_with_2_domains_filter.dat", "r") as file:
    for line in file:
        line = line.strip().split("\t")
        data.append((line[0], line[1], line[2], line[3], int(line[4]), int(line[5])))

df = pd.DataFrame(data, columns=['UniProt ID', 'Domain_ID', 'Name', 'SSF', 'start', 'stop'])
df['Distance'] = df.stop - df.start
df["Distance"] = pd.to_numeric(df["Distance"], errors="coerce")
grouped_df = df.groupby("Domain_ID")['Distance']

# Define custom aggregation functions
agg_funcs = {
    'Mean_Distance': lambda x: np.nanmean(x),
    'Std_Distance': lambda x: np.nanstd(x),
    'Count': lambda x: (~np.isnan(x)).sum(),
    'Std/Mean_Distance': lambda x: np.nanstd(x) / np.nanmean(x) if np.nanmean(x) != 0 else np.nan
}

# Calculate the mean, standard deviation, and count for each group, excluding NaN values
grouped_stats = grouped_df.agg(**agg_funcs)

# Rename the columns for clarity
grouped_stats.columns = ["Mean_Distance", "Std_Distance", "Count", "Std/Mean_Distance"]

# Reset index to make the group keys accessible as columns
grouped_stats.reset_index(inplace=True)

#grouped_stats.to_csv("domain_stats_pdb.tsv", sep="\t", index=False)
grouped_stats = grouped_stats[grouped_stats["Std/Mean_Distance"]<= 0.2]
df = df.drop(columns=['Distance'])
m = df.Domain_ID.isin(grouped_stats.Domain_ID)
merged_df = df[m]
merged_df.to_csv('ipr_with_2_domains_filter_02.dat', sep='\t', index=False, header=False)
    
        