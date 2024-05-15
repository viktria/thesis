#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from sklearn.linear_model import LinearRegression
import math

# Histogram of Standard Deviations
grouped_stats_pdb_1 = pd.read_csv("grouped_stats_pdb.tsv", sep="\t")
grouped_stats_pdb = grouped_stats_pdb_1[grouped_stats_pdb_1['Mean_Distance'] != 0]
grouped_stats_af_1 = pd.read_csv("grouped_stats_af.tsv", sep="\t")
grouped_stats_af = grouped_stats_af_1[grouped_stats_af_1['Mean_Distance'] != 0]
plt.figure(figsize=(8, 6))
plt.hist(grouped_stats_af["Std_Distance"], bins=4, color='indianred', alpha=0.7, label='AlphaFold')
plt.hist(grouped_stats_pdb["Std_Distance"], bins=10, color='cornflowerblue', alpha=0.6, label='PDB')
plt.xlabel('Standard Deviation')
plt.ylabel('Frequency')
plt.title('Histogram of Standard Deviations')
plt.legend()
plt.grid(True)
plt.savefig('/bigdisk/users/kelvi/results/combined_std_histogram.png', dpi=300)
plt.close()

# Standard Deviation of Distances
df1 = grouped_stats_pdb[grouped_stats_pdb['Mean_Distance'] != 0]
df2 = grouped_stats_af
df_distances = pd.read_csv('best_pdb_full.txt', sep=' ', header=None, index_col=False,
                        names=['SP_PRIMARY','PDB','CHAIN','Domain_1','Start_1','Stop_1','Domain_2','Start_2','Stop_2'])
# Merge the distance DataFrame with the first two based on common identifier
merged_df = pd.merge(df_distances, df1, on=['Domain_1', 'Domain_2'], how='inner')
merged_df = pd.merge(merged_df, df2, on=['Domain_1', 'Domain_2'], how='inner')
# Calculate distances between domains
merged_df['Distance'] = merged_df['Start_2'] - merged_df['Stop_1']
filtered_df = merged_df[(merged_df['Std_Distance_x'] > 4) & (merged_df['Std_Distance_y'] < 2.5)]
# Save filtered data to a file
filtered_df.to_csv('/bigdisk/users/kelvi/results/b_pdb_s_af_distances.tsv', sep='\t', index=False)
# Extract relevant columns containing standard distances
std_distances_x = merged_df['Std_Distance_x']
std_distances_y = merged_df['Std_Distance_y']
std_distances_x = std_distances_x.to_numpy()
std_distances_y = std_distances_y.to_numpy()
# Reshape std_distances_x to a 2D array
std_distances_x = std_distances_x.reshape(-1, 1)
# Fit a line using Linear Regression from Scikit-learn
model = LinearRegression()
model.fit(std_distances_x, std_distances_y)
# Plot the scatterplot and the fitted line
plt.figure(figsize=(8, 6))
plt.scatter(std_distances_x, std_distances_y, c=abs(merged_df['Distance']), cmap='coolwarm', alpha=0.5, norm=LogNorm())
plt.colorbar(label='Distance')
plt.plot([0, std_distances_x.max()], [0, std_distances_x.max() * model.coef_[0]], color='indianred', label=f'Slope: {model.coef_[0]:.2f}')
plt.axis([0, 13, 0, 13])
plt.title('Standard Deviation of Distances')
plt.xlabel('PDB')
plt.ylabel('AlphaFold')
plt.grid(True)
plt.legend()
plt.savefig('/bigdisk/users/kelvi/results/std_compare.png', dpi=300)
plt.close()

# Distance Comparison
df3_a = pd.read_csv("distance_pdb.tsv", sep="\t")
df3 = df3_a[df3_a['Distance'] != 'Not Found']
df4 = pd.read_csv("distance_af.tsv", sep="\t")
m_df = pd.merge(df3, df4, on=["UniProt ID", "PDB ID","Merge_ID"], suffixes=('_x', '_y'))
domain_counts = {}
for _, row in m_df.iterrows():
    key = (row['Domain_1_x'], row['Domain_2_x'])
    domain_counts[key] = domain_counts.get(key, 0) + 1
color_dict = {}
color_index = 0
for key, count in domain_counts.items():
    if count == 1:  # Domain pair appears only once, assign special color
        color_dict[key] = -1  # Assign a special value for special color
    else:
        color_dict[key] = color_index
        color_index += 1
# Create a new column in the DataFrame to store colors
m_df['Color'] = [color_dict[(row['Domain_1_x'], row['Domain_2_x'])] for _, row in m_df.iterrows()]
m_df.sort_values(by='Distance_x', inplace=True)
# Plotting
plt.figure(figsize=(8, 6))
for color, group in m_df.groupby('Color'):
    if color == -1:  # Gray color for domain pairs appearing only once
        plt.scatter(group['Distance_x'].astype(float), group['Distance_y'].astype(float), label='Single Occurrence', color='red',marker='x', alpha=0.7)
    else:
        plt.scatter(group['Distance_x'].astype(float), group['Distance_y'].astype(float), cmap='plasma', alpha=0.7)
plt.xlabel("Distances from the PDB calculations")
plt.ylabel("Distances from the AlphaFold calculations")
plt.title("Comparison of Distances")
plt.legend()
plt.grid(True)
plt.tick_params(axis='x', rotation=90)
plt.tight_layout()
plt.savefig('/bigdisk/users/kelvi/results/distance_compare.png', dpi=300)
plt.close()