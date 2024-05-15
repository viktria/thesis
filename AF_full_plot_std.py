#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, kruskal, ks_2samp

# Load data
af_df = pd.read_csv("distance_full_af.tsv", sep="\t")
af_df["Distance"] = pd.to_numeric(af_df["Distance"], errors="coerce")
# Group the data by Domain_1 and Domain_2
grouped_df = af_df.groupby(["Domain_1", "Domain_2"])
# Calculate the mean and standard deviation for each group, excluding NaN values
grouped_stats = grouped_df.agg({"Distance": [lambda x: np.nanmean(x), lambda x: np.nanstd(x), lambda x: (~np.isnan(x)).sum()]})
# Rename the columns for clarity
grouped_stats.columns = ["Mean_Distance", "Std_Distance", "Count"]
# Reset index to make the group keys accessible as columns
grouped_stats.reset_index(inplace=True)
# Write the grouped statistics to a new TSV file
grouped_stats.to_csv("grouped_stats_af_full.tsv", sep="\t", index=False)

grouped_stats_af_full = pd.read_csv("grouped_stats_af_full.tsv", sep="\t")
grouped_stats_af_full = grouped_stats_af_full[grouped_stats_af_full['Count']!=1]
# Group by Std_Distance and sum the counts
grouped_full = grouped_stats_af_full.groupby('Std_Distance')['Count'].sum().reset_index()
# Calculate total count
total_count_full = grouped_full['Count'].sum()
# Calculate relative frequencies
grouped_full['Relative_Frequency'] = grouped_full['Count'] / total_count_full
grouped_full.sort_values(by=['Relative_Frequency'])
grouped_full= grouped_full[grouped_full['Std_Distance'] <= 10]
grouped_full= grouped_full[grouped_full['Relative_Frequency'] <= 0.025]

grouped_stats_af = pd.read_csv("grouped_stats_af.tsv", sep="\t")
grouped_stats_af = grouped_stats_af[grouped_stats_af['Count']!=1]
# Group by Std_Distance and sum the counts
grouped = grouped_stats_af.groupby('Std_Distance')['Count'].sum().reset_index()
# Calculate total count
total_count = grouped['Count'].sum()
# Calculate relative frequencies
grouped['Relative_Frequency'] = grouped['Count'] / total_count
grouped.sort_values(by=['Relative_Frequency'])
grouped= grouped[grouped['Std_Distance'] <= 10]

grouped_stats_pdb = pd.read_csv("grouped_stats_pdb.tsv", sep="\t")
grouped_stats_pdb = grouped_stats_pdb[grouped_stats_pdb['Count']!=1]
# Group by Std_Distance and sum the counts
grouped_pdb = grouped_stats_pdb.groupby('Std_Distance')['Count'].sum().reset_index()
# Calculate total count
total_count_pdb = grouped_pdb['Count'].sum()
# Calculate relative frequencies
grouped_pdb['Relative_Frequency'] = grouped_pdb['Count'] / total_count_pdb
grouped_pdb.sort_values(by=['Relative_Frequency'])
grouped_pdb= grouped_pdb[grouped_pdb['Std_Distance'] <= 10]

fig, ax = plt.subplots(figsize=(8, 6))

# First violin plot
p1 = ax.violinplot(grouped_full['Std_Distance'], positions=[1], widths=0.8,
              showmeans=True, showextrema=True, showmedians=True)
for pc in p1['bodies']:
    pc.set_facecolor('cornflowerblue')
    pc.set_edgecolor('cornflowerblue')
p1['cmedians'].set_colors('cornflowerblue')
p1['cmeans'].set_colors('cornflowerblue')
p1['cmins'].set_colors('cornflowerblue')
p1['cmaxes'].set_colors('cornflowerblue')
p1['cbars'].set_colors('cornflowerblue')
# Second violin plot
p2 = ax.violinplot(grouped['Std_Distance'], positions=[2], widths=0.8,
              showmeans=True, showextrema=True, showmedians=True)
for pc in p2['bodies']:
    pc.set_facecolor('indianred')
    pc.set_edgecolor('indianred')
p2['cmedians'].set_colors('indianred')
p2['cmeans'].set_colors('indianred')
p2['cmins'].set_colors('indianred')
p2['cmaxes'].set_colors('indianred')
p2['cbars'].set_colors('indianred')
# Third violin plot
p3 = ax.violinplot(grouped_pdb['Std_Distance'], positions=[3], widths=0.8,
              showmeans=True, showextrema=True, showmedians=True)
for pc in p3['bodies']:
    pc.set_facecolor('mediumseagreen')
    pc.set_edgecolor('mediumseagreen')
p3['cmedians'].set_colors('mediumseagreen')
p3['cmeans'].set_colors('mediumseagreen')
p3['cmins'].set_colors('mediumseagreen')
p3['cmaxes'].set_colors('mediumseagreen')
p3['cbars'].set_colors('mediumseagreen')

ax.set_title('Violin plot of Standard Deviations of two AlphaFold and the PDB sets', fontsize=10)

# Set x-axis ticks and labels
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['AF full set', 'AF PDB set', 'PDB set'])

# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig('/bigdisk/users/kelvi/results/violin_plots_3.png', dpi=300)

# Show the plots
plt.show()

#statistics
#mw

# Mann-Whitney U test
statistic, p_value = ks_2samp(grouped_pdb['Std_Distance'],grouped['Std_Distance'])

# Eredmények kiíratása
print("Statisztika_mw_af_pdb:", statistic)
print("P-érték_mw:", p_value)

statistic, p_value = ks_2samp(grouped_full['Std_Distance'],grouped['Std_Distance'])

# Eredmények kiíratása
print("Statisztika_mw_af:", statistic)
print("P-érték_mw:", p_value)


#kruskal

# Kruskal-Wallis H próba végrehajtása
statistic, p_value = kruskal(grouped['Std_Distance'], grouped_pdb['Std_Distance'],grouped_full['Std_Distance'])

# Eredmények kiíratása
print("Statisztika_kw:", statistic)
print("P-érték_kw:", p_value)
