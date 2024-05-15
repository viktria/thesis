#!/usr/bin/python3

import matplotlib.pyplot as plt
import pandas as pd

# histogram_single
df1 = pd.read_csv('ipr_supradomains.dat', sep='\t', header=None)
df1.columns = ['A', 'B', 'C', 'D', 'E', 'F']

# Count occurrences of domain values
df2 = df1.A.value_counts().reset_index(name='Sum_of_domains')

# Group by the sum of domains
protein_dom_num = df2.groupby('Sum_of_domains').size()
plt.figure(figsize=(8, 6))
protein_dom_num.plot(kind='bar',color='mediumturquoise')
plt.yscale('log')
plt.title('Number of Proteins by Domain Count', fontsize=14)
plt.xlabel('Domain Count', fontsize=12)
plt.ylabel('Number of Proteins', fontsize=12)
plt.tick_params(axis='x', labelsize=10)
plt.savefig('/bigdisk/users/kelvi/results/histogram_single.png', dpi=300)
plt.close()

# histogram_pair
df = pd.read_csv('ipr_id_dom1_dom2.txt', sep=",", header=None)
df.columns = ["a", "b", "c"]
df3 = df.groupby(['b', 'c']).size().reset_index(name='Sum')
pair_dom_number = df3.groupby('Sum').size()
max_value = pair_dom_number.max()
fig, axs = plt.subplots(1, 2, figsize=(8, 6), gridspec_kw={'width_ratios': [30, 1]}, facecolor='w')

# Plot data for the first part
pair_dom_number.iloc[:30].plot(kind='bar', ax=axs[0],color='cornflowerblue')
axs[0].set_ylabel('Number of Proteins',fontsize=12)
axs[0].tick_params(right=False, labelright=False)  # Hide the y ticks and labels on the right side of the first plot
axs[0].set_ylim(0, max_value)
axs[0].set_xlabel(' ')  

# Plot data for the second part
pair_dom_number.iloc[-1:].plot(kind='bar', ax=axs[1],color='cornflowerblue')
axs[1].tick_params(left=False, labelleft=False)  # Hide the y ticks and labels on the left side of the second plot
axs[1].set_ylabel(' ')  # Remove the y-label for the second plot
axs[1].set_ylim(0, max_value)  
axs[1].set_xlabel(' ')
axs[1].get_yaxis().set_visible(False)

# Hide the spines between ax and ax2
axs[0].spines['right'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[0].yaxis.tick_left()
axs[0].tick_params(labelright=False)
axs[1].yaxis.tick_right()
fig.text(0.5, 0.02, 'Domain Pair Count', ha='center', fontsize=12)
fig.suptitle('Number of Domain Pair Occurrences', y=0.95, fontsize=14)
plt.subplots_adjust(wspace=0.05)
plt.savefig('/bigdisk/users/kelvi/results/histogram_pair.png', dpi=300)
plt.show()
