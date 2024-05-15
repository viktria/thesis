#!/usr/bin/python3

import pandas as pd
from glob import glob

filename = '/bigdisk/users/kelvi/distance_full_af.tsv'

with open(filename, 'a') as singleFile:
    first_tsv = True
    for tsv in glob('/bigdisk/users/kelvi/results/AF_full/*.tsv'):
        if tsv == filename:
            pass
        else:
            header = True
            for line in open(tsv, 'r'):
                if first_tsv and header:
                    singleFile.write(line)
                    first_tsv = False
                    header = False
                elif header:
                    header = False
                else:
                    singleFile.write(line)
    singleFile.close()

df = pd.read_csv('distance_full_af.tsv', sep="\t",skiprows=1, 
                 names=['UniProt_ID', 'Domain_1', 'Domain_2', 'Distance', 'Merge_ID','empty'])
df.drop(columns=['Merge_ID','empty'], inplace=True)
df.to_csv('distance_full_af.tsv', sep='\t', index=False)


