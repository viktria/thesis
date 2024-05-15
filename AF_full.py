#!/usr/bin/python3
import os
import gemmi
import math
import csv
import pandas as pd
import sys
import time

def select_atoms(folder, unp, d1_name, domain1_start, domain1_end, d2_name, domain2_start, domain2_end):
    domain1_mean = gemmi.Vec3(0.0, 0.0, 0.0)
    domain2_mean = gemmi.Vec3(0.0, 0.0, 0.0)
    domain1_count = 0
    domain2_count = 0
    for pdb_id_file in os.listdir(folder):
        if pdb_id_file.startswith("AF-" + unp):
            mmCIF_file_path = os.path.join(folder, pdb_id_file)
            print(f"Reading file: {mmCIF_file_path}")
            try:
                doc = gemmi.cif.read(mmCIF_file_path)
            except FileNotFoundError:
                print(f"File not found: {mmCIF_file_path}")
                continue  # Continue to the next iteration if file not found
            for row in doc.sole_block().find("_atom_site.",
                                              ["group_PDB", "label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
                                              "label_seq_id"]):
                start_stop_str = row[5]
                if start_stop_str.isdigit():
                    start_stop = int(start_stop_str)
                    if (row[0] == "ATOM" and domain1_start <= start_stop <= domain1_end):
                            domain1_mean.x += float(row[2])
                            domain1_mean.y += float(row[3])
                            domain1_mean.z += float(row[4])
                            domain1_count += 1
                    elif (row[0] == "ATOM" and domain2_start <= start_stop <= domain2_end):
                            domain2_mean.x += float(row[2])
                            domain2_mean.y += float(row[3])
                            domain2_mean.z += float(row[4])
                            domain2_count += 1
                else:
                     # Handle cases where start_stop is not a valid integer
                    continue 
    if domain1_count == 0 or domain2_count == 0:
        return None  # Return None if no atoms were found for either domain
    domain1_mean /= domain1_count
    domain2_mean /= domain2_count
    
    # Calculate the distance between mean coordinates of two domains
    distance = math.sqrt((domain1_mean.x - domain2_mean.x)**2 + (domain1_mean.y - domain2_mean.y)**2 + (domain1_mean.z - domain2_mean.z)**2)
    
    return [d1_name,d2_name,distance]

start_time = time.time()
input_file = sys.argv[1]
directory, filename = os.path.split(input_file)
filename_without_extension, extension = os.path.splitext(filename)
prefix = filename_without_extension.split('_')[0]
pdb_folder = '/bigdisk/databases/AlphaFold/data/swissprot_cif/'
output_file = f'/bigdisk/users/kelvi/results/AF_full/distance_{prefix}_af.tsv'
o2 = f'/bigdisk/users/kelvi/results/AF_full/not_found/not_found_id_af_{prefix}.txt'

# Read input file using pandas
data = pd.read_csv(input_file, sep=',', header=None, index_col=False,
                        names=['SP_PRIMARY','Domain_1','Domain_2', 'Start_1','Stop_1','Start_2','Stop_2'])

# Open and write to output file
with open(output_file, 'a', newline='') as f_output, open(o2, 'a') as f_not_found:
    writer = csv.writer(f_output, delimiter='\t')
    writer.writerow(['UniProt ID','PDB ID','CHAIN','Domain_1','Domain_2', 'Distance','Merge_ID'])
    
    for _, row in data.iterrows():
        unp = row['SP_PRIMARY']
        domain1_name = row['Domain_1']
        domain1_start = row['Start_1']
        domain1_end = row['Stop_1']
        domain2_name = row['Domain_2']
        domain2_start = row['Start_2']
        domain2_end = row['Stop_2']
        
        # Call select_atoms function for each PDB ID and domain ranges
        result = select_atoms(pdb_folder, unp, domain1_name, domain1_start, domain1_end, domain2_name, domain2_start, domain2_end)
        
        if result is not None and result != "two"and result != "one" :
            writer.writerow([unp,result[0],result[1],result[2],domain1_start])
        else:
           #if everything went well, only those which don't have a .cif file downloaded
           f_not_found.write(unp + '\t' +  '\t' + domain1_name + '\t' + domain2_name + '\n')  # Write the PDB ID to the "not_found_id.txt" file 

end_time = time.time()
runtime = end_time - start_time
print("Script runtime:", runtime, "seconds")
print("Results written to:", output_file)
print("IDs not found written to:", o2) 
