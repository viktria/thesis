#!/usr/bin/python3
import os
import gemmi
import math
import csv
import pandas as pd

def select_atoms(folder, pdb_id, chain_id, unp, d1_name, domain1_start, domain1_end, d2_name, domain2_start, domain2_end):
    subfolder_name = pdb_id[1:3]  # Extract the middle two characters of the PDB ID
    subfolder_path = os.path.join(folder, subfolder_name)
    if not os.path.isdir(subfolder_path):
        return None
    
    domain1_mean = gemmi.Vec3(0.0, 0.0, 0.0)
    domain2_mean = gemmi.Vec3(0.0, 0.0, 0.0)
    domain1_count = 0
    domain2_count = 0
    chid=chain_id
    for pdb_id_file in os.listdir(subfolder_path):
        if pdb_id_file.startswith(pdb_id):
            mmCIF_file_path = os.path.join(subfolder_path, pdb_id_file)
            print(f"Reading file: {mmCIF_file_path}")
            doc = gemmi.cif.read(mmCIF_file_path)
            for row in doc.sole_block().find("_atom_site.",
                                              ["group_PDB", "label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
                                               "pdbx_sifts_xref_db_name", "pdbx_sifts_xref_db_acc","pdbx_sifts_xref_db_num"]):
                start_stop_str = row[7]
                if start_stop_str.isdigit():
                    start_stop = int(start_stop_str)
                    if (row[0] == "ATOM" and row[5] == "UNP" and row[6] == unp and  
                        domain1_start <= start_stop <= domain1_end):
                        if chid == '-':
                            chid = row[1][0]
                        if chid == row[1][0]:
                            domain1_mean.x += float(row[2])
                            domain1_mean.y += float(row[3])
                            domain1_mean.z += float(row[4])
                            domain1_count += 1
                    if (row[0] == "ATOM" and row[5] == "UNP" and row[6] == unp and  
                            domain2_start <= start_stop <= domain2_end):
                        if chid == '-':
                            chid = row[1][0]
                        if chid == row[1][0]:
                            domain2_mean.x += float(row[2])
                            domain2_mean.y += float(row[3])
                            domain2_mean.z += float(row[4])
                            domain2_count += 1
                else:
                     # Handle cases where start_stop is not a valid integer
                    continue             
    if domain1_count == 0 or domain2_count == 0 :
        return None  # Return None if no atoms were found for either domain
    
    
    domain1_mean /= domain1_count
    domain2_mean /= domain2_count
    # Calculate the distance between mean coordinates of two domains
    distance = math.sqrt((domain1_mean.x - domain2_mean.x)**2 + (domain1_mean.y - domain2_mean.y)**2 + (domain1_mean.z - domain2_mean.z)**2)
    
    return [d1_name,d2_name,distance]

input_file = '/bigdisk/users/kelvi/unp_pdb_dom.txt'
pdb_folder = '/home/data/wwPDB/data/structures/divided/updated_mmcif/'
output_file = 'distance_pdb.tsv'
o2 = 'not_found_pdb.txt'

# Read input file using pandas
data = pd.read_csv(input_file, sep=' ')

# Open and write to output file
with open(output_file, 'w', newline='') as f_output, open(o2, 'w') as f_not_found:
    writer = csv.writer(f_output, delimiter='\t')
    writer.writerow(['UniProt ID','PDB ID','CHAIN','Domain_1','Domain_2', 'Distance','Merge_ID'])
    for _, row in data.iterrows():
        pdb_id = row['PDB']
        chain_id = row['CHAIN']
        unp = row['SP_PRIMARY']
        domain1_name = row['Domain_1']
        domain1_start = row['Start_1']
        domain1_end = row['Stop_1']
        domain2_name = row['Domain_2']
        domain2_start = row['Start_2']
        domain2_end = row['Stop_2']
        
        result = select_atoms(pdb_folder, pdb_id, chain_id, unp, domain1_name, int(domain1_start), int(domain1_end), domain2_name, int(domain2_start), int(domain2_end))
        
        if result is not None:
            writer.writerow([unp,pdb_id,chain_id,result[0],result[1],result[2],domain1_start])
        else:
            f_not_found.write(unp + '\t' + pdb_id + '\t' + domain1_name + '\t' + domain2_name + '\n')

print("Results written to:", output_file)