#!/usr/bin/python3

import os
import gemmi
import csv

def select_best_pdb_entries(input_file, pdb_folder):
    print("Selecting the best PDB ID for each UniProt ID based on atom count...")
    best_atom_count_found = 0
    entries_processed = 0
    best_pdb_entries = {}

    with open('best_pdb_full.txt', 'a') as output_file:
        with open(input_file, 'r') as f:
            for line in f:
                uniprot_id, pdb_id, chain, d1, start1, stop1, d2, start2, stop2 = line.strip().split()
                keys = [uniprot_id,d1,start1]
                best_pdb_entries.setdefault(tuple(keys), [])
                best_pdb_entries[uniprot_id,d1,start1].append((pdb_id, chain, d1, start1, stop1, d2, start2, stop2))

        for key, values in best_pdb_entries.items():
            uniprot_id = key[0]
            print(f"Processing UniProt ID: {uniprot_id}")
            pdb_ids = []
            chain = []
            d1 = []
            start1 = []
            stop1 = []
            d2 = [] 
            start2 = []
            stop2 = []
            for pdb_entry in values:
                pdb_id, chain_val, d1_val, start1_val, stop1_val, d2_val, start2_val, stop2_val = pdb_entry
                pdb_ids.append(pdb_id)
                chain.append(chain_val)
                d1.append(d1_val)
                start1.append(start1_val)
                stop1.append(stop1_val)
                d2.append(d2_val)
                start2.append(start2_val)
                stop2.append(stop2_val)
                # Check if only one PDB ID is associated with the UniProt ID
            if len(values) == 1:
                best_atom_count_found += 1  # Increment the count as we're adding this single PDB ID directly
                output_file.write(f"{uniprot_id} {pdb_ids[0]} {chain[0]} {d1[0]} {int(start1[0])} {int(stop1[0])} {d2[0]} {int(start2[0])} {int(stop2[0])}\n")
            else:
                processed_atom_counts, processed_entries = process_uniprot_entries(uniprot_id, pdb_ids, chain, d1, start1, stop1, d2, start2, stop2, pdb_folder)
                best_atom_count_found += processed_atom_counts
                entries_processed += processed_entries

        print(f"Found best atom counts for {best_atom_count_found}/{entries_processed} entries.")

def process_uniprot_entries(uniprot_id, pdb_ids, chain, d1, start1s, stop1s, d2, start2s, stop2s, pdb_folder):
    best_atom_count = 0
    best_pdb_id = None
    processed_atom_counts = 0
    processed_entries = 0

    with open('best_pdb_full.txt', 'a') as output_file:
        for i in range(len(pdb_ids)):
            pdb_id = pdb_ids[i]
            chin = chain[i]
            d_1 = d1[i]
            start1 = int(start1s[i])
            stop1 = int(stop1s[i])
            d_2 = d2[i]
            start2 = int(start2s[i])
            stop2 = int(stop2s[i])
            print(f"Processing PDB ID: {pdb_id}")
            subfolder_name = pdb_id[1:3]
            subfolder_path = os.path.join(pdb_folder, subfolder_name)
            if not os.path.isdir(subfolder_path):
                continue
            domain1_count = 0
            domain2_count = 0
            atom_count = 0
            for pdb_id_file in os.listdir(subfolder_path):
                if pdb_id_file.startswith(pdb_id):
                    mmCIF_file_path = os.path.join(subfolder_path, pdb_id_file)
                    print(f"Reading file: {mmCIF_file_path}")
                    doc = gemmi.cif.read(mmCIF_file_path)
                    atom_count = 0
                    domain1_count = 0
                    domain2_count = 0
                    for row in doc.sole_block().find("_atom_site.",
                                                     ["group_PDB", "label_asym_id","pdbx_sifts_xref_db_name",
                                                      "pdbx_sifts_xref_db_acc","pdbx_sifts_xref_db_num"]):
                        start_stop_str = row[4]
                        if start_stop_str.isdigit():
                            start_stop = int(start_stop_str)
                            if (row[0] == "ATOM" and row[2] == "UNP" and row[3] == uniprot_id and 
                                start1 <= start_stop <= stop1):
                                if chin == row[1][0]:
                                    domain1_count += 1
                                else:
                                    continue
                            else:
                                continue    
                            if(row[0] == "ATOM" and row[2] == "UNP" and row[3] == uniprot_id and 
                                start2 <= start_stop <= stop2):
                                if chin == row[1][0]:
                                    domain2_count += 1
                                else:
                                    continue
                            else:
                                continue
                        else:
                         # Handle cases where start_stop is not a valid integer
                            continue    
                    if domain1_count != 0 and domain2_count != 0:
                        atom_count = domain1_count + domain2_count
                        processed_entries += 1
                    else:
                        atom_count = 0       
            if atom_count > best_atom_count:
                best_atom_count = atom_count
                print(f"best_atom_count {best_atom_count}")
                best_pdb_id = pdb_id
                best_data = [chin,d_1,start1,stop1,d_2,start2,stop2]
                processed_atom_counts = 1
        if best_pdb_id is not None:
            output_file.write(f"{uniprot_id} {best_pdb_id} {best_data[0]} {best_data[1]} {best_data[2]} {best_data[3]} {best_data[4]} {best_data[5]} {best_data[6]}\n")

    return processed_atom_counts, processed_entries

input_file = '/bigdisk/users/kelvi/all_pdb.txt'
pdb_folder = '/home/data/wwPDB/data/structures/divided/updated_mmcif/'

select_best_pdb_entries(input_file, pdb_folder)
