#!/bin/bash

download_dir="/bigdisk/users/kelvi/data/pdb_af"
not_found_file="/bigdisk/users/kelvi/not_found_cif.tsv"

# Loop through each line in codes.tsv
while IFS= read -r i; do
    # Check if the file exists
    if wget --spider "https://alphafold.ebi.ac.uk/files/AF-$i-F1-model_v4.cif" 2>&1 | grep -q '404 Not Found'; then
        # Check if the code is already in the not found file
        if ! grep -q "^$i$" "$not_found_file"; then
            echo -e "$i" >> "$not_found_file"
        fi
        echo "File not found for code $i. Skipping..."
    else
        # Download the file if it exists
        wget -P "$download_dir" "https://alphafold.ebi.ac.uk/files/AF-$i-F1-model_v4.cif"
    fi
done < /bigdisk/users/kelvi/codes.tsv

