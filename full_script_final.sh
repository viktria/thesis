#! /bin/bash

# Count how many proteins contain the same number of domains
cut -f 1 ipr_supradomains.dat | sort | uniq -c | sort -n | cut -c 1-8 | uniq -c > ipr_domains_per_protein.txt

# Transform it so that it contains uniprot_id, domain_1, domain_2
python3 one.py ipr_supradomains.dat > ipr_start_stop.txt

# Transformation of ipr_start_stop.txt to include only the uniprot_id and the two domain names
cut -d ',' -f 1,2,3 ipr_start_stop.txt > ipr_id_dom1_dom2.txt

# Joining data with pdb_id using the uniprot_id: filtered_pdb_chain_uniprot.tsv
# Initial files:
# ipr_start_stop.txt
# pdb_chain_uniprot.tsv
python3 2_pdb_filter.py

# Filtering ipr_supradomains_2_or_more.dat based on filtered_pdb_chain_uniprot.tsv
awk -F'\t' '{print $3 "," $4}' filtered_pdb_chain_uniprot_domain.tsv | sort | uniq > uniprot_filter_part1.txt
grep -f uniprot_filter_part1.txt ipr_id_dom1_dom2.txt > ipr_id_dom1_dom2_filtered.txt
awk -F'\t' '{print $3 "\t" $4}' filtered_pdb_chain_uniprot_domain.tsv | sort | uniq > uniprot_filter_part1.txt
awk -F'\t' '{print $3 "\t" $7}' filtered_pdb_chain_uniprot_domain.tsv | sort | uniq > uniprot_filter_part2.txt
sort uniprot_filter_part1.txt uniprot_filter_part2.txt | uniq > uniprot_filter.txt; rm uniprot_filter_*.txt
grep -f uniprot_filter.txt ipr_supradomains.dat > ipr_with_2_domains_filter.dat

# Filtering ipr_with_2_domains_filter.dat:
# ipr_with_2_domains_filter_02.dat filters the original file based on variance/average, removing domains with values worse than 0.2
python3 3_domain_filter.py

rm ipr_with_2_domains_filter.dat

# Filtering proteins appearing only once
cut -f 1 ipr_with_2_domains_filter_02.dat | sort | uniq -c | grep -v ' 1 ' | cut -c 9- > temp
grep -f temp ipr_with_2_domains_filter_02.dat > ipr_supradomains_filtered.dat; rm ipr_with_2_domains_filter_02.dat
{ head -n 1 filtered_pdb_chain_uniprot.tsv && grep -f temp filtered_pdb_chain_uniprot.tsv; } > refiltered_pdb_chain_uniprot.tsv ; rm filtered_pdb_chain_uniprot.tsv
grep -f temp uniprot_filter.txt > tmp

# Histogram of domain frequencies: histogram_single.png
python3 4_histograms.py

# Other processes
python3 5_filter.py
python3 6_unp_pdb_dom.py

# Calculating "center of mass" in the pdb data based on unp_pdb_dom.txt: distance_pdb.tsv
python3 7_count_pdb.py

# Selecting identifiers where "Not Found" is not the result
grep -v "Not Found" distance_pdb.tsv | awk '{print $1 " " $2 " " $3 " " $4}' | uniq -c | cut -c 9- > sel.txt
grep -f sel.txt unp_pdb_dom.txt > all_pdb.txt

# Selecting the best pdb_id-s associated with a given uniprot_id: best_pdb_full.txt based on all_pdb.txt
python3 8_gemmi.py

cut -d ' ' -f 1 best_pdb_full.txt| uniq -c | cut -c 9- > codes.tsv
cut -d ' ' -f 1,2,3,4 best_pdb_full.txt | sort | uniq -c | cut -c 9- | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}'  > temp
grep -f temp distance_pdb.tsv > test

# Downloading AF cif files corresponding to codes.tsv and writing any not found to not_found_cif.tsv
bash 9_download_script.sh

# Calculating "center of mass" in the AF data: distance_af.tsv
python3 10_count_af.py

# Calculating average and standard deviation for pdb data - domain pairs
# based on best_pdb_full.txt filtered using distance_pdb.tsv and not_found_cif.tsv
# resulting in grouped_stats_pdb.tsv
python3 11_mean_sdev_pdb.py

# Calculating average and standard deviation for AF data - domain pairs
# based on best_pdb_full.txt filtered using distance_af.tsv
# resulting in grouped_stats_af.tsv
python3 12_mean_sdev_af.py

# Creating plots using grouped_stats_pdb.tsv and grouped_stats_af.tsv
# as well as calculating distances between domain pairs for further analysis
# b_pdb_s_af_distances.tsv - PDB with variance greater than 4 but AF less than 2.5
# combined_std_histogram.png - standard deviation for both datasets
# std_compare.png - comparison of distances between datasets
python3 13_plot_std.py

# Calculating "center of mass" in the AF data using ipr_start_stop.txt: distance_af_full.tsv
python3 AF_full_run.py

# Merging intermediate results
python3 concat_tsv.py

# Comparing full AF with filtered AF: combined_std_histogram_af.png
python3 AF_full_plot_std.py
