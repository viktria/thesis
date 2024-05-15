#!/usr/bin/python3
import sys
import csv

def print_pairs(p, pairs):
    for pair in pairs:
        print(p + "," + pair[0] + "," + pair[3] + "," + str(pair[1]) + "," + str(pair[2])+ "," + str(pair[4]) + "," + str(pair[5]))

def process_data(filename):
    protein_domain_map = {}
    with open(filename, "r") as data_file:
        tsv_reader = csv.reader(data_file, delimiter='\t')

        for row in tsv_reader:
            protein_id = row[0]
            domain_id = row[1]
            domain_start = int(row[4])
            domain_end = int(row[5])

            if protein_id not in protein_domain_map:
                protein_domain_map[protein_id] = []

            protein_domain_map[protein_id].append((domain_id, domain_start, domain_end))
    for protein_id, domains in protein_domain_map.items():
        pairs = []
        sorted_domains = sorted(domains, key=lambda x: x[1])# Sort domains based on start position
        for i in range(len(sorted_domains) - 1):
            pairs.append((sorted_domains[i][0],sorted_domains[i][1],sorted_domains[i][2] , sorted_domains[i + 1][0],sorted_domains[i + 1][1], sorted_domains[i + 1][2]))
        print_pairs(protein_id, pairs)

process_data(sys.argv[1])
