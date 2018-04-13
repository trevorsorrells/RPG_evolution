#!/usr/bin/env python
################################################################################
# ASE.py
#
# modifies large files on the server containing alignment files
# writes orthology file based on YGAP annotations
# will use bowtie2 mapped reads to ORFs rather than transcript level estimates
#
################################################################################

import re

project_dir = '/Volumes/RAID/trevor/sequencing/allele_specific_expression/'
file_dir = '/Volumes/RAID/trevor/sequencing/allele_specific_expression/output_files/aln_to_orfs/'

KlacKwic_file_names = ['TS03_S7','TS03_S8','TS03_S12','TS04_S13','TS04_S14','TS04_S15','TS04_S16','TS04_S17','TS04_S18','TS04_S22','TS04_S23','TS04_S24']
KlacKmar_file_names = ['TS03_S9','TS03_S10','TS03_S11','TS04_S19','TS04_S20','TS04_S21']


def truncate_file(input_file_name, length):
    #truncates a sam file to do faster tests with methods
    outFile = open(file_dir + input_file_name + '_truncated.sam','w')
    with open(file_dir + input_file_name + '_KlacKwic_genes_unique.sam','r') as inFile:
        i = 0
        for line in inFile:
            if i < length:
                 outFile.write(line)
                 i += 1
            else:
                break
    outFile.close()


def count_reads(input_file_name, genome_file_path):
    #counts reads in each orf in genome_file_path
    #input_file_path is a sam file with multi-mapping reads removed
    gene_names = []
    gene_read_counts = []
    print input_file_name
    with open(genome_file_path,'r') as inFile:
        for line in inFile:
            if line[0] == '>':
                gene_names.append(line[1:].strip())
                gene_read_counts.append(0)
    with open(file_dir + input_file_name , 'r') as inFile:
        for line in inFile:
            line_split = line.split()
            for info in line_split:
                try:
                    if info[:4] in ['KLLA','KWIC','KM02']:
                        gene_read_counts[gene_names.index(info)] += 1
                except IndexError:
                    continue
    with open(file_dir + 'counts_' + input_file_name + '.txt' , 'w') as outFile:
        for name, count in zip(gene_names, gene_read_counts):
            outFile.write(name + '\t' + str(count) + '\n')

for file_name in KlacKwic_file_names:
    count_reads(file_name + '_KlacKwic_genes_unique.sam', project_dir + 'genomes_9_21_2015/KlacKwic_genes.txt')

for file_name in KlacKmar_file_names:
    count_reads(file_name + '_KlacKmar_genes_unique.sam', project_dir + 'genomes_9_21_2015/KlacKmar_genes.txt')