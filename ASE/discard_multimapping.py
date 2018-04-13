#!/usr/bin/env python
################################################################################
# discard_multimapping.py
#
# takes sam files and removes reads that map equally well to 2 or more locations
# (to avoid ambiguous mapping)
#
################################################################################


file_dir = '/Volumes/RAID/trevor/sequencing/allele_specific_expression/output_files/aln_to_orfs/'

KlacKwic_file_names = ['TS03_S7','TS03_S8','TS03_S12','TS04_S13','TS04_S14','TS04_S15','TS04_S16','TS04_S17','TS04_S18','TS04_S22','TS04_S23','TS04_S24']
KlacKmar_file_names = ['TS03_S9','TS03_S10','TS03_S11','TS04_S19','TS04_S20','TS04_S21']

for file_name in KlacKwic_file_names:
    saved_read_count = 0
    discarded_read_count = 0
    outFile = open(file_dir + file_name + '_KlacKwic_genes_unique.sam','w')
    with open(file_dir + file_name + '_KlacKwic_genes.sam','r') as inFile:
        for line in inFile:
            if line[0] == '@':
                 outFile.write(line)
                 continue
            line_split = line.split()
            qual_score = 300
            second_best = -100
            for info in line_split:
                try:
                    if info[:5] == 'AS:i:':
                        qual_score = int(info[5:])
                    elif info[:5] == 'XS:i:':
                        second_best = int(info[5:])
                except IndexError:
                    continue
            if qual_score == 300:
                continue
            if qual_score > second_best:
                outFile.write(line)
                saved_read_count += 1
            else:
                discarded_read_count += 1
    outFile.close()
    print file_name + ' KlacKwic saved ' + str(saved_read_count) + ' reads, discarded ' + str(discarded_read_count) + ' or ' + str(float(discarded_read_count)/(saved_read_count+discarded_read_count)*100) + '%'

for file_name in KlacKmar_file_names:
    saved_read_count = 0
    discarded_read_count = 0
    outFile = open(file_dir + file_name + '_KlacKmar_genes_unique.sam','w')
    with open(file_dir + file_name + '_KlacKmar_genes.sam','r') as inFile:
        for line in inFile:
            line_split = line.split()
            second_best = -100
            for info in line_split:
                try:
                    if info[:5] == 'AS:i:':
                        qual_score = int(info[5:])
                    elif info[:5] == 'XS:i:':
                        second_best = int(info[5:])
                except IndexError:
                    continue
            if qual_score > second_best:
                outFile.write(line)
                saved_read_count += 1
            else:
                discarded_read_count += 1
    outFile.close()
    print file_name + ' KlacKmar saved ' + str(saved_read_count) + ' reads, discarded ' + str(discarded_read_count) + ' or ' + str(float(discarded_read_count)/(saved_read_count+discarded_read_count)*100) + '%'
