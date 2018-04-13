#!/usr/bin/env python
################################################################################
# aln_reads.py
#
# Aligns sequencing reads to single species or hybrid genome
# Navigate to directory containing folders corresponding to sequencing lanes, run on server
#
################################################################################

import subprocess

all_files = ['TS03-ATCACG_S1_L001_R1_001.fastq.gz','TS03-CGATGT_S2_L001_R1_001.fastq.gz','TS03-TTAGGC_S3_L001_R1_001.fastq.gz','TS03-TGACCA_S4_L001_R1_001.fastq.gz','TS03-ACAGTG_S5_L001_R1_001.fastq.gz','TS03-GCCAAT_S6_L001_R1_001.fastq.gz','TS03-CAGATC_S7_L001_R1_001.fastq.gz','TS03-ACTTGA_S8_L001_R1_001.fastq.gz','TS03-GATCAG_S9_L001_R1_001.fastq.gz','TS03-TAGCTT_S10_L001_R1_001.fastq.gz','TS03-GGCTAC_S11_L001_R1_001.fastq.gz','TS03-CTTGTA_S12_L001_R1_001.fastq.gz']

aln_to_Kmar = ['TS03-GATCAG_S9_L001_R1_001.fastq.gz','TS03-TAGCTT_S10_L001_R1_001.fastq.gz','TS03-GGCTAC_S11_L001_R1_001.fastq.gz']
aln_to_Kwic = ['TS03-CAGATC_S7_L001_R1_001.fastq.gz','TS03-CTTGTA_S12_L001_R1_001.fastq.gz']
aln_to_Kaes = ['TS03-ACTTGA_S8_L001_R1_001.fastq.gz','TS04-TGACCA_S16_L002_R1_001.fastq.gz','TS04-ACAGTG_S17_L002_R1_001.fastq.gz','TS04-GCCAAT_S18_L002_R1_001.fastq.gz']

aln_to_KlacKwic = ['TS03-CAGATC_S7_L001_R1_001.fastq.gz','TS03-ACTTGA_S8_L001_R1_001.fastq.gz','TS03-CTTGTA_S12_L001_R1_001.fastq.gz']
aln_to_KlacKmar = ['TS03-GATCAG_S9_L001_R1_001.fastq.gz','TS03-TAGCTT_S10_L001_R1_001.fastq.gz','TS03-GGCTAC_S11_L001_R1_001.fastq.gz']

RNA_KlacKmar = ['TS04-CAGATC_S19_L002_R1_001.fastq.gz','TS04-ACTTGA_S20_L002_R1_001.fastq.gz','TS04-GATCAG_S21_L002_R1_001.fastq.gz']
RNA_KlacKwic = ['TS04-ATCACG_S13_L002_R1_001.fastq.gz','TS04-CGATGT_S14_L002_R1_001.fastq.gz','TS04-TTAGGC_S15_L002_R1_001.fastq.gz','TS04-TGACCA_S16_L002_R1_001.fastq.gz','TS04-ACAGTG_S17_L002_R1_001.fastq.gz','TS04-GCCAAT_S18_L002_R1_001.fastq.gz','TS04-TAGCTT_S22_L002_R1_001.fastq.gz','TS04-GGCTAC_S23_L002_R1_001.fastq.gz','TS04-CTTGTA_S24_L002_R1_001.fastq.gz']


def bowtieClineNoMismatch(read_file, genome_name):
    #no mismatches
    run_name = read_file[:4]
    sample_name = read_file.split('_')[1]
    sam_file = run_name + '_' + sample_name + '_' + genome_name +'.sam'
    details_file = 'details_' + run_name + '_' + sample_name + '_' + genome_name + '.txt'
    return 'bowtie2 -x /Users/Shared/sequencing_analysis/indexes/' + genome_name + ' --score-min L,0,0 -U ' + run_name + '/' + read_file + ' -S ' + sam_file + ' 2> ' + details_file

def bowtieClineDefault(read_file, genome_name):
    #default mismatches
    run_name = read_file[:4]
    sample_name = read_file.split('_')[1]
    sam_file = run_name + '_' + sample_name + '_' + genome_name +'.sam'
    details_file = 'details_' + run_name + '_' + sample_name + '_' + genome_name + '.txt'
    return 'bowtie2 -x /Users/Shared/sequencing_analysis/indexes/' + genome_name + ' -U ' + run_name + '/' + read_file + ' -S ' + sam_file + ' 2> ' + details_file

for read_file in aln_to_KlacKmar:
    cline = bowtieClineDefault(read_file, 'KlacKmar_genes')
    print cline
    return_code = subprocess.call(cline, shell=True)
    print return_code
    
for read_file in aln_to_KlacKwic:
    cline = bowtieClineDefault(read_file, 'KlacKwic_genes')
    print cline
    return_code = subprocess.call(cline, shell=True)
    print return_code
    
for read_file in RNA_KlacKmar:
    cline = bowtieClineDefault(read_file, 'KlacKmar_genes')
    print cline
    return_code = subprocess.call(cline, shell=True)
    print return_code
    
for read_file in RNA_KlacKwic:
    cline = bowtieClineDefault(read_file, 'KlacKwic_genes')
    print cline
    return_code = subprocess.call(cline, shell=True)
    print return_code