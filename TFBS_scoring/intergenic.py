#!/usr/bin/env python
################################################################################
# intergenic.py
#
# Makes files of the intergenic of each gene from YGOB, CGOB, and JGI genome repositories
#
################################################################################

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import string
from operator import itemgetter
import gzip
import subprocess

YGOB_species = ['Cglabrata','Ecymbalariae','Egossypii','Kafricana','Klactis','Knaganishii','Lkluyveri','Lthermotolerans','Lwaltii','Ncastellii','Ndairenensis','Scerevisiae','Skudriavzevii','Smikatae','Suvarum','Tblattae','Tdelbrueckii','Tphaffii','Vpolyspora','Zrouxii']
CGOB_species = ['Calb','Cdub','Clus','Cort','Cpar','Cten','Ctro','Dhan','Lelo','Mgui','Spas','Ssti','Wo-1']
JGI_species = ['Arub','Bino','Cara','Ccas','Ctan','Cten','Dbru','Hpol','Hval','Hbur','Lsta','Mbic','Nful','Ptan','Pmem','Spas','Wano']
JGI_folder_names = ['Ascru1','Babin1','Canar1','Canca1','Canta1','Cante1','Dekbr2','Hanpo2','Hanva1_1','Hypbu1','Lipst1_1','Metbi1','Nadfu1','Pacta1_2','Picme2','Spapa3','Wican1']

#source folder for genome files
genome_info = '/Users/trevorsorrells/Documents/UCSF/Projects/genome_info/' 

def removeNNs(sequence):
    #truncates intergenic regions at the first two consecutive Ns
    for i in range(len(sequence)-2):
        if sequence[len(sequence)-i-1] == 'N' and sequence[len(sequence)-i-2] == 'N':
            return sequence[len(sequence)-i:]
    return sequence

def intergenic_CGOB(species):
    #extracts intergenic regions from CGOB sourced genomes
    tab_file = open(genome_info + 'genomes_cgob/gene_order/' + species + '.tab.txt','r')
    tabs = []
    [tabs.append(line.split('\t')) for line in tab_file]
    tab_file.close()
    for chromosome in SeqIO.parse(genome_info + 'genomes_cgob/genomes/' + species + '.fsa.txt','fasta'):
        chrom_num = chromosome.id.split('_')[2]
        last_end = 0
        for i, gene in enumerate(tabs):
            if gene[5] == chrom_num:
                if gene[1] == '1':
                    stop = int(gene[2])-1
                    start = last_end
                    gene.append(removeNNs(chromosome.seq[start:stop]))
                else:
                    start = int(gene[3])
                    if i == len(tabs)-1:
                        stop = len(chromosome)
                    elif tabs[i+1][5] != chrom_num:
                        stop = len(chromosome)
                    else:
                        stop = int(tabs[i+1][2])-1
                    gene.append(removeNNs(chromosome.seq[start:stop].reverse_complement()))
                last_end = int(gene[3])-1
    seq_records = []
    for gene in tabs:
        seq_records.append(SeqRecord(gene[-1], id=gene[0], description=''))
    SeqIO.write(seq_records, genome_info+ 'intergenics/' + species[:4] + '_intergenics.faa','fasta')

def intergenic_YGOB(species):
    #extracts intergenic regions from YGOB sourced genomes
    tab_file = open(genome_info + 'genomes_ygob/' + species + '_genome.tab.txt','r')
    tabs = []
    [tabs.append(line.split('\t')) for line in tab_file]
    tab_file.close()
    for chromosome in SeqIO.parse(genome_info + 'genomes_ygob/' + species + '_sequence.fsa.txt','fasta'):
        chrom_num = chromosome.id.split('_')[2]
        last_end = 0
        for i, gene in enumerate(tabs):
            if gene[5] == chrom_num:
                if gene[1] == '1':
                    stop = int(gene[2])-1
                    start = last_end
                    gene.append(removeNNs(chromosome.seq[start:stop]))
                else:
                    start = int(gene[3])
                    if i == len(tabs)-1:
                        stop = len(chromosome)
                    elif tabs[i+1][5] != chrom_num:
                        stop = len(chromosome)
                    else:
                        stop = int(tabs[i+1][2])-1
                    gene.append(removeNNs(chromosome.seq[start:stop].reverse_complement()))
                last_end = int(gene[3])-1
    seq_records = []
    for gene in tabs:
        seq_records.append(SeqRecord(gene[9], id=gene[0], description=''))
    SeqIO.write(seq_records, genome_info+ 'intergenics/' + species[:4] + '_intergenics.faa','fasta')


def intergenic_JGI(species, clade_name, folder_name, gff_file_name, genome_file_name):
    #extracts intergenic regions from JGI sourced genomes
    tab_file = gzip.open(genome_info + 'genomes_jgi/' + clade_name + '/' + folder_name + '/' + gff_file_name, 'r')
    tabs = []
    [tabs.append(line.split('\t')) for line in tab_file]
    tab_file.close()
    #create a data structure by chromosome then by gene, containing all the gff lines for each gene
    chromosome_order = []
    tabs_ordered = []
    current_chromosome_id = tabs[0][0]
    current_chromosome_tabs = []
    current_gene_id = tabs[0][8].split('"')[1]
    current_gene_lines = []
    for line in tabs:
        if current_gene_id == line[8].split('"')[1]:
            current_gene_lines.append(line)
        else:
            current_chromosome_tabs.append(current_gene_lines)
            current_gene_id = line[8].split('"')[1]
            current_gene_lines = [line]
            if line[0] != current_chromosome_id:
                chromosome_order.append(current_chromosome_id)
                tabs_ordered.append(current_chromosome_tabs)
                current_chromosome_id = line[0]
                current_chromosome_tabs = []
    chromosome_order.append(current_chromosome_id)
    tabs_ordered.append(current_chromosome_tabs)
    #create a list containing the [GeneName, strand, start, stop] for each gene coding sequence in each chromosome
    tab_summaries = []
    [tab_summaries.append([]) for x in chromosome_order]
    for chromosome in SeqIO.parse(genome_info + 'genomes_jgi/' + clade_name + '/' + folder_name + '/' + genome_file_name, 'fasta'):
        if chromosome.id not in chromosome_order:
            continue
        chrom_index = chromosome_order.index(chromosome.id)
        for current_gene_lines in tabs_ordered[chrom_index]:
            CDSs = []
            for gene_line in current_gene_lines:
                if gene_line[2] == 'CDS':
                    CDSs.append(gene_line)
            if len(CDSs) == 1:
                current_gene_beginning = int(CDSs[0][3])-1
                current_gene_end = int(CDSs[0][4])
            else:
                first_exon = CDSs[0]
                last_exon = CDSs[0]
                i=1
                while i < len(CDSs):
                    if CDSs[i][8].split()[-1] < first_exon[8].split()[-1]:
                        first_exon = CDSs[i]
                    if CDSs[i][8].split()[-1] > last_exon[8].split()[-1]:
                        last_exon = CDSs[i]
                    i += 1
                if first_exon[6] == '+':
                    current_gene_beginning = int(first_exon[3])-1
                    current_gene_end = int(last_exon[4])
                else:
                    current_gene_beginning = int(last_exon[3])-1
                    current_gene_end = int(first_exon[4])
            tab_summaries[chrom_index].append([current_gene_lines[0][8].split('"')[1], current_gene_lines[0][6], int(current_gene_beginning), int(current_gene_end)])
    #list of intergenic regions
    intergenics_for_output = []
    for chromosome in SeqIO.parse(genome_info + 'genomes_jgi/' + clade_name + '/' + folder_name + '/' + genome_file_name, 'fasta'):
        if chromosome.id not in chromosome_order:
            continue
        chrom_index = chromosome_order.index(chromosome.id)
        last_gene_end = 0
        for i, current_gene_summary in enumerate(tab_summaries[chrom_index]):
            if current_gene_summary[1] == '+':
                start = last_gene_end
                stop = current_gene_summary[2]
                intergenics_for_output.append(SeqRecord(removeNNs(chromosome.seq[start:stop]), id=current_gene_summary[0], description=''))
            else:
                start = current_gene_summary[3]
                if i < len(tab_summaries[chrom_index])-1:
                    stop = tab_summaries[chrom_index][i+1][2]
                else:
                    stop = len(chromosome)
                current_seq = removeNNs(chromosome.seq[start:stop].reverse_complement())
                seq_record = SeqRecord(current_seq, id=current_gene_summary[0], description='')
                intergenics_for_output.append(seq_record)
            last_gene_end = start
    SeqIO.write(intergenics_for_output, genome_info+ 'intergenics/' + species[:4] + '_intergenics.faa','fasta')

def intergenic_gff(species, genome_folder_name, clade_name, folder_name, gff_file_name, genome_file_name):
    #extracts intergenic regions for genomes with a gff file
    tab_file = open(genome_info + genome_folder_name + '/' + clade_name + '/' + folder_name + '/' + gff_file_name, 'r')
    tabs = []
    [tabs.append(line.split('\t')) for line in tab_file]
    tab_file.close()
    print tabs[1]
    #create a data structure by chromosome then by gene, containing all the gff lines for each gene
    chromosome_order = []
    tabs_ordered = []
    current_chromosome_id = tabs[0][0]
    current_chromosome_tabs = []
    current_gene_id = tabs[0][8].split('"')[1]
    current_gene_lines = []
    for line in tabs:
        if line[2] != 'CDS':
            continue
        if current_gene_id == line[8].split('"')[1]:
            current_gene_lines.append(line)
        else:
            current_chromosome_tabs.append(current_gene_lines)
            current_gene_id = line[8].split('"')[1]
            current_gene_lines = [line]
            if line[0] != current_chromosome_id:
                chromosome_order.append(current_chromosome_id)
                tabs_ordered.append(current_chromosome_tabs)
                current_chromosome_id = line[0]
                current_chromosome_tabs = []
    chromosome_order.append(current_chromosome_id)
    tabs_ordered.append(current_chromosome_tabs)
    #create a list containing the [GeneName, strand, start, stop] for each gene coding sequence in each chromosome
    tab_summaries = []
    [tab_summaries.append([]) for x in chromosome_order]
    for chromosome in SeqIO.parse(genome_info + genome_folder_name + '/' + clade_name + '/' + folder_name + '/' + genome_file_name, 'fasta'):
        if chromosome.id not in chromosome_order:
            continue
        chrom_index = chromosome_order.index(chromosome.id)
        for current_gene_lines in tabs_ordered[chrom_index]:
            CDSs = []
            for gene_line in current_gene_lines:
                if gene_line[2] == 'CDS':
                    CDSs.append(gene_line)
            if len(CDSs) == 0:
                continue
            if len(CDSs) == 1:
                current_gene_beginning = int(CDSs[0][3])-1
                current_gene_end = int(CDSs[0][4])
            else:
                print tabs_ordered[chrom_index]
                first_exon = CDSs[0]
                last_exon = CDSs[0]
                i=1
                while i < len(CDSs):
                    if CDSs[i][8].split()[-1] < first_exon[8].split()[-1]:
                        first_exon = CDSs[i]
                    if CDSs[i][8].split()[-1] > last_exon[8].split()[-1]:
                        last_exon = CDSs[i]
                    i += 1
                if first_exon[6] == '+':
                    current_gene_beginning = int(first_exon[3])-1
                    current_gene_end = int(last_exon[4])
                else:
                    current_gene_beginning = int(last_exon[3])-1
                    current_gene_end = int(first_exon[4])
            tab_summaries[chrom_index].append([current_gene_lines[0][8].split('"')[1], current_gene_lines[0][6], int(current_gene_beginning), int(current_gene_end)])
    #list of intergenic regions
    print tab_summaries[:10]
    intergenics_for_output = []
    for chromosome in SeqIO.parse(genome_info + genome_folder_name + '/' + clade_name + '/' + folder_name + '/' + genome_file_name, 'fasta'):
        if chromosome.id not in chromosome_order:
            continue
        chrom_index = chromosome_order.index(chromosome.id)
        last_gene_end = 0
        for i, current_gene_summary in enumerate(tab_summaries[chrom_index]):
            if current_gene_summary[1] == '+':
                start = last_gene_end
                stop = current_gene_summary[2]
                intergenics_for_output.append(SeqRecord(removeNNs(chromosome.seq[start:stop]), id=current_gene_summary[0], description=''))
            else:
                start = current_gene_summary[3]
                if i < len(tab_summaries[chrom_index])-1:
                    stop = tab_summaries[chrom_index][i+1][2]
                else:
                    stop = len(chromosome)
                current_seq = removeNNs(chromosome.seq[start:stop].reverse_complement())
                seq_record = SeqRecord(current_seq, id=current_gene_summary[0], description='')
                intergenics_for_output.append(seq_record)
            last_gene_end = start
    SeqIO.write(intergenics_for_output, genome_info+ 'intergenics/' + species[:4] + '_intergenics.faa','fasta')


def intergenic_Kdob():
    #extracts intergenic regions for the Kluyveromyces dobzhanskii genome in embl format
    #first load contig and gene information
    gene_features_list = [] 
    contig_features_list = [] 
    gene_features = [] #order to be strand, start, stop, type, name
    contig_features = [] #order to be strand, start, stop, name
    info_file = open(genome_info + 'genomes_other/Kdob/gsASS001Lm.aug.trna.KL.rRNA2.embl','r')
    for line in info_file:
        line_split = line.split()
        features_to_save = ['CDS','tRNA','rRNA']
        if line_split[0] == 'FT':
            if line_split[1] in features_to_save and len(line_split) == 3:
                if line_split[2][0] == 'c':
                    gene_features.append('-')
                else:
                    gene_features.append('+')
                start_stop = line_split[2].strip()
                start_stop = start_stop.replace(',','..').split('..')
                if len(start_stop) < 2:
                    gene_features = gene_features[:-1]
                    continue
                #this stuff is for removing non-digits from the start & stops
                all = string.maketrans('','')
                nodigs = all.translate(all, string.digits)
                #a number of genes are split between introns so need to catch that:
                gene_features.extend([start_stop[0].translate(all, nodigs), start_stop[-1].translate(all, nodigs), line_split[1]])
                #include end of stop codon in CDSs
                if line_split[1]=='CDS':
                    if gene_features[0] == '+':
                        gene_features[2] = str(int(gene_features[2])+3)
                    else:
                        gene_features[1] = str(int(gene_features[1])-3)
            elif line_split[1] == 'misc_feature':
                if line_split[2][0] == 'c':
                    contig_features.append('-')
                else:
                    contig_features.append('+')
                start_stop = line_split[2].split('..')
                #this stuff is for removing non-digits from the start & stops
                all = string.maketrans('','')
                nodigs = all.translate(all, string.digits)
                contig_features.extend([start_stop[0].translate(all, nodigs), start_stop[1].translate(all, nodigs)])
            elif len(gene_features) == 4:
                if line_split[1][:10] == '/locus_tag':
                    name_split = line_split[1].split('"')
                    gene_features.append(name_split[1])
                    gene_features_list.append(gene_features)
                    gene_features = []
            elif len(contig_features) == 3:
                if line_split[1][:6] == '/label':
                    name_split = line_split[1].split('=')
                    contig_features.append(name_split[1])
                    contig_features_list.append(contig_features)
                    contig_features = []
    #gene features look like this: ['-', '13172', '13891', 'CDS', 'KLDOg8']
    #contig features look like this: ['-', '1', '161517', 'contig00028']
    #this is written assuming the contigs and genes all in order
    i = 0 #gene number
    j = 0 #contig number
    for contig in SeqIO.parse(genome_info + 'genomes_other/Kdob/Kdob_genome.faa.txt','fasta'):
        contig_start = int(contig_features_list[j][1])-1
        last_gene_end = contig_start
        while i < len(gene_features_list) and int(gene_features_list[i][1]) < int(contig_features_list[j][2]):
            if gene_features_list[i][0] == '+':
                intergenic_start = last_gene_end-contig_start
                intergenic_stop = int(gene_features_list[i][1])-1-contig_start
                if contig_features_list[j][0] == '+':
                    gene_features_list[i].append(removeNNs(contig.seq[intergenic_start:intergenic_stop]))
                else:
                    gene_features_list[i].append(removeNNs(contig.seq.reverse_complement()[intergenic_start:intergenic_stop]))
            else:
                intergenic_start = int(gene_features_list[i][2])-contig_start
                if i == len(gene_features_list)-1 or int(gene_features_list[i+1][1]) >= int(contig_features_list[j][2]):
                    intergenic_stop = len(contig)
                else:
                    intergenic_stop = int(gene_features_list[i+1][1])-1-contig_start
                if contig_features_list[j][0] == '+':
                    gene_features_list[i].append(removeNNs(contig.seq[intergenic_start:intergenic_stop].reverse_complement()))
                else:
                    gene_features_list[i].append(removeNNs(contig.seq.reverse_complement()[intergenic_start:intergenic_stop].reverse_complement()))
            last_gene_end = int(gene_features_list[i][2])
            i += 1
        j += 1
    intergenic_records = []
    for gene in gene_features_list:
        if len(gene) <6:
            print gene
        intergenic_records.append(SeqRecord(gene[5], id=gene[4], description=''))
    SeqIO.write(intergenic_records, genome_info+'intergenics/Kdob_intergenics.faa','fasta')

def intergenic_embl_format(species_name, genome_folder_name, species_folder_name, embl_file_name, scaffold_file_name ):
    #extracts intergenic regions for generic embl format genomes
    #first load contig and gene information
    gene_features_list = [] 
    contig_features_list = [] 
    gene_features = [] #order to be strand, start, stop, type, name, contig_number
    contig_features = [] #order to be strand, start, stop, name
    info_file = open(genome_info + genome_folder_name +'/'+ species_folder_name +'/'+ embl_file_name,'r')
    chr_number = 0
    for line in info_file:
        line_split = line.split()
        features_to_save = ['CDS','tRNA','rRNA']
        if line_split[0] == 'FT':
            if line_split[1] in features_to_save and len(line_split) == 3:
                if line_split[2][0] == 'c':
                    gene_features.append('-')
                else:
                    gene_features.append('+')
                start_stop = line_split[2].strip()
                start_stop = start_stop.replace(',','..').split('..')
                if len(start_stop) < 2:
                    gene_features = gene_features[:-1]
                    continue
                #this stuff is for removing non-digits from the start & stops
                all = string.maketrans('','')
                nodigs = all.translate(all, string.digits)
                #a number of genes are split between introns so need to catch that:
                gene_features.extend([start_stop[0].translate(all, nodigs), start_stop[-1].translate(all, nodigs), line_split[1]])
                #include end of stop codon in CDSs
                if line_split[1]=='CDS':
                    if gene_features[0] == '+':
                        gene_features[2] = str(int(gene_features[2])+3)
                    else:
                        gene_features[1] = str(int(gene_features[1])-3)
            elif line_split[1] == 'source':
                chr_number += 1
                if line_split[2][0] == 'c':
                    contig_features.append('-')
                else:
                    contig_features.append('+')
                start_stop = line_split[2].split('..')
                #this stuff is for removing non-digits from the start & stops
                all = string.maketrans('','')
                nodigs = all.translate(all, string.digits)
                contig_features.extend([start_stop[0].translate(all, nodigs), start_stop[1].translate(all, nodigs)])
            elif len(gene_features) == 4:
                if line_split[1][:10] == '/locus_tag':
                    name_split = line_split[1].split('"')
                    gene_features.extend([name_split[1].split('_')[1], chr_number])
                    gene_features_list.append(gene_features)
                    gene_features = []
            elif len(contig_features) == 3:
                if line_split[1][:5] == '/note':
                    contig_features.append(line_split[2].strip('"'))
                    contig_features_list.append(contig_features)
                    contig_features = []
    #gene features look like this: ['-', '13172', '13891', 'CDS', 'KLDOg8',1]
    #contig features look like this: ['-', '1', '161517', 'contig00028',1]
    #this is written assuming the contigs and genes all in order
    i = 0 #gene number
    j = 0 #contig number
    print gene_features_list[:10]
    print contig_features_list[:10]
    for contig in SeqIO.parse(genome_info + genome_folder_name + '/' + species_folder_name + '/' + scaffold_file_name,'fasta'):
        contig_start = int(contig_features_list[j][1])-1
        last_gene_end = contig_start
        while i < len(gene_features_list) and int(gene_features_list[i][1]) < int(contig_features_list[j][2]):
            if gene_features_list[i][5] != j + 1:
                print gene_features_list[i]
                print j
                break 
            if gene_features_list[i][0] == '+':
                intergenic_start = last_gene_end-contig_start
                intergenic_stop = int(gene_features_list[i][1])-1-contig_start
                if contig_features_list[j][0] == '+':
                    gene_features_list[i].append(removeNNs(contig.seq[intergenic_start:intergenic_stop]))
                else:
                    gene_features_list[i].append(removeNNs(contig.seq.reverse_complement()[intergenic_start:intergenic_stop]))
            else:
                intergenic_start = int(gene_features_list[i][2])-contig_start
                if i == len(gene_features_list)-1 or int(gene_features_list[i+1][1]) >= int(contig_features_list[j][2]):
                    intergenic_stop = len(contig)
                else:
                    intergenic_stop = int(gene_features_list[i+1][1])-1-contig_start
                if contig_features_list[j][0] == '+':
                    gene_features_list[i].append(removeNNs(contig.seq[intergenic_start:intergenic_stop].reverse_complement()))
                else:
                    gene_features_list[i].append(removeNNs(contig.seq.reverse_complement()[intergenic_start:intergenic_stop].reverse_complement()))
            last_gene_end = int(gene_features_list[i][2])
            i += 1
        j += 1
    intergenic_records = []
    for gene in gene_features_list:
        if len(gene) <7:
            print gene
        intergenic_records.append(SeqRecord(gene[6], id=gene[4], description=''))
    SeqIO.write(intergenic_records, genome_info+'intergenics/' + species_name + '_intergenics.faa','fasta')

def intergenic_KaesKwic(species):
    #extracts intergenic regions from the Kluyveromyces aestuarii and wickerhamii genomes
    tab_file = open(genome_info+'Kwic_Kaes/'+species+'/all_orfs.ffn','r')
    gene_features = []
    gene_features_list = [] #order is orf_name, contig_number, strand, start, stop
    for line in tab_file:
        if line[0] == '>':
            line_split = line[1:].split('_')
            gene_features.extend([line[1:].strip(),line_split[0][:5],line_split[1]])
            indices = line_split[2].split('-')
            if gene_features[2] == '+':
                gene_features.extend([int(indices[0]),int(indices[1])])
            else:
                gene_features.extend([int(indices[1]),int(indices[0])])
            gene_features_list.append(gene_features)
            gene_features = []
    #gene features look like this ['00001.1_+_2974-3763_[K_aestuarii]', '00001', '+',2974,3763]
    #sort by start index, then by contig
    gene_features_list = sorted(gene_features_list, key=itemgetter(3))
    gene_features_list = sorted(gene_features_list, key=itemgetter(1))
    for i in range(20):
        print gene_features_list[i]

    i=0 #index of current gene
    for contig in SeqIO.parse(genome_info+'Kwic_Kaes/'+species+'/all_nt_sequences.fna','fasta'):
        contig_number = contig.id[6:11]
        while i<len(gene_features_list) and gene_features_list[i][1] == contig_number:
            if gene_features_list[i][2] == '+':
                #find intergenic start: only if previous element is completely before current element
                j=-1
                intergenic_start = 0
                while i+j >= 0 and gene_features_list[i+j][1] == contig_number and intergenic_start == 0:
                    if gene_features_list[i+j][4] < gene_features_list[i][3]:
                        if gene_features_list[i+j][2] == '+':
                            intergenic_start = gene_features_list[i+j][4]+3
                        else:
                            intergenic_start = gene_features_list[i+j][4]
                    else:
                        j += -1
                intergenic_stop = gene_features_list[i][3]
                gene_features_list[i].append(removeNNs(contig.seq[intergenic_start:intergenic_stop]))
            else:
                intergenic_start = gene_features_list[i][4]
                #find intergenic stop: only if previous element is completely after current element
                j=1
                intergenic_stop = len(contig)
                while i+j < len(gene_features_list) and gene_features_list[i+j][1] == contig_number and intergenic_stop == len(contig):
                    if gene_features_list[i+j][3] > gene_features_list[i][4]:
                        if gene_features_list[i+j][2] == '+':
                            intergenic_stop = gene_features_list[i+j][3]
                        else:
                            intergenic_stop = gene_features_list[i+j][3]-3
                    else:
                        j += 1
                gene_features_list[i].append(removeNNs(contig.seq[intergenic_start:intergenic_stop].reverse_complement()))
            last_gene_end = gene_features_list[i][3]
            i += 1
    intergenic_records = []
    for gene in gene_features_list:
        intergenic_records.append(SeqRecord(gene[5], id=gene[0],description=''))
    SeqIO.write(intergenic_records, genome_info+'intergenics/'+species+'_intergenics.faa','fasta')


def gunzip_file(species_to_add, genome_source_folder,clade_folder_name, species_folder_name, proteins_file_name):
    cline1 = 'gunzip ' + genome_info + genome_source_folder+'/' +clade_folder_name+'/'+species_folder_name+'/'+proteins_file_name
    print subprocess.call(cline1, shell=True)


def fix_Huva_chrom_names():
    #alter chromosome names so they can be used to extract intergenic regions
    new_chromosomes = []
    for chromosome in SeqIO.parse(genome_info + '/genomes_other/Huva/APLS01000001-APLS01000335.fasta.txt','fasta'):
        chromosome.id = chromosome.description.split()[5].strip(',')
        chromosome.description = ''
        new_chromosomes.append(chromosome)
    SeqIO.write(new_chromosomes, genome_info + '/genomes_other/Huva/fixed_chromosomes.fasta.txt','fasta')


for specie, folder in zip(JGI_species, JGI_folder_names):
    intergenic_JGI(specie, folder)

for specie in CGOB_species:
    intergenic_CGOB(specie)

for specie in YGOB_species:
    intergenic_YGOB(specie)

intergenic_Kdob()

intergenic_KaesKwic('Kwic')

fix_Huva_chrom_names()
intergenic_gff('Huva', 'genomes_other', '', 'Huva', 'Huva.gff.txt', 'fixed_chromosomes.fasta.txt')

