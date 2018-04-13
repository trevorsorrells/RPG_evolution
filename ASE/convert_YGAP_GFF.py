#!/usr/bin/env python
################################################################################
# convert_YGAP_GFF.py
#
# Converts YGAP output to gff format, also outputs fasta file of genes
# Mates genomes in silico to form interspecies hybrids
# 
################################################################################

project_dir = '/Users/trevorsorrells/Documents/UCSF/Projects/Expression/allele_specific_expression/'

def test_chr_size(annotation_path, genome_path):
    #to test whether ygob file column 5 is chr, Scaffold, or order number
    inFile1 = open(genome_path, 'r')
    genome_scaffold_n = []
    genome_chr_n = []
    chr_len_list = []
    curr_chr = ''
    line1 = inFile1.readline()
    genome_chr_n.append(int(line1[4:].split()[0]))
    genome_scaffold_n.append(int(line1.split()[1][9:]))
    for line in inFile1:
        if line[0] == '>':
            genome_chr_n.append(int(line[4:].split()[0]))
            genome_scaffold_n.append(int(line.split()[1][9:]))
            chr_len_list.append(len(curr_chr))
            curr_chr = ''
        else: 
            curr_chr = curr_chr + line.strip()
    chr_len_list.append(len(curr_chr))
    inFile1.close()
    print len(genome_scaffold_n)
    print genome_scaffold_n
    print len(genome_chr_n)
    print genome_chr_n
    print len(chr_len_list)
    print chr_len_list
    
    order_bool = True
    chr_bool = True
    scaffold_bool = True
    
    inFile2 = open(annotation_path, 'r')
    for line in inFile2:
        line_split = line.split()
        col_5_number = int(line_split[5])
        coordinates = line_split[7]
        if 'complement' in coordinates:
            complement = True
            coordinates = coordinates[10:]
        else:
            complement = False
        coordinates = coordinates.split(',')
        for exon in coordinates:
            start_stop = exon.split('..')
            stop = int(start_stop[1].strip('()'))
            if stop > chr_len_list[col_5_number-1]:
                order_bool = False
            if stop > chr_len_list[genome_chr_n.index(col_5_number)]:
                chr_bool = False
            if stop > chr_len_list[genome_scaffold_n.index(col_5_number)]:
                scaffold_bool = False
    print 'order ' + str(order_bool)
    print 'chr ' + str(chr_bool)
    print 'scaffold ' + str(scaffold_bool)

def check_gff_with_genome(genome_path,annotation_path):
    #to test whether hybrid genome gff file matches genome chromosomes
    inFile1 = open(genome_path, 'r')
    genome_chr_n = []
    chr_len_list = []
    curr_chr = ''
    genome_chr_n.append(inFile1.readline()[1:].split()[0])
    for line in inFile1:
        if line[0] == '>':
            genome_chr_n.append(line[1:].split()[0])
            chr_len_list.append(len(curr_chr))
            curr_chr = ''
        else: 
            curr_chr = curr_chr + line.strip()
    chr_len_list.append(len(curr_chr))
    inFile1.close()
    print len(genome_chr_n)
    print genome_chr_n
    print len(chr_len_list)
    print chr_len_list

    inFile2 = open(annotation_path, 'r')
    for line in inFile2:
        line_split = line.split()
        len_chromosome = chr_len_list[genome_chr_n.index(line_split[0])]
        if int(line_split[4]) > len_chromosome:
            print line_split[9] + ' end ' + line_split[4] + ' longer than ' + line_split[0] + ' ' + str(len_chromosome)

def alter_gff_kwic(genome_path,annotation_path, gff_outpath):
    #YGAP annotated the ends of genes going beyond the length of a chromosome
    #this method fixes this
    inFile1 = open(genome_path, 'r')
    genome_chr_n = []
    chr_len_list = []
    curr_chr = ''
    genome_chr_n.append(inFile1.readline()[1:].split()[0])
    for line in inFile1:
        if line[0] == '>':
            genome_chr_n.append(line[1:].split()[0])
            chr_len_list.append(len(curr_chr))
            curr_chr = ''
        else: 
            curr_chr = curr_chr + line.strip()
    chr_len_list.append(len(curr_chr))
    inFile1.close()
    print len(genome_chr_n)
    print genome_chr_n
    print len(chr_len_list)
    print chr_len_list

    outFile = open(gff_outpath, 'w')
    inFile2 = open(annotation_path, 'r')
    for line in inFile2:
        line_split = line.split('\t')
        len_chromosome = chr_len_list[genome_chr_n.index(line_split[0])]
        if int(line_split[4]) > len_chromosome:
            print line_split[8] + ' end ' + line_split[4] + ' longer than ' + line_split[0] + ' ' + str(len_chromosome)
            line_split[4] = str(int(line_split[4])-3)
            outFile.write('\t'.join(line_split))
        else:
            outFile.write(line)
    inFile2.close()
    outFile.close()


def convert_to_gff(annotation_path, genome_path, output_path):
    #first have to get scaffold names dictionary for the genome
    #original genomes with orresponding old chr identifiers are located in genome_info folder
    inFile1 = open(genome_path, 'r')
    genome_chr_list = []
    genome_scaffold_list = []
    for line in inFile1:
        if line[0] == '>':
            genome_chr_list.append(line[1:].split()[0])
            genome_scaffold_list.append(line.split()[1].strip())
    inFile1.close()
    inFile2 = open(annotation_path, 'r')
    gff_lines = []
    for line in inFile2:
        line_split = line.split()
        coordinates = line_split[7]
        if 'complement' in coordinates:
            complement = True
            coordinates = coordinates[10:]
        else:
            complement = False
        coordinates = coordinates.split(',')
        for exon in coordinates:
            gff_list = ['chr' + line_split[5], 'YGOB/YGAP','CDS']
            start_stop = exon.split('..')
            gff_list.extend([start_stop[0].strip('()'), start_stop[1].strip('()'), '.'])
            if complement:
                gff_list.append('-')
            else:
                gff_list.append('+')
            gff_list.extend(['.', 'gene_id "' + line_split[0]+ '"; transcript_id "' + line_split[0]+ '"'])
            gff_lines.append('\t'.join(gff_list))
    inFile2.close()
    outFile = open(output_path, 'w')
    outFile.write('\n'.join(gff_lines))
    outFile.close()

def mate_genomes(genome_file1, genome_file2, gff_file1, gff_file2, genome_out, gff_out):
    #creates hybrid genomes from original genomes
    #adds n to all the chromosome names in the 2nd genome when creating the file
    #to keep unique chromosome identifiers (where n is the number of chromosomes
    #in genome1)
    inFile = open(project_dir + genome_file1, 'r')
    genomeOutFile = open(project_dir + genome_out, 'w')
    chromosome_n = 0
    for line in inFile:
        if line[0] == '>':
            chromosome_n += 1
            line_split = line.split()
            genomeOutFile.write(line_split[0]+ '\n')
        else:
            genomeOutFile.write(line)
    inFile.close()
    inFile = open(project_dir + genome_file2, 'r')
    for line in inFile:
        if line[0] == '>':
            line_split = line.split()
            chromosome = int(line_split[0][4:])+chromosome_n
            genomeOutFile.write(line_split[0][:4]+str(chromosome) + '\n')
        else:
            genomeOutFile.write(line)
    inFile.close()
    genomeOutFile.close()
    inFile = open(project_dir + gff_file1, 'r')
    gffOutFile = open(project_dir + gff_out, 'w')
    gffOutFile.write(inFile.read())
    gffOutFile.write('\n')
    inFile.close()
    inFile = open(project_dir + gff_file2, 'r')
    for line in inFile:
        line_split = line.split('\t')
        chromosome = int(line_split[0][3:])+chromosome_n
        line_split[0] = line_split[0][:3]+str(chromosome)
        gffOutFile.write('\t'.join(line_split))
    inFile.close()
    gffOutFile.close()


convert_to_gff(project_dir + 'genomes_9_21_2015/km02.annotation_final.ygap.txt', project_dir + 'genomes_9_21_2015/km02.genome.ygap.txt',project_dir + 'genomes_9_21_2015/Kmar.gff.txt')
convert_to_gff( project_dir + 'genomes_9_21_2015/kwic.annotation_final.ygap.txt', project_dir + 'genomes_9_21_2015/kwic.genome.ygap.txt',project_dir + 'genomes_9_21_2015/Kwic.gff.txt')
alter_gff_kwic(project_dir +'genomes_9_21_2015/kwic.genome.ygap.txt', project_dir + 'genomes_9_21_2015/Kwic.gff.txt',project_dir + 'genomes_9_21_2015/Kwic_fixed.gff.txt')
convert_to_gff(project_dir + 'genomes_9_21_2015/Klactis_genome.tab.txt', project_dir + 'genomes_9_21_2015/Klactis_sequence.fsa.txt',project_dir + 'genomes_9_21_2015/Klac.gff.txt')

mate_genomes('genomes_9_21_2015/Klactis_sequence.fsa.txt', 'genomes_9_21_2015/kwic.genome.ygap.txt', 'genomes_9_21_2015/Klac.gff.txt', 'genomes_9_21_2015/Kwic_fixed.gff.txt', 'genomes_9_21_2015/Klac_Kwic.fsa.txt', 'genomes_9_21_2015/Klac_Kwic.gff.txt')
mate_genomes('genomes_9_21_2015/Klactis_sequence.fsa.txt', 'genomes_9_21_2015/km02.genome.ygap.txt', 'genomes_9_21_2015/Klac.gff.txt', 'genomes_9_21_2015/Kmar.gff.txt', 'genomes_9_21_2015/Klac_Kmar.fsa.txt', 'genomes_9_21_2015/Klac_Kmar.gff.txt')

test_chr_size(project_dir + 'genomes_9_21_2015/km02.annotation_final.ygap.txt', project_dir + 'genomes_9_21_2015/km02.genome.ygap.txt')
test_chr_size(project_dir + 'genomes_9_21_2015/kwic.annotation_final.ygap.txt', project_dir + 'genomes_9_21_2015/kwic.genome.ygap.txt')

check_gff_with_genome(project_dir +'genomes_9_21_2015/Klactis_sequence.fsa.txt', project_dir + 'genomes_9_21_2015/Klac.gff.txt')
check_gff_with_genome(project_dir +'genomes_9_21_2015/km02.genome.ygap.txt', project_dir + 'genomes_9_21_2015/Kmar.gff.txt')
check_gff_with_genome(project_dir +'genomes_9_21_2015/kwic.genome.ygap.txt', project_dir + 'genomes_9_21_2015/Kwic_fixed.gff.txt')
check_gff_with_genome(project_dir +'genomes_9_21_2015/Klac_Kmar.fsa.txt', project_dir + 'genomes_9_21_2015/Klac_Kmar.gff.txt')
check_gff_with_genome(project_dir +'genomes_9_21_2015/Klac_Kwic.fsa.txt', project_dir + 'genomes_9_21_2015/Klac_Kwic.gff.txt')
