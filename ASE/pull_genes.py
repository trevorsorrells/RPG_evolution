#!/usr/bin/env python
################################################################################
# pull_genes.py
#
# creates a gene file for bowtie alignment from each of the genes annotated by YGAP
#
################################################################################

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

project_dir = '/Users/trevorsorrells/Documents/UCSF/Projects/Expression/allele_specific_expression/genomes_9_21_2015/'

def grab_seqs(genome_file, gff_file, output_name):
    #uses gff files and genome files for individual species in genomes_9_21 folder
    chromosomes = []
    for seq_record in SeqIO.parse(project_dir + genome_file, 'fasta'):
        chromosomes.append(seq_record)
    print chromosomes
    gene_seqs = []
    gene_names = []
    inFile = open(project_dir + gff_file, 'r')
    #add first exon because needs to test for multiple exons
    line_split = inFile.readline().split('\t')
    gene_name = line_split[8].split('"')[1]
    gene_names.append(gene_name)
    for chrom in chromosomes:
        if chrom.id == line_split[0]:
            if line_split[6] == '+':
                gene_seqs.append(chrom[int(line_split[3])-1:int(line_split[4])])
            else:
                gene_seqs.append(chrom[int(line_split[3])-1:int(line_split[4])].reverse_complement())
            break
    for line in inFile:
        line_split = line.split('\t')
        gene_name = line_split[8].split('"')[1]
        for chrom in chromosomes:
            if chrom.id == line_split[0]:
                if gene_names[-1] != gene_name:
                    gene_names.append(gene_name)
                    if line_split[6] == '+':
                        gene_seqs.append(chrom[int(line_split[3])-1:int(line_split[4])])
                    else:
                        gene_seqs.append(chrom[int(line_split[3])-1:int(line_split[4])].reverse_complement())
                else:
                    if line_split[6] == '+':
                        gene_seqs[-1] = gene_seqs[-1] + chrom[int(line_split[3])-1:int(line_split[4])]
                    else:
                       gene_seqs[-1] = chrom[int(line_split[3])-1:int(line_split[4])].reverse_complement() + gene_seqs[-1]
                break
    for seq, name in zip(gene_seqs, gene_names):
        seq.id = name
        seq.description = ''
    print len(gene_seqs)
    protein_seqs = []
    for gene in gene_seqs:
        protein_seqs.append(SeqRecord(gene.seq.translate(), id = gene.id, description=gene.description))
    #SeqIO.write(gene_seqs, project_dir + output_name, 'fasta')
    SeqIO.write(protein_seqs, project_dir + output_name + 'proteins.txt', 'fasta')

grab_seqs('Klactis_sequence.fsa.txt','Klac.gff.txt','Klac_genes.fsa.txt')
grab_seqs('kwic.genome.ygap.txt', 'Kwic_fixed.gff.txt', 'Kwic_genes.fsa.txt')
grab_seqs('km02.genome.ygap.txt','Kmar.gff.txt','Kmar_genes.fsa.txt')

grab_seqs('Klac_Kwic.fsa.txt','Klac_Kwic.gff.txt','KlacKwic_genes.txt')
grab_seqs('Klac_Kmar.fsa.txt','Klac_Kmar.gff.txt','KlacKmar_genes.txt')