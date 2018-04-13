#!/usr/bin/env python
################################################################################
# ASE.py
#
# looks for allele specific expression in various ways, can be used on local laptop computer
# writes orthology file based on YGAP annotations
# will use bowtie2 mapped reads to ORFs rather than transcript level estimates
#
################################################################################

import re
from scipy import stats
import math

project_dir = '/Users/trevorsorrells/Documents/UCSF/Projects/Expression/allele_specific_expression/'

KlacKwic_file_names = ['TS03_S7','TS03_S8','TS03_S12','TS04_S13','TS04_S14','TS04_S15','TS04_S16','TS04_S17','TS04_S18','TS04_S22','TS04_S23','TS04_S24']
#these read counts are after filtering out reads that map multiple times equally well
KlacKwic_read_counts = [15607048, 14203091, 20928491, 20170155, 22025225, 19099931, 21397747, 18298962, 30345095, 22426132, 23293393, 17734184]

KlacKmar_file_names = ['TS03_S9','TS03_S10','TS03_S11','TS04_S19','TS04_S20','TS04_S21']
KlacKmar_read_counts = [22817082, 30233079, 28136083, 29175477, 19652398, 22540374]

KlacKwic_gDNA = ['TS03_S7','TS03_S8','TS03_S12']
KlacKmar_gDNA = ['TS03_S9','TS03_S10','TS03_S11']
KlacKmar_RNA = ['TS04_S19','TS04_S20','TS04_S21']
KlacKwic_RNA = ['TS04_S13','TS04_S14','TS04_S15','TS04_S16','TS04_S17','TS04_S18','TS04_S22','TS04_S23','TS04_S24']
KlacKwic_good_RNA = ['TS04_S13','TS04_S14','TS04_S15','TS04_S16','TS04_S17','TS04_S22','TS04_S23']
KlacKwic_best_RNA = ['TS04_S13','TS04_S15','TS04_S16','TS04_S17']

average_number_reads = 22115774

def write_orthology_file():
    #takes the annotations for Kmar, and Kwic, and outputs a file with each gene name in a 
    #different column with the first, second column the Scerevisiae ortholog, third is Klac ortholog
    inFile = open('/Users/trevorsorrells/Documents/UCSF/Projects/genome_info/genomes_ygob/Pillars.tab.txt','r')
    lines_for_output = []
    for line in inFile:
        line_split = line.split()
        if line_split[15][:3] == 'KLL':
            output_line = ['.','.',line_split[15],'.','.']
            if line_split[11][0] == 'Y' or line_split[11][:3] == 'Sce':
                output_line[0] = line_split[11]
            if line_split[21][0] == 'Y' or line_split[21][:3] == 'Sce':
                output_line[1] = line_split[21]
            if output_line[0] == '.' and output_line[1] == '.':
                continue
            lines_for_output.append(output_line)
    inFile.close()
    inFile = open(project_dir + 'genomes_9_21_2015/km02.annotation_final.ygap.txt','r')
    for line in inFile:
        line_split = line.split('\t')
        if len(line_split[8]) > 1:
            gene_names = re.split(' |\|', line_split[8])
            for output_line in lines_for_output:
                for gene_name in gene_names:
                    if output_line[0] == gene_name or output_line[1] == gene_name:
                        output_line[3] = line_split[0]
                        break
    inFile.close()
    inFile = open(project_dir + 'genomes_9_21_2015/kwic.annotation_final.ygap.txt','r')
    for line in inFile:
        line_split = line.split('\t')
        if len(line_split[8]) > 1:
            gene_names = re.split(' |\|', line_split[8])
            for output_line in lines_for_output:
                for gene_name in gene_names:
                    if output_line[0] == gene_name or output_line[1] == gene_name:
                        output_line[4] = line_split[0]
                        break
    inFile.close()
    print len(lines_for_output)
    outFile = open(project_dir + 'genomes_9_21_2015/orthology_list.txt','w')
    for output_line in lines_for_output:
        outFile.write('\t'.join(output_line) + '\n')
    outFile.close()

def normalize_reads_by_experiment():
    #this method takes read counts for each gene from the ASE_server.py script and normalizes
    #based on the total reads in each experiment
    all_files = []
    all_files.extend(KlacKwic_file_names)
    all_files.extend(KlacKmar_file_names)
    all_read_counts = []
    all_read_counts.extend(KlacKwic_read_counts)
    all_read_counts.extend(KlacKmar_read_counts)
    for file_name, read_counts in zip(all_files, all_read_counts):
        normalization_factor = float(average_number_reads)/read_counts
        if file_name in KlacKwic_file_names:
            file_path = project_dir + 'output_files/read_counts/counts_' + file_name + '_KlacKwic_genes_unique.sam.txt'
            out_file_path = project_dir + 'output_files/read_counts_normalized/counts_' + file_name + '_KlacKwic.txt'
        else:
            file_path = project_dir + 'output_files/read_counts/counts_' + file_name + '_KlacKmar_genes_unique.sam.txt'
            out_file_path = project_dir + 'output_files/read_counts_normalized/counts_' + file_name + '_KlacKmar.txt'
        with open(out_file_path, 'w') as outFile:
            with open(file_path, 'r') as inFile:
                for line in inFile:
                    line_split = line.split()
                    norm_read_count = float(line_split[1])*normalization_factor
                    outFile.write(line_split[0] + '\t' + str(norm_read_count) + '\n')
                    
def normalize_reads_by_gDNA_KlacKmar(cutoff=30):
    #This method takes mRNA read counts and normalizes them by gDNA reads per gene
    #cutoff is the minimum number of reads in gDNA below which the gene is discarded from the output due to the concern
    #of dividing by a small number
    KlacKmar_gene_list = []
    KlacKmar_gDNA_bygene = []
    with open(project_dir + 'output_files/read_counts_normalized/counts_' + KlacKmar_gDNA[0] + '_KlacKmar.txt','r') as inFile:
        for line in inFile:
            KlacKmar_gene_list.append(line.split()[0])
    for replicate in KlacKmar_gDNA:
        with open(project_dir + 'output_files/read_counts_normalized/counts_' + replicate + '_KlacKmar.txt','r') as inFile:
            current_read_counts = []
            for line in inFile:
                current_read_counts.append(float(line.split()[1]))
            KlacKmar_gDNA_bygene.append(current_read_counts)
    KlacKmar_mean_gDNA_read_counts = []
    for i in range(len(KlacKmar_gDNA_bygene[0])):
        KlacKmar_mean_gDNA_read_counts.append((KlacKmar_gDNA_bygene[0][i]+KlacKmar_gDNA_bygene[1][i]+KlacKmar_gDNA_bygene[2][i])/3)
    print KlacKmar_mean_gDNA_read_counts[:20]
    for replicate in KlacKmar_RNA:
        print replicate
        total_genes = 0
        with open(project_dir + 'output_files/read_counts_normalized/counts_' + replicate + '_KlacKmar.txt','r') as inFile:
            with open(project_dir + 'output_files/read_counts_over_gDNA/counts_' + replicate + '_KlacKmar.txt','w') as outFile:
                for line, gene_name, gDNA_value in zip(inFile, KlacKmar_gene_list, KlacKmar_mean_gDNA_read_counts):
                    if gDNA_value > cutoff:
                        outFile.write(gene_name + '\t' + str(float(line.split()[1])/gDNA_value) + '\n')
                        total_genes += 1
        print total_genes

def normalize_reads_by_gDNA_KlacKwic(cutoff=30):
    #This method takes mRNA read counts and normalizes them by gDNA reads per gene
    #cutoff is the minimum number of reads in gDNA below which the gene is discarded from the output due to the concern
    #of dividing by a small number
    KlacKwic_gene_list = []
    KlacKwic_gDNA_bygene = []
    with open(project_dir + 'output_files/read_counts_normalized/counts_' + KlacKwic_gDNA[0] + '_KlacKwic.txt','r') as inFile:
        for line in inFile:
            KlacKwic_gene_list.append(line.split()[0])
    for replicate in KlacKwic_gDNA:
        with open(project_dir + 'output_files/read_counts_normalized/counts_' + replicate + '_KlacKwic.txt','r') as inFile:
            current_read_counts = []
            for line in inFile:
                current_read_counts.append(float(line.split()[1]))
            KlacKwic_gDNA_bygene.append(current_read_counts)
    KlacKwic_mean_gDNA_read_counts = []
    for i in range(len(KlacKwic_gDNA_bygene[0])):
        KlacKwic_mean_gDNA_read_counts.append((KlacKwic_gDNA_bygene[0][i]+KlacKwic_gDNA_bygene[1][i]+KlacKwic_gDNA_bygene[2][i])/3)
    print KlacKwic_mean_gDNA_read_counts[:20]
    for replicate in KlacKwic_RNA:
        print replicate
        total_genes = 0
        with open(project_dir + 'output_files/read_counts_normalized/counts_' + replicate + '_KlacKwic.txt','r') as inFile:
            with open(project_dir + 'output_files/read_counts_over_gDNA/counts_' + replicate + '_KlacKwic.txt','w') as outFile:
                for line, gene_name, gDNA_value in zip(inFile, KlacKwic_gene_list, KlacKwic_mean_gDNA_read_counts):
                    if gDNA_value > cutoff:
                        outFile.write(gene_name + '\t' + str(float(line.split()[1])/gDNA_value) + '\n')
                        total_genes += 1
        print total_genes

    

def calculate_ASE_individual_reps(replicate_list, second_species):
    #uses the YGAP input as the ortholog mapping, then outputs Scer, other gene names, gene-pair's average Klac/2nd species ratio for each replicate
    #second species is either 'Kwic' or 'Kmar'
    gene_list = []
    with open(project_dir + 'output_files/read_counts_over_gDNA/counts_' + replicate_list[0] + '_Klac'+ second_species + '.txt','r') as inFile:
        for line in inFile:
            gene_list.append(line.split()[0])
    replicate_values = []
    for replicate in replicate_list:
        with open(project_dir + 'output_files/read_counts_over_gDNA/counts_' + replicate + '_Klac'+ second_species + '.txt','r') as inFile:
            current_read_counts = []
            for line in inFile:
                current_read_counts.append(float(line.split()[1]))
            replicate_values.append(current_read_counts)
    print replicate_values[0][:20]
    ortholog_list = []
    with open(project_dir + '/genomes_9_21_2015/orthology_list.txt','r') as inFile:
        for line in inFile:
            ortholog_list.append(line.split())
    print ortholog_list[:10]
    #go through ortholog list, and for each pair of genes 
    if second_species == 'Kmar':
        sec_species_ortho_index = 3
    elif second_species == 'Kwic':
        sec_species_ortho_index = 4
    total_gene_pairs = 0
    with open(project_dir + 'output_files/ASE_Klac'+ second_species + '.txt','w') as outFile:
        for orthogroup in ortholog_list:
            if orthogroup[2] in gene_list and orthogroup[sec_species_ortho_index] in gene_list:
                total_gene_pairs += 1
                Klac_values = []
                Klac_index = gene_list.index(orthogroup[2])
                for replicate in replicate_values:
                    Klac_values.append(replicate[Klac_index])
                sec_species_values = []
                sec_species_gene_index = gene_list.index(orthogroup[sec_species_ortho_index])
                for replicate in replicate_values:
                    sec_species_values.append(replicate[sec_species_gene_index])
                ASE_ratios = []
                for Klac_value, sec_species_value in zip(Klac_values, sec_species_values):
                    if sec_species_value != 0.0:
                        ASE_ratios.append(Klac_value/sec_species_value)
                    else: ASE_ratios.append('undefined')
                outFile.write('\t'.join(orthogroup[:3]) + '\t' + orthogroup[sec_species_ortho_index] +'\t' + '\t'.join(map(str, ASE_ratios)) + '\n')
    print 'total ortholog pairs examined: ' + str(total_gene_pairs)


def perform_RP_enrichment():
    #calculate enrichment of the RP subunits to test for differential expression
    gene_list_names = ['large_subunit','small_subunit']
    gene_list_path = '/Users/trevorsorrells/Documents/UCSF/Projects/TFBS_scoring/ribosome_project/go_terms/'
    gene_list = []
    for gene_file in gene_list_names:
        handle = open(gene_list_path + gene_file + '.txt', 'r')
        for line in handle:
            current_gene = line.split()[1]
            if current_gene not in gene_list:
                gene_list.append(current_gene)
        handle.close()
    ASE_file_data = []
    with open(project_dir + 'output_files/ASE/ASE_KlacKmar.txt','r') as inFile:
        for line in inFile:
            line_split = line.split()
            line_split[4] = float(line_split[4])
            line_split[5] = float(line_split[5])
            ASE_file_data.append(line_split)
    calc_enrichment_geneset(gene_list, 'RPs', ASE_file_data, 1.1, .9090909)


def calc_enrichment_geneset(Scer_gene_list, list_name, ASE_data, min_ratio, max_ratio):
    #NEW METHOD USES SIMPLIFIED ASE_data AND LOG VALUES INSTEAD OF DIRECT RATIOS
    #uses hypergeometric test to determine whether the list of orthologs in Scer_gene_list
    #is shows allele-specific expression in Klac relative to Kwic or Kmar
    #min_ratio and max_ratio are the minimum difference in expression to be considered
    #numbers in 2x2 table used in fisher's exact test are: list|not list
    #                                                up:     a | b
    #                                                not up: c | d
    #first calculate mean of ratios
    ASE_values = []
    for ortho_pair in ASE_data:
        if ortho_pair[0] in Scer_gene_list or ortho_pair[1] in Scer_gene_list:
            ASE_values.append(ortho_pair[2])
    if len(ASE_values) == 0:
        return 'no orthologs'
    mean = sum(ASE_values)/len(ASE_values)
    #first calculate  P-value for up in Klac
    a,b,c,d = 0,0,0,0
    for ortho_pair in ASE_data:
        if ortho_pair[0] in Scer_gene_list or ortho_pair[1] in Scer_gene_list:
            if ortho_pair[3] == 'significant' and ortho_pair[2] > min_ratio:
                a += 1
            else:
                c += 1
        else:
            if ortho_pair[3] == 'significant' and ortho_pair[2] > min_ratio:
                b += 1
            else:
                d += 1
    oddsratio1, pvalue1 = stats.fisher_exact([[a, b], [c, d]])
    #next calculate P-value for down in Klac
    e,f,g,h = 0,0,0,0
    for ortho_pair in ASE_data:
        if ortho_pair[0] in Scer_gene_list or ortho_pair[1] in Scer_gene_list:
            if ortho_pair[3] == 'significant' and ortho_pair[2] < max_ratio:
                e += 1
            else:
                g += 1
        else:
            if ortho_pair[3] == 'significant' and ortho_pair[2] < max_ratio:
                f += 1
            else:
                h += 1
    oddsrati2o, pvalue2 = stats.fisher_exact([[e, f], [g, h]])
    #only report up if mean is higher, and down if mean is lower
    if mean > 0:
        return [pvalue1, mean, min_ratio, a, b, c, d, list_name]
    else:
        return [pvalue2, mean, max_ratio, e, f, g, h, list_name]


def perform_GO_enrichment(ASE_file_name, min_ratio1, max_ratio1):
    #NOTE THIS METHOD USES LOG2(ASE) VALUES NOT RAW RATIOS!!!
    #this takes as input the output of FDR calculation using individual reps
    #takes all Scer GO categories and tests for up or down regulation in Klac relative to 2nd species
    #then ranks all of these categories and outputs the most up or down regulated. 
    #first make a list of all the GO terms to analyze and make lists of genes that correspond to them
    GO_IDs = []
    GO_gene_lists = []
    with open('/Users/trevorsorrells/Documents/UCSF/Projects/genome_info/GO_analysis/gene_association.sgd','r') as inFile:
        for line in inFile:
            if line[0] != '!':
                line_split = line.split('\t')
                GO_ID = line_split[4]
                gene_name = line_split[10].split('|')[0]
                if GO_ID in GO_IDs:
                    GO_gene_lists[GO_IDs.index(GO_ID)].append(gene_name)
                else:
                    GO_IDs.append(GO_ID)
                    GO_gene_lists.append([gene_name])
    print str(len(GO_IDs)) + ' total GO terms'
    #load ASE file of interest
    ASE_file_data = []
    with open(project_dir + 'output_files/ASE/'+ASE_file_name,'r') as inFile:
        for line in inFile:
            line_split = line.split()
            #calculate average ASE value
            ASE_values = []
            for ASE_value in line_split[3:len(line_split)-1]:
                ASE_values.append(float(ASE_value))
            ASE_file_data.append([line_split[0],line_split[1],sum(ASE_values)/len(ASE_values), line_split[-1]])
    #grab GO descriptions from other file
    GO_dict = {}
    with open('/Users/trevorsorrells/Documents/UCSF/Projects/genome_info/GO_analysis/go_terms.tab.txt','r') as inFile:
        for line in inFile:
            line_split = line.split('\t')
            while len(line_split[0]) < 7:
                line_split[0] = '0' + line_split[0]
            line_split[0] = 'GO:' + line_split[0]
            GO_dict[line_split[0]] = line_split[1]
    #perform GO tests
    all_GO_tests = []
    for GO_ID, gene_list in zip(GO_IDs, GO_gene_lists):
        enrichment = calc_enrichment_geneset(gene_list, GO_ID, ASE_file_data, min_ratio1, max_ratio1)
        if enrichment == 'no orthologs':
            continue
        all_GO_tests.append(enrichment)
        all_GO_tests[-1].append(GO_dict[GO_ID])
    sorted_GO_tests = sorted(all_GO_tests)
    m = len(sorted_GO_tests)
    print str(m) + ' total tests performed'
    k = 1
    with open(project_dir + 'output_files/GO_tests/GO_'+ASE_file_name,'w') as outFile:
        outFile.write('\t'.join(['pvalue', 'mean', 'min_or_max_ratio', 'a', 'b', 'c', 'd', 'GO_ID','description'])+'\n')
        for test in sorted_GO_tests:
            if float(k)/m*0.05 > test[0]:
                test.append('significant')
            else:
                test.append('not')
            outFile.write('\t'.join(map(str, test)) + '\n')
            k += 1
    #also test RPs
    gene_list_names = ['large_subunit','small_subunit']
    gene_list_path = '/Users/trevorsorrells/Documents/UCSF/Projects/TFBS_scoring/ribosome_project/go_terms/'
    gene_list = []
    for gene_file in gene_list_names:
        handle = open(gene_list_path + gene_file + '.txt', 'r')
        for line in handle:
            current_gene = line.split()[1]
            if current_gene not in gene_list:
                gene_list.append(current_gene)
        handle.close()
    print '\t'.join(map(str,calc_enrichment_geneset(gene_list, GO_ID, ASE_file_data, min_ratio1, max_ratio1)))

def define_genes_up_Kwic(rep_indices, FDR_rate):
    #loads the default file for Kmar, Kwic with 3 samples (aka different hybrids) or replicates
    gene_ASE_values = [] #will be a list of genes, with each genes having fpkm values for 3 replicates
    with open(project_dir + '/output_files/ASE/ASE_KlacKwic_individual_reps.txt', 'r') as inFile:
        for line in inFile:
            line_split = line.split('\t')
            current_gene = []
            for value in line_split[4:]:
                try:
                    current_gene.append(math.log(float(value),2))
                except ValueError:
                    current_gene.append(0.0)
            current_gene.append(line_split[0])
            current_gene.append(line_split[1])
            gene_ASE_values.append(current_gene)
    i = 0
    while i < len(gene_ASE_values):
        values_to_use = []
        for rep_index in rep_indices:
            values_to_use.append(gene_ASE_values[i][rep_index])
        output = stats.ttest_1samp(values_to_use, 0.0)
        if math.isnan(output[1]):
            gene_ASE_values.remove(gene_ASE_values[i])
        else:
            gene_ASE_values[i] = [output[1], gene_ASE_values[i][-2], gene_ASE_values[i][-1]]
            gene_ASE_values[i].extend(values_to_use)
            i += 1
    sorted_values = sorted(gene_ASE_values)
    #now perform Benjamini Hochberg correction
    m = len(sorted_values)
    k = 1
    significant = True
    print str(m) + ' total orthologs examined'
    significant_count = 0
    for gene_ASE_value in sorted_values:
        if gene_ASE_value[0] < float(k)/m*FDR_rate and significant:
             gene_ASE_value.append('significant')
             significant_count += 1
        else:
            significant = False
            gene_ASE_value.append('not')
        k += 1
    print str(significant_count) + ' genes showed ASE'
    print sorted_values[:10]
    with open(project_dir + '/output_files/ASE/ASE_KlacKwic'+ str(len(rep_indices))+'_BH_FDR.txt', 'w') as outFile:
        for gene_ASE_value in sorted_values:
            outFile.write('\t'.join(map(str,gene_ASE_value)) + '\n')


def define_genes_up(sec_species_name, FDR_rate):
    #loads the default file for Kmar, Kwic with 3 samples (aka different hybrids) or replicates
    gene_ASE_values = [] #will be a list of genes, with each genes having fpkm values for 3 replicates
    with open(project_dir + '/output_files/ASE/ASE_Klac'+sec_species_name+'_individual_reps.txt', 'r') as inFile:
        for line in inFile:
            line_split = line.split('\t')
            current_gene = []
            for value in line_split[4:]:
                try:
                    current_gene.append(math.log(float(value),2))
                except ValueError:
                    current_gene.append(0.0)
            current_gene.append(line_split[0])
            current_gene.append(line_split[1])
            gene_ASE_values.append(current_gene)
    i = 0
    while i < len(gene_ASE_values):
        output = stats.ttest_1samp(gene_ASE_values[i][:3], 0.0)
        if math.isnan(output[1]):
            gene_ASE_values.remove(gene_ASE_values[i])
        else:
            gene_ASE_values[i] = [output[1], gene_ASE_values[i][3],gene_ASE_values[i][4], gene_ASE_values[i][0], gene_ASE_values[i][1], gene_ASE_values[i][2]]
            i += 1
    sorted_values = sorted(gene_ASE_values)
    #now perform Benjamini Hochberg correction
    m = len(sorted_values)
    k = 1
    significant = True
    print str(m) + ' total orthologs examined'
    significant_count = 0
    for gene_ASE_value in sorted_values:
        if gene_ASE_value[0] < float(k)/m*FDR_rate and significant:
             gene_ASE_value.append('significant')
             significant_count += 1
        else:
            significant = False
            gene_ASE_value.append('not')
        k += 1
    print str(significant_count) + ' genes showed ASE'
    print sorted_values[:10]
    with open(project_dir + '/output_files/ASE/ASE_Klac'+sec_species_name+'_BH_FDR.txt', 'w') as outFile:
        for gene_ASE_value in sorted_values:
            outFile.write('\t'.join(map(str,gene_ASE_value)) + '\n')
    with open(project_dir + '/output_files/ASE/ASE_Klac'+sec_species_name+'_BH_FDR_average.txt', 'w') as outFile:
        for gene_ASE_value in sorted_values:
            outFile.write('\t'.join(map(str,gene_ASE_value)) + '\n')


write_orthology_file()
normalize_reads_by_experiment()
normalize_reads_by_gDNA_KlacKmar()
normalize_reads_by_gDNA_KlacKwic()

calculate_ASE_individual_reps(KlacKwic_RNA, 'Kwic')
calculate_ASE_individual_reps(KlacKmar_RNA, 'Kwic')

define_genes_up('Kmar', 0.05)
define_genes_up_Kwic([0,2,3,4], 0.05)

perform_GO_enrichment('ASE_KlacKmar_BH_FDR.txt', 0.13750352375, -0.13750352375)
perform_GO_enrichment('ASE_KlacKwic7_BH_FDR.txt', 0.13750352375, -0.13750352375)
