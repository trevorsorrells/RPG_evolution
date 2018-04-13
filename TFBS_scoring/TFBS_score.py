#!/usr/bin/env python
################################################################################
# TFBS_score.py
#
# A set of tools for scoring TFBSs in a set of upstream regions 
#
################################################################################

import time
import math
from scipy import stats
import gzip

project_folder = '/Users/trevorsorrells/Documents/UCSF/Projects/TFBS_scoring/ribosome_project/'
genome_info = '/Users/trevorsorrells/Documents/UCSF/Projects/genome_info'

#list of ribosomal protein genes
inFile = open('ribosome_project/lists/RPs.txt','r')
gene_list = inFile.read().split()
gene_list = sorted(gene_list)
inFile.close()

#dicitonary for accounting for changing genus names over time
old_new_dict = {'Kpol':'Vpol','Kthe':'Lthe','Kwal':'Lwal','Sklu':'Lklu','Agos':'Egos','Scas':'Ncas'}

#for reverse complement
base_complement = {'a':'t','t':'a','g':'c','c':'g','A':'T','T':'A', 'G':'C','C':'G','N':'N','W':'W','M':'K','K':'M','R':'Y','Y':'R','B':'V','V':'B','S':'S','D':'H','H':'D'}

#relevant subsets of species
post_WGD_species = ['Scer','Suva','Cgla','Scas','Vpol','Ndai','Ncas','Knag','Kafr','Tpha','Tbla']
kluyveromyces_zygo_species = ['Zrou','Klac','Kwic','Egos','Lthe','Lklu','Lwal','Tdel','Ecym']
CTG_species = ['Calb','Cdub','Ctro','Cpar','Lelo','Clus','Cgui','Dhan','Psti','Mbic','Cten','Hbur','Spas','Ctan']
albicans_clade = ['Calb','Cdub','Ctro','Cpar','Cort','Lelo','Spas']

#list of species for analysis
with open('ribosome_project/lists/all_species_in_order.txt','r') as inFile:
    species_list_long = inFile.read().split('\n')
species_output_order = []
for long_species in species_list_long:
    name_split = long_species.split()
    species_output_order.append(name_split[0][0] + name_split[1][:3])
print len(species_output_order)
print species_output_order

def reverseComplement(sequence):
    #returns the reverse complement of input sequence
    letters = list(sequence)
    letters.reverse()
    letters = [base_complement[base] for base in letters]
    return ''.join(letters)

def generateMotifDict(motifSourc):
    #returns a tuple of dictionaries where each dictionary corresponds to a position in the motif.
    #Also returns the reverse complement motif for scoring the same DNA sequence
    # + and - strands at the same time
    motif_file = open(motifSourc, 'r')
    weight_list = []
    for line in motif_file:
        weight_list.append(line.split())
    motif_positions = []
    motif_positions_rc = []
    for i in range(len(weight_list)):
        motif_positions.append({'A':float(weight_list[i][0]),'C':float(weight_list[i][1]),'G':float(weight_list[i][2]),'T':float(weight_list[i][3]),'R':-5.0,'N':-5.0,'W':-5.0,'M':-5.0,'K':-5.0,'B':-5.0,'V':-5.0,'Y':-5.0,'S':-5.0,'D':-5.0,'H':-5.0})
        motif_positions_rc.append({'T':float(weight_list[-1-i][0]),'G':float(weight_list[-1-i][1]),'C':float(weight_list[-1-i][2]),'A':float(weight_list[-1-i][3]),'R':-5.0,'N':-5.0,'W':-5.0,'M':-5.0,'K':-5.0,'B':-5.0,'V':-5.0,'Y':-5.0,'S':-5.0,'D':-5.0,'H':-5.0})
    motif_file.close()
    max_score = 0.0
    for position in motif_positions:
        max_score += max(position.values())
    print 'max score ' + str(max_score)
    return tuple(motif_positions), tuple(motif_positions_rc)

def scoreNew(motif, motif_rc, seq_representation, cutoff=0):
    #returns the top score or the number of scores above a given cutoff for a certain intergenic region
    #or if cutoff=0, returns whole list of sequences (default)
    #is written only to be used with log-odds matrices!!
    #location information is midpoint of each motif 
    orf_scores=[]
    size = len(motif)
    for intergenic in seq_representation:
        score_list = []
        [score_list.append((-100.0,0,'','')) for x in range(100)]
        for i in range(len(intergenic)-size):
            motif_score = 0.0
            motif_rc_score = 0.0
            for j in range(size):
                motif_score += motif[j][intergenic[i+j]]
                motif_rc_score += motif_rc[j][intergenic[i+j]]
            if motif_score > score_list[0][0]:
                motif_hit_seq = intergenic[i:i+size]
                score_list.append((motif_score,i+size/2,''.join(motif_hit_seq),'+'))
                score_list = sorted(score_list)[1:]
            if motif_rc_score > score_list[0][0]:
                motif_hit_seq = intergenic[i:i+size]
                score_list.append((motif_rc_score,i+size/2,reverseComplement(''.join(motif_hit_seq)),'-'))
                score_list = sorted(score_list)[1:]
        orf_scores.append(score_list)
    if cutoff==0:
        return orf_scores
    elif cutoff == 'max':
        return orf_scores[-1][0]
    else:
        n_hits = 0
        n=-1
        while orf_scores[n][0] > cutoff:
            n_hits += 1
            n += -1
        return n_hits

def scoreForward(motif, motif_rc, seq_representation, cutoff=0):
    #returns the top score or the number of scores above a given cutoff for a certain intergenic region
    #or if cutoff=0, returns whole list of sequences (default)
    #only returns directional binding sites in the forward direction (requires motif & intergenic region to be directional)
    #the part that is new is that the scoring is scaled by bits of information given background frequencies
    #is written only to be used with log-odds matrices!!
    #location information is midpoint of each motif 
    orf_scores=[]
    size = len(motif)
    for intergenic in seq_representation:
        score_list = []
        [score_list.append((0.0,0,'')) for x in range(60)]
        for i in range(len(intergenic)-size):
            motif_score = 0.0
            for j in range(size):
                motif_score += motif[j][intergenic[i+j]]
            if motif_score > score_list[0][0]:
                motif_hit_seq = intergenic[i:i+size]
                score_list.append((motif_score,i+size/2,''.join(motif_hit_seq),'+'))
                score_list = sorted(score_list)[1:]
        orf_scores.append(score_list)
    if cutoff==0:
        return orf_scores
    elif cutoff == 'max':
        return orf_scores[-1][0]
    else:
        n_hits = 0
        n=-1
        while orf_scores[n][0] > cutoff:
            n_hits += 1
            n += -1
        return n_hits

def scoreAndOutputAllGenes(intergenic_folder, maxLen, motifSource, output_path, cutoff1=0, startPosition=0):
    #scores all promoters in the genome and outputs information for all individual genes
    if cutoff1 == 0:
        out_file = open(output_path + '_max_score.txt','w')
    else:
        out_file = open(output_path + str(cutoff1)+ 'cutoff.txt','w')
    motif_dict, motif_dict_rc = generateMotifDict(motifSource)
    for species in species_output_order:
        print species
        intergenics_file = open(intergenic_folder +species + '_intergenics.faa','r')
        seq = ''
        last_gene_info = [intergenics_file.readline()[1:].strip()]
        for line in intergenics_file:
            if line[0] == '>':
                if len(seq) > maxLen:
                    seq = seq[-maxLen:]
                last_gene_info.append(tuple(seq.upper()))
                out_file.write(species + '\t' + last_gene_info[0] + '\t')
                #score last_gene_info
                if len(seq) > 0:
                    cur_gen_scores = scoreNew(motif_dict, motif_dict_rc, [last_gene_info[1]], 0)
                    if cutoff1 == 0:
                        #either write all information, or just top location if file to be used for relative locations
                        #if cur_gen_scores[0][-1][0] > 0.0:
                            #out_file.write(str(cur_gen_scores[0][-1][1]))
                        out_file.write(str(cur_gen_scores[0][-1][0]) +'\t'+ str(cur_gen_scores[0][-1][1]) +'\t'+ cur_gen_scores[0][-1][2] +'\t'+ cur_gen_scores[0][-1][3])
                    else:
                        n = 1
                        scores = []
                        while cur_gen_scores[0][-n][0] > cutoff1:
                            if cur_gen_scores[0][-n][1] not in scores:
                                scores.append(cur_gen_scores[0][-n][1])
                            n += 1
                            if n == 100:
                                print "more than 100 sites:"
                                print last_gene_info[0]
                                break
                        out_file.write('\t'.join(map(str, scores)))
                out_file.write('\n')
                seq = ''
                last_gene_info = [line[1:].strip()]
            else:
                seq += line.strip()
        #account for last gene in file
        if len(seq) > maxLen:
            seq = seq[-maxLen:]
        last_gene_info.append(tuple(seq.upper()))
        out_file.write(species + '\t' + last_gene_info[0] + '\t')
        #score last_gene_info
        if len(seq) > 0:
            cur_gen_scores = scoreNew(motif_dict, motif_dict_rc, [last_gene_info[1]], 0.0)
            if cutoff1 == 0:
                if cur_gen_scores[0][-1][0] > 0.0:
                    #either write all information, or just top location if file to be used for relative locations
                    #out_file.write(str(cur_gen_scores[0][-1][0]) +'\t'+ str(cur_gen_scores[0][-1][1]) +'\t'+ cur_gen_scores[0][-1][2] +'\t'+ cur_gen_scores[0][-1][3])
                    out_file.write(str(cur_gen_scores[0][-1][1]))
            else:
                n = 1
                while cur_gen_scores[0][-n][0] > cutoff1:
                    if cur_gen_scores[0][-n][1] not in scores:
                        scores.append(cur_gen_scores[0][-n][1])
                    n += 1
                out_file.write('\t'.join(map(str, scores)))
        out_file.write('\n')
    out_file.close()

#START HERE WITH CONTINUED EDITING


def genomeWideMotifCounts(genome_wide_score_output):
    #takes a set of genes and whole-genome scores for a TF, then returns two lists containing the 
    #number of motifs found for each gene in the gene set and genome-wide

    #retrieve gene list (ribosomal protein) names
    gene_list_by_species = []
    [gene_list_by_species.append([]) for species in species_output_order]
    
    for species_name, current_gene_list in zip(species_output_order, gene_list_by_species):
        with open(project_folder + '/intergenics/by_species/' + species_name, 'r') as inFile:
            for line in inFile:
                if line[0] == '>':
                    if line.split()[2] in gene_list:
                        current_gene_list.append(line.split()[1].split('|')[-1])

    #go through genome-wide scores counting the number of motifs in each gene
    #place each count into either the gene set, or all other gene list of the given species
    gene_set_table = []
    other_genes_table = []
    [gene_set_table.append([]) for species in species_output_order]
    [other_genes_table.append([]) for species in species_output_order]

    species_scored_but_not_output = []
    with open(project_folder + '/outputs/whole_genome_scores_1000/' + genome_wide_score_output + '.txt', 'r') as inFile:
        for line in inFile:
            line_split = line.split()
            if len(line_split) < 2:
                print line_split
            if line_split[0] not in species_output_order: 
                if line_split[0] not in species_scored_but_not_output:
                    species_scored_but_not_output.append(line_split[0])
                continue
            current_species_index = species_output_order.index(line_split[0])
            if line_split[1][:3] == 'BN1':
                line_split[1] = line_split[1][6:]
            if line_split[1] in gene_list_by_species[current_species_index]:
                gene_set_table[current_species_index].append(len(line_split) - 2)
            else:
                other_genes_table[current_species_index].append(len(line_split) - 2)
    print 'species scored but not in output list:'
    print species_scored_but_not_output
    return gene_set_table, other_genes_table

def genomeWideEnrichment(gene_set_counts, other_gene_counts, proportion_cutoff=1, enrichment_cutoff=None):
    #takes output from genomeWideMotifCounts() and calculates the proportion of genes in the gene set with
    #at least proportion_cutoff motifs, and then calculates the maximum enrichment of motifs relative to the rest of the genome
    #and the number of motifs at which the maximum enrichment is found
    #if cutoff is given, then 
    #output is a list, so multiple different motifs and scoring cutoffs can be output into a single file
    genes_in_set = []
    proportions_gene_set = []
    proportions_other_genes = []
    enrichment_p_values = []
    cutoffs = []
    for species_name, gene_set_species, other_gene_species in zip(species_output_order, gene_set_counts, other_gene_counts):
        
        gene_set_species.sort()
        other_gene_species.sort()
        genes_in_set.append(len(gene_set_species))

        #first calculate proportions of gene set and other genes that have a motif over proportion cutoff
        if proportion_cutoff not in gene_set_species:
            proportions_gene_set.append(0.0)
        else:
            proportion_index = gene_set_species.index(proportion_cutoff)
            proportions_gene_set.append(float(len(gene_set_species)-proportion_index)/len(gene_set_species))
        if proportion_cutoff not in other_gene_species:
             proportions_other_genes.append(0.0)
        else:
            proportion_index = other_gene_species.index(proportion_cutoff)
            proportions_other_genes.append(float(len(other_gene_species)-proportion_index)/len(other_gene_species))
        
        #next calculate enrichment if cutoff is fixed
        if enrichment_cutoff != None:
            if enrichment_cutoff not in gene_set_species:
                a = 0
                c = len(gene_set_species)
            else:
                a = gene_set_species.index(enrichment_cutoff)
                c = len(gene_set_species)-gene_set_species.index(enrichment_cutoff)
            if enrichment_cutoff not in other_gene_species:
                b = 0
                c = len(other_gene_species)
            else:
                b = other_gene_species.index(enrichment_cutoff)
                d = len(other_gene_species)-other_gene_species.index(enrichment_cutoff)
            cutoffs.append(enrichment_cutoff)
            enrichment_p_values.append(-math.log10(stats.fisher_exact([[a,b],[c,d]])[1]))
        else:
            #calculate enrichment if cutoff changes to maximize log10(P)
            max_count = gene_set_species[-1]
            current_cutoff = 1
            max_log10_p = 0.0
            chosen_cutoff = 0
            while current_cutoff <= max_count:
                if current_cutoff not in gene_set_species:
                    a = 0
                    c = len(gene_set_species)
                else:
                    a = gene_set_species.index(current_cutoff)
                    c = len(gene_set_species)-gene_set_species.index(current_cutoff)
                if current_cutoff not in other_gene_species:
                    b = 0
                    c = len(other_gene_species)
                else:
                    b = other_gene_species.index(current_cutoff)
                    d = len(other_gene_species)-other_gene_species.index(current_cutoff)
                current_p_value = -math.log10(stats.fisher_exact([[a,b],[c,d]])[1])
                if current_p_value > max_log10_p:
                    max_log10_p = current_p_value
                    chosen_cutoff = current_cutoff
                current_cutoff += 1
            cutoffs.append(chosen_cutoff)
            enrichment_p_values.append(max_log10_p)
                
    return genes_in_set, proportions_gene_set, proportions_other_genes, enrichment_p_values, cutoffs


def scoreAndOutputAllEnrichment(motif_output_file_names,proportion_cutoff_list, enrichment_cutoff_list, output_name):
    #outputs the enrichment for motifs relative to the rest of the genome for each species as a single file
    with open(project_folder + '/outputs/whole_genome_scores_1000/' + output_name + '.txt', 'w') as outFile:
        outFile.write('species\t' + '\t'.join(species_output_order) + '\n\n')
        for output, prop_cut, enrich_cut in zip(motif_output_file_names, proportion_cutoff_list, enrichment_cutoff_list):
            outFile.write(output + '\n')
            gene_counts,other_counts = genomeWideMotifCounts(output)
            n_genes, proportions_set, proportions_other, p_values, cutoff_list = genomeWideEnrichment(gene_counts, other_counts, prop_cut, enrich_cut)
            outFile.write('n_genes\t' + '\t'.join(map(str, n_genes)) + '\n')
            outFile.write('proportions_RPs\t' + '\t'.join(map(str, proportions_set)) + '\n')
            outFile.write('proportions_others\t' + '\t'.join(map(str, proportions_other)) + '\n')
            outFile.write('hypergeo\t' + '\t'.join(map(str, p_values)) + '\n')
            outFile.write('cutoff\t' + '\t'.join(map(str, cutoff_list)) + '\n')


#ribosomal protein project analyses

#score 1000bp upstream of  all genes for all motifs at a cutoff
motif_list = ['CGACAAC_Bcin','cbf1_Spas','dot6_scertf','fhl1_scertf','hmo1_Lthe','mcm1_scertf','rap1_scertf','rim101_scertf','rrn7_Aade','sfp1_scertf','stb3_scertf','tbf1_Arub','tbf1_Calb','tbf1_Ptan','tbf1_Sjap']
for motif in motif_list:
    current_cutoff = 6.0 #or 8.0 for more information-rich motifs
    start = time.time()
    scoreAndOutputAllGenes('../genome_info/intergenics/', 1000, 'ribosome_project/motifs/'+ motif +'.txt', 'ribosome_project/outputs/whole_genome_scores_1000/'+ motif + '_', current_cutoff)
    print 'Scoring took ' + str(time.time() - start) + ' seconds'

#ouput the enrichment for all motifs in ribosomal protein genes vs the rest of the genome
motif_output_list = ['CGACAAC_Bcin_6.0cutoff','CGACAAC_Bcin_8.0cutoff','cbf1_Spas_6.0cutoff','cbf1_Spas_8.0cutoff','dot6_scertf_6.0cutoff','dot6_scertf_8.0cutoff','fhl1_scertf_6.0cutoff','fhl1_scertf_8.0cutoff','hmo1_Lthe_6.0cutoff','hmo1_Lthe_8.0cutoff','mcm1_scertf_6.0cutoff','mcm1_scertf_8.0cutoff','rap1_scertf_6.0cutoff','rap1_scertf_8.0cutoff','rim101_scertf_6.0cutoff','rim101_scertf_8.0cutoff','rrn7_Aade_6.0cutoff','rrn7_Aade_8.0cutoff','sfp1_scertf_6.0cutoff','sfp1_scertf_8.0cutoff','stb3_scertf_6.0cutoff','stb3_scertf_8.0cutoff','tbf1_Arub_6.0cutoff','tbf1_Arub_8.0cutoff','tbf1_Calb_6.0cutoff','tbf1_Calb_8.0cutoff','tbf1_Ptan_6.0cutoff','tbf1_Ptan_8.0cutoff','tbf1_Sjap_6.0cutoff','tbf1_Sjap_8.0cutoff']
prop_cutoffs = [1]*30
enrich_cutoffs = [1]*30
scoreAndOutputAllEnrichment(motif_output_list, prop_cutoffs, enrich_cutoffs, 'all_scoring_summary')

#comment out reverse complement scoring in scoreAndOutputAllGenes() to score only Rap1 sites pointing toward transcription start site
start = time.time()
scoreAndOutputAllGenes('/Users/trevorsorrells/Documents/UCSF/Projects/genome_info/intergenics/', 1000, project_folder + '/motifs/rap1_scertf.txt', project_folder + 'outputs/whole_genome_scores_1000/Rap1ForwardMax_', 0.0)
print 'Scoring took ' + str(time.time() - start) + ' seconds'
