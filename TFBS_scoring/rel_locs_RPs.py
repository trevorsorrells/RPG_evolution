#!/usr/bin/env python
################################################################################
# rel_locs_RPs.py
#
# Takes genome wide locations and converts into relative locations graphs
#
################################################################################

import pylab as P
import matplotlib.pyplot as plt
import numpy as np

#directory for project
project_dir = '/Users/trevorsorrells/Documents/UCSF/Projects/TFBS_scoring/ribosome_project/'

#load list of species whose intergenics have been extacted
with open(project_dir + '/lists/all_species_in_order.txt','r') as inFile:
    species_list_long = inFile.read().split('\n')
species_output_order = []
#create abbreviated species list
for long_species in species_list_long:
    name_split = long_species.split()
    species_output_order.append(name_split[0][0] + name_split[1][:3])

#print species_output_order

#load the list of ribosomal protein gene names
inFile = open(project_dir + '/lists/RPs.txt','r')
gene_list = inFile.read().split()
gene_list = sorted(gene_list)
inFile.close()
#add systematic gene IDs to gene_list, because whole genome information uses IDs
for species in species_output_order:
    inFile = open(project_dir + 'intergenics/by_species/' + species, 'r')
    for line in inFile:
        line_split = line.split()
        if len(line_split) > 2:
            if line_split[2] in gene_list:
                if line_split[1][:3] == 'jgi':
                    gene_list.append(species + ' ' + line_split[1].split('|')[-1])
                else:
                    gene_list.append(species + ' ' + line_split[1])


def simple_relative_locs(input_file_path, species_divisions, gene_list='none'): 
    #reads file containing a single set of genes and divides up by species. Returns a vector
    #containing data for each species division
    inFile = open(input_file_path, 'r')
    data_table = []
    weight_table = []
    for sp_dv in species_divisions:
        data_table.append([])
        weight_table.append(0)
    for line in inFile:
        line_split = line.split()
        if gene_list != 'none':
            if line_split[1] not in gene_list:
                continue
        i=2
        while i < len(line_split):
            line_split[i] = int(line_split[i])
            i += 1
        i=0
        for sp_dv in species_divisions:
            if line_split[0] in sp_dv:
                #remove duplicate items by list(set())
                data_table[i].extend(list(set(line_split[2:])))
                #remove 0 if necessary
                #if 0 in data_table[i]:
                #    data_table[i].remove(0)
                weight_table[i] += 1
            i += 1
    inFile.close()
    #repeat weight for each value in the data table
    return data_table, weight_table


def loadYGOBorthos(gene_list_systematic):
    #loads the orthologs of  genes using YGOB for the purposes of studying
    #evolution of Rap1-Mcm1 binding sites in Kazachstania naganishii  at genes other than
    #the ribosomal protein genes
    YGOB_path = '/Users/trevorsorrells/Documents/UCSF/Projects/genome_info/genomes_ygob/'
    #YGOB_gene_species_dict = {'Kpol':'Vpol','TPHA':'Tpha','Tpha':'Tpha','TBLA':'Tbla','Tbla':'Tbla','NDAI':'Ndai','Ndai':'Ndai','NCAS':'Ncas','Ncas':'Ncas','KNAG':'Knag','Knag':'Knag','KAFR':'Kafr','Kafr':'Kafr','CAGL':'Cgla','Cgla':'Cgla','Suva':'Suva','Skud':'Skud','Smik':'Smik','Scer':'Scer','ZYRO':'Zrou','Zrou':'Zrou','TDEL':'Tdel','Tdel':'Tdel','KLLA':'Klac','Klac':'Klac','Ecym':'Ecym','SAKL':'Lklu','Lklu':'Lklu','Sklu':'Lklu','KLTH':'Lthe','Kthe':'Lthe','Kwal':'Lwal'}
    #not_in_species_dict = ['Scer','Egos','Anc']
    
    ortho_list = []
    [ortho_list.append([]) for x in gene_list_systematic]
    genes_in_ortho_list = []
    ortho_file = open(YGOB_path + 'Pillars.tab.txt', 'r')
    for line in ortho_file:
        split_line = line.split()
        #while '---' in split_line:
        #    split_line.remove('---')
        for gene in split_line:
            if gene in gene_list_systematic:
                i = gene_list_systematic.index(gene)
                genes_in_ortho_list.append(gene)
                ortho_list[i] = split_line
    ortho_file.close()
    print 'number of genes in gene list:'
    print len(gene_list)
    print 'number of orthogroups found in YGOB:'
    print len(genes_in_ortho_list)
    
    for orthogroup in ortho_list:
        if len(orthogroup) < 2:
            print orthogroup
            print 'was not in YGOB'
            orthogroup.append('not in YGOB')
    return ortho_list

def output_values_YGOB_species(gene_list, genes_by_orthos):
    #outputs the position weight matrix score of Rap1 and Mcm1 binding sites at non-ribosomal protein genes
    #in post-whole genome hybridization species 

    YGOB_species = ['Cgla','Ecym','Egos','Kafr','Klac','Knag','Lklu','Lthe','Lwal','Ncas','Ndai','Scer','Skud','Smik','Suva','Tbla','Tdel','Tpha','Vpol','Zrou']

    output_lines1 = []
    all_genes_list = []
    for ortho_group in genes_by_orthos:
        current_structure = []
        for gene in ortho_group:
            current_structure.append('')
            all_genes_list.append(gene)
        output_lines1.append(current_structure)
    inFile1 = open(project_dir + 'outputs/whole_genome_scores_500/rap1_scertf__max_score.txt', 'r')
    locs1 = []
    [locs1.append(line.split()) for line in inFile1]
    inFile.close()
    for score in locs1:
        if score[1] in all_genes_list:
            for i, ortho_group in enumerate(genes_by_orthos):
                for j, gene in enumerate(ortho_group):
                    if gene == score[1]:
                        if len(score)>2:
                            output_lines1[i][j] = score[2]

    output_lines2 = []
    for ortho_group in genes_by_orthos:
        current_structure = []
        for gene in ortho_group:
            current_structure.append('')
            all_genes_list.append(gene)
        output_lines2.append(current_structure)
    inFile2 = open(project_dir + 'outputs/whole_genome_scores_500/mcm1_scertf__max_score.txt', 'r')
    locs2 = []
    [locs2.append(line.split()) for line in inFile2]
    inFile.close()
    for score in locs2:
        if score[1] in all_genes_list:
            for i, ortho_group in enumerate(genes_by_orthos):
                for j, gene in enumerate(ortho_group):
                    if gene == score[1]:
                        if len(score)>2:
                            output_lines2[i][j] = score[2]
    outFile = open(project_dir + 'outputs/Knag_site_evolution/Rap1_Mcm1_scores_500.txt','w')
    outFile.write('Rap1\n')
    for gene, scores in zip(gene_list, output_lines1):
        outFile.write(gene + '\t' + '\t'.join(scores) + '\n')
    outFile.write('Mcm1\n')
    for gene, scores in zip(gene_list, output_lines2):
        outFile.write(gene + '\t' + '\t'.join(scores) + '\n')
    outFile.close()

def generateGraph(data_vect, weight_vect, min_range, max_range, output_path, output_name, heatmap=False, font_size=10):
    #outputs heatmaps or histograms showing the relative spacing between regulator binding sites
    #in the regulatory regions of sets of genes
    
    #remove site locations that are out of range of graph
    for column in data_vect:
        i=0
        while i < len(column):
            if column[i] >max_range or column[i] < min_range:
                column.remove(column[i])
            else:
                i += 1
    
    #add pseudocounts to avoid errors in plotting heatmaps
    for i in range(len(data_vect)):
        if len(data_vect[i]) == 0:
            data_vect[i].append(0)
            weight_vect[i] = max(weight_vect)
    #print data_vect
    #print weight_vect

    #calculate weights to result in values per promoter
    weight_input = []
    for i in range(len(data_vect)):
        weight_input.append([])
        for entry in data_vect[i]:
            weight_input[i].append(1.0/float(weight_vect[i]))
    
    #calculate proportions for each group
    proportions = []
    for i in range(len(data_vect)):
        if len(data_vect[i]) == 0:
            print i
            proportions.append(0.0)
        else:
            proportions.append(len(data_vect[i])*weight_input[i][0])

    #for i in range(len(proportions)):
    #    print proportions[-i-1]
    
    if not heatmap:
        n, bins, patches = P.hist(data_vect[::-1], 50, normed=False, weights = weight_input[::-1], histtype='step', label=sp_dv_names[::-1], stacked=False)
        P.legend(loc=1)
        #P.show()
        P.savefig(output_path + '/histograms/' + output_name + '_hist.pdf', format='pdf')
        P.close()
    
    if heatmap:
        n, bins, patches = P.hist(data_vect[::-1], 20, normed=False, weights = weight_input[::-1], histtype='step', label=sp_dv_names[::-1], stacked=False)
        P.legend(loc=1)
        #convert n into correct datastructure
        heatdata = np.array(np.ravel(n))
        max_density = max(heatdata)
        heatdata = np.reshape(heatdata, (len(n), len(bins)-1))
        print heatdata
        #output heatmap
        fig, ax = plt.subplots()
        htmap = ax.pcolor(heatdata, cmap = plt.cm.Blues)
        plt.setp(ax.get_yticklabels(), fontsize=font_size)
        #ax.set_xticks(np.arange(len(bins)-1), minor=False)
        ax.set_yticks(np.arange(len(sp_dv_names)), minor=False)
        ax.set_ybound(0, len(sp_dv_names))
        ax.set_yticklabels(sp_dv_names[::-1], minor=False)
        ax.set_xticklabels(((np.arange(5))*(max_range-min_range)/4+min_range).astype(int), minor=False)
        plt.title('max density: ' + str(max_density))
        P.savefig(output_path + '/heatmaps/' + output_name + '_heat.pdf', format='pdf')
        P.close()

def locsToRelative(input1, input2, output_name):
    #takes the data of the locations for 2 different TFs and calculates their locations of the 2nd relative to the 1st
    inFile1 = open(project_dir + 'outputs/whole_genome_scores_1000/' + input1, 'r')
    locs1 = []
    [locs1.append(line.split()) for line in inFile1]
    inFile.close()
    inFile2 = open(project_dir + 'outputs/whole_genome_scores_1000/' + input2, 'r')
    locs2 = []
    [locs2.append(line.split()) for line in inFile2]
    inFile2.close()
    output_lines = []
    for loc1_rec, loc2_rec in zip(locs1, locs2):
        current_output = []
        current_output.extend(loc1_rec[:2])
        if len(loc1_rec) > 2 and len(loc2_rec) > 2:
            loc1_locs = list(set(loc1_rec[2:]))
            loc2_locs = list(set(loc2_rec[2:]))
            for loc1 in loc1_locs:
                for loc2 in loc2_locs:
                    try: 
                        current_output.append(str(int(loc2)-int(loc1)))
                    except ValueError:
                        continue
        output_lines.append(current_output)
    outFile = open(project_dir + 'outputs/wg_rel_locs/'+output_name, 'w')
    for record in output_lines:
        outFile.write('\t'.join(record) + '\n')
    outFile.close()

def locsToRelativeRPsOnly(input1, input2, output_name):
    #takes the data of the locations for 2 different TFs and calculates their locations of the 2nd relative to the 1st
    #outputs these results only for the ribosomal protein genes
    inFile1 = open(project_dir + 'outputs/whole_genome_scores_1000/' + input1, 'r')
    locs1 = []
    [locs1.append(line.split()) for line in inFile1]
    inFile.close()
    inFile2 = open(project_dir + 'outputs/whole_genome_scores_1000/' + input2, 'r')
    locs2 = []
    [locs2.append(line.split()) for line in inFile2]
    inFile2.close()
    output_lines = []
    for loc1_rec, loc2_rec in zip(locs1, locs2):
        if loc1_rec[0] + ' ' + loc1_rec[1] in gene_list:
            current_output = []
            current_output.extend(loc1_rec[:2])
            if len(loc1_rec) > 2 and len(loc2_rec) > 2:
                loc1_locs = list(set(loc1_rec[2:]))
                loc2_locs = list(set(loc2_rec[2:]))
                for loc1 in loc1_locs:
                    for loc2 in loc2_locs:
                        try: 
                            current_output.append(str(int(loc2)-int(loc1)))
                        except ValueError:
                            continue
            output_lines.append(current_output)
    outFile = open(project_dir + 'outputs/relative_locations/relative_locations_1000/'+output_name, 'w')
    for record in output_lines:
        outFile.write('\t'.join(record) + '\n')
    outFile.close()

#analyses

#create a list of species divisioins for making graphs
species_divs = []
[species_divs.append([sp]) for sp in species_output_order]
sp_dv_names = species_output_order

#take locations of binding sites and calculate the distance between them
outputs = ['fhl1_scertf_6.0cutoff','fhl1_scertf_8.0cutoff','rap1_scertf_6.0cutoff','rap1_scertf_8.0cutoff','rrn7_Aade_6.0cutoff','rrn7_Aade_8.0cutoff','tbf1_Arub_6.0cutoff','tbf1_Arub_8.0cutoff']
for regulator in outputs:
    in_name = regulator + '.txt'
    out_name1 = 'mcm1_max_' + regulator + '.txt'
    locsToRelativeRPsOnly('mcm1_scertf__max_score.txt', in_name, out_name1)
    #output graphs of the relative locations of binding sites
    in_name = out_name1
    data_vector, weight_vector = simple_relative_locs(project_dir + 'outputs/relative_locations/relative_locations_1000/' + in_name ,species_divs)
    generateGraph(data_vector, weight_vector, -200, 200, project_dir + 'outputs/relative_locations/graphs/',in_name, heatmap=True, font_size=5)

#evolution of non-RP Rap1-Mcm1 sites
locsToRelative('Rap1ForwardMax_6.0', 'Mcm1_6.0', 'Rap1ForwardMax_6.0_Mcm1_6.0')
genes_by_ortho = loadYGOBorthos(gene_list)
output_values_YGOB_species(gene_list, genes_by_ortho)