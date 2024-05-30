#from line_profiler_pycharm import profile

from collections import OrderedDict,defaultdict
import copy
import math
import sys
import argparse

try:
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *


def coinfinder_output(options):
    print("placeholder")


def get_cores(options,genome_dict):
    ##Calculate core groups
    groups = OrderedDict()
    cores = OrderedDict()
    prev_top = len(genome_dict)
    first = True
    for group in options.core_groups.split(','):
        calculated_floor = math.floor(int(group) / 100 * len(genome_dict))
        if first == False:
            groups[group] = (calculated_floor,prev_top -1)
        else:
            groups[group] = (calculated_floor, prev_top)
            first = False
        prev_top = calculated_floor
        first_core_group = 'first_core_' + group
        cores[first_core_group] = 0
        if options.reclustered != None:
            extended_core_group = 'extended_core_' + group
            cores[extended_core_group] = 0
            combined_core_group = 'combined_core_' + group
            cores[combined_core_group] = 0
            second_core_group = 'second_core_' + group
            cores[second_core_group] = 0
            only_second_core_group = 'only_second_core_' + group
            cores[only_second_core_group] = 0
    return cores, groups

#@profile
def calc_pep_only_core(pep_num, groups, cores):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= pep_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['first_core_'+family_group] +=1

#@profile
def calc_single_pep_extended_StORF_only_core(pep_num, groups, cores, second_num): # Count gene families extended with StORFs
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= pep_num+second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['extended_core_' + family_group] += 1


#@profile
def calc_multi_pep_extended_StORF_only_core(pep_num, groups, cores, second_num): # Count seperately those gene families extended with StORF_Reporter but combined >1 PEP
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= pep_num+second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['combined_core_' + family_group] += 1


#@profile
def calc_StORF_only_core(groups, cores, second_num):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['second_core_' + family_group] += 1

#@profile
def calc_only_StORF_only_core(groups, cores, second_num): # only count the true storf onlies
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['only_second_core_' + family_group] += 1





#@profile
def combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered):
    num_clustered_First = defaultdict(list)
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    for cluster, pep_genomes in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster
        Com_PEP_Genomes = 0
        Seconds = 0
        seen_Seconds = []
        added_Second_genomes = 0
        try:  # get the cluster from the storf clusters which contains this rep
            clustered_combined = combined_pangenome_clusters_First_Second_clustered[rep]  # Not true clusters - I put a PEP as key myself
            seen_clust_Genomes = []
            num_clustered_First[cluster].append(rep + '_' + str(len(pep_genomes)))
            for clust in clustered_combined:
                if options.sequence_tag not in clust:  # Not good enough at the moment
                    clust_Genome = clust.split('|')[0]
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in pep_genomes:
                            Com_PEP_Genomes += 1
                    num_clustered_First[cluster].append(clust + '_' + str(reps[clust][1]))
                elif options.sequence_tag in clust:
                    Seconds += 1
                    clust_Genome = clust.split('|')[0]
                    if clust_Genome not in seen_Seconds:
                        seen_Seconds.append(clust_Genome)
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in pep_genomes:
                            added_Second_genomes += 1
                else:
                    sys.exit("Error: looking for sequence_tag")

            size_of_pep_clusters = []
            peps = num_clustered_First[cluster]
            for pep in peps:
                pep = pep.rsplit('_', 1)
                size_of_pep_clusters.append(int(pep[1]))
            pangenome_clusters_Type[cluster] = [len(num_clustered_First[cluster]), sum(size_of_pep_clusters),
                                                size_of_pep_clusters, added_Second_genomes, Seconds, len(seen_Seconds)]

        except KeyError:
            ###Singleton
            num_pep_genomes = [len(pep_genomes)]
            pangenome_clusters_Type[cluster] = [1, len(pep_genomes), num_pep_genomes, added_Second_genomes, Seconds,
                                                len(seen_Seconds)]

    return pangenome_clusters_Type

#@profile
def single_clustering_counting(options, pangenome_clusters_First, reps):
    num_clustered_PEP = defaultdict(list)
    recorded_PEP = []
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    for cluster, pep_genomes in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster

        try:  # get the cluster from the storf clusters which contains this rep
            num_clustered_PEP[cluster].append(rep + '_' + str(len(pep_genomes)))
            size_of_pep_clusters = []
            peps = num_clustered_PEP[cluster]
            for pep in peps:
                pep = pep.rsplit('_', 1)
                size_of_pep_clusters.append(int(pep[1]))
                recorded_PEP.append(pep[0])
            pangenome_clusters_Type[cluster] = [len(num_clustered_PEP[cluster]), sum(size_of_pep_clusters),
                                                size_of_pep_clusters, 0, 0, 0]

        except KeyError:
            ###Singleton
            num_pep_genomes = [len(pep_genomes)]
            pangenome_clusters_Type[cluster] = [1, len(pep_genomes), num_pep_genomes, 0, 0, 0]

    return pangenome_clusters_Type



#@profile
def combined_clustering_CDHIT(options, genome_dict):
    unique_genomes = []
    Second_in = open(options.reclustered, 'r')
    combined_pangenome_clusters_First = OrderedDict()
    combined_pangenome_clusters_First_sequences = OrderedDict()
    combined_pangenome_clusters_Second = OrderedDict()
    combined_pangenome_clusters_Second_sequences = OrderedDict()
    combined_pangenome_clusters_First_Second_clustered = OrderedDict()

    not_Second_only_cluster_ids = []
    already_seen_PEP = []
    Combined_clusters = OrderedDict()
    Combined_reps = OrderedDict()
    first = True
    for line in Second_in:
        if line.startswith('>'):
            if first == False:
                cluster_size = len(Combined_clusters[cluster_id])
                Combined_reps.update({rep: cluster_size})
                for pep in combined_pangenome_clusters_First_sequences[cluster_id]:
                    if pep != []:
                        if pep in already_seen_PEP:
                            continue
                        else:
                            already_seen_PEP.append(pep)
                if len(combined_pangenome_clusters_Second_sequences[cluster_id]) > 0 and len(combined_pangenome_clusters_First_sequences[cluster_id]) > 0:
                    if len(combined_pangenome_clusters_First_sequences[cluster_id]) > 1:  # If we have clustered >1 PEP family, we need to record 1 as key and all others are val
                        all_but_first = combined_pangenome_clusters_First_sequences[cluster_id][1:]
                        storfs_clustered = combined_pangenome_clusters_Second_sequences[cluster_id]
                        VALUE = all_but_first + storfs_clustered
                    else:
                        VALUE = combined_pangenome_clusters_Second_sequences[cluster_id]
                    KEY = combined_pangenome_clusters_First_sequences[cluster_id][0]
                    combined_pangenome_clusters_First_Second_clustered.update({KEY: VALUE})
            cluster_id = line.strip('>')
            cluster_id = cluster_id.strip('\n')
            cluster_id = cluster_id.split(' ')[1]
            Combined_clusters.update({cluster_id: []})
            combined_pangenome_clusters_First.update({cluster_id: []})
            combined_pangenome_clusters_First_sequences.update({cluster_id: []})
            combined_pangenome_clusters_Second.update({cluster_id: []})
            combined_pangenome_clusters_Second_sequences.update({cluster_id: []})

            first = False
        else:
            clustered = line.split('\t')[1]
            clustered = clustered.split('>')[1]
            clustered = clustered.split('...')[0]
            genome = clustered.split('|')[0]
            genome_dict[genome] += 1
            if '*' in line:
                rep = clustered
                Combined_reps.update({rep: 0})
            if first == False:
                Combined_clusters[cluster_id].append(clustered)
                clustered_genome = clustered.split('|')[0]
                if options.sequence_tag in line:
                    if clustered_genome not in combined_pangenome_clusters_Second[cluster_id]:
                        combined_pangenome_clusters_Second[cluster_id].append(clustered_genome)
                    combined_pangenome_clusters_Second_sequences[cluster_id].append(clustered)
                else:
                    if cluster_id not in not_Second_only_cluster_ids:
                        not_Second_only_cluster_ids.append(cluster_id)  # Tell us which StORF_Reporter clustered are unmatched to a PEP
                    if clustered_genome not in combined_pangenome_clusters_First[cluster_id]:
                        combined_pangenome_clusters_First[cluster_id].append(clustered_genome)
                    combined_pangenome_clusters_First_sequences[cluster_id].append(clustered)


    return combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids, combined_pangenome_clusters_Second, unique_genomes

def cluster_MMSEQS2(options):
    cluster_id = 0
    last_rep = ''
    first = True
    First_in = open(options.clusters, 'r')
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    genome_dict = defaultdict(int)
    reps = OrderedDict()
    for line in First_in:
        rep, child = line.strip().split('\t')
        child_genome = child.split('|')[0]  # Extracting the genome identifier from the child sequence

        # Counting occurrences of genomes
        genome_dict[child_genome] += 1

        if first == True:
            pangenome_clusters_First[0] = []
            pangenome_clusters_First_sequences[0] = []
            first = False

        if rep != last_rep and last_rep != '':
            cluster_id +=1
            pangenome_clusters_First[cluster_id] = []
            pangenome_clusters_First_sequences[cluster_id] = []
            cluster_size = len(pangenome_clusters_First_sequences[cluster_id-1])
            reps.update({last_rep: [cluster_size, len(pangenome_clusters_First[cluster_id-1])]})
            pangenome_clusters_First[cluster_id] = []
            pangenome_clusters_First_sequences[cluster_id] = []

        if child_genome not in pangenome_clusters_First[cluster_id]:
            pangenome_clusters_First[cluster_id].append(child_genome)

        pangenome_clusters_First_sequences[cluster_id].append(child)
        last_rep = rep
        cluster_size = len(pangenome_clusters_First_sequences[cluster_id])
        reps.update({rep: [cluster_size, len(pangenome_clusters_First[cluster_id])]})


    return genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps



def cluster_CDHIT(options):
    First_in = open(options.clusters, 'r')
    clusters = OrderedDict()
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    first = True
    genome_dict = defaultdict(int)
    reps = OrderedDict()
    ## Load in all data for easier reuse later
    for line in First_in:
        if line.startswith('>'):
            if first == False:
                cluster_size = len(clusters[cluster_id])
                reps.update({rep: [cluster_size, len(pangenome_clusters_First[cluster_id])]})
            cluster_id = line.strip('>')
            cluster_id = cluster_id.strip('\n')
            cluster_id = cluster_id.split(' ')[1]
            clusters.update({cluster_id: []})
            pangenome_clusters_First.update({cluster_id: []})
            pangenome_clusters_First_sequences.update({cluster_id: []})

            first = False
        else:
            clustered = line.split('\t')[1]
            clustered = clustered.split('>')[1]
            clustered = clustered.split('...')[0]
            genome = clustered.split('|')[0]
            genome_dict[genome] += 1
            if '*' in line:
                rep = clustered
                reps.update({rep: [0, 0]})
            if first == False:
                clusters[cluster_id].append(clustered)
                clustered_genome = clustered.split('|')[0]
                if clustered_genome not in pangenome_clusters_First[cluster_id]:
                    pangenome_clusters_First[cluster_id].append(clustered_genome)
                pangenome_clusters_First_sequences[cluster_id].append(clustered)
    return genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps

#@profile
def cluster(options):

    if options.format == 'CD-HIT':
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps = cluster_CDHIT(options)
    elif options.format == 'MMseqs2':
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps = cluster_MMSEQS2(options)

    ######################################
    cores, groups = get_cores(options, genome_dict)
    ###

    if options.reclustered != None:
        if options.format == 'CD-HIT':
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second,\
                unique_genomes = combined_clustering_CDHIT(options, genome_dict)
        pangenome_clusters_Type = combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered)
    else:
        pangenome_clusters_Type = single_clustering_counting(options, pangenome_clusters_First, reps)


    counter = 0
    Number_Of_StORF_Extending_But_Same_Genomes = 0


    print("Running")
    for cluster, numbers in pangenome_clusters_Type.items():
    ############################### Calc PEP only
        ######### TO fix - Only loop through the first 1's to get the baseline pep numbs?
        if numbers[0] == 1 and numbers[1] >=2: # If StORFs did not combine PEP reps
            calc_pep_only_core(numbers[1],groups,cores)#,numbers[3])
            counter +=1
        elif numbers[0] >1 and numbers[1] >=2: # IF StORFs combined multiple PEP
            calc_pep_only_core(numbers[2][0],groups,cores)
            counter += 1

    ############################# Calc PEP and StORF_Reporter - M
        if numbers[0] == 1 and numbers[3] >= 1: # If StORFs did not combine PEP reps
            calc_single_pep_extended_StORF_only_core(numbers[1],groups, cores, numbers[3])
        elif numbers[0] >1 and numbers[3] >= 1: # IF unique StORFs combined multiple PEP
            #grouped_pep = sum(numbers[2])
            #for num in numbers[2]:
            calc_multi_pep_extended_StORF_only_core(numbers[1], groups, cores, numbers[3])

        elif numbers[4] >= 1:
            Number_Of_StORF_Extending_But_Same_Genomes +=1
            # for num in numbers[2]:

            ############## Typing for the StORF_Reporter-Data
    if options.reclustered != None:
        combined_pangenome_clusters_ONLY_Second_Type = defaultdict(list)
        combined_pangenome_clusters_Second_Type = defaultdict(list)
        for cluster, genomes in combined_pangenome_clusters_Second.items():
            if cluster in not_Second_only_cluster_ids:
                combined_pangenome_clusters_Second_Type[cluster] = [cluster, len(genomes)]
            else:
                combined_pangenome_clusters_ONLY_Second_Type[cluster] = [cluster, len(genomes)]
                #     multi_PEP_Combined_By_StORFs_num_of_PEP_Clusters +=1
        for cluster, data in combined_pangenome_clusters_Second_Type.items():
            # if data[1] >= 2:
            calc_StORF_only_core(groups, cores, data[1])  # ,numbers[3])multi_PEP_Combined_By_StORFs
        for cluster, data in combined_pangenome_clusters_ONLY_Second_Type.items():
            if data[1] >= 2:

                calc_only_StORF_only_core(groups, cores, data[1])  # ,numbers[3])
    ###########################
    print("End")
    key_order = ['first_core_', 'extended_core_', 'combined_core_', 'second_core_','only_second_core_']
    print("Gene Family Groups:")
    for key_prefix in key_order:
        for key, value in cores.items():
            if key.startswith(key_prefix):
                print(f"{key}: {value}")

    if options.coinfinder_out != None:
        print("")


def main():

    parser = argparse.ArgumentParser(description='PyamilySeq ' + PyamilySeq_Version + ': PyamilySeq Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-c', action='store', dest='clusters', help='Clustering output file from CD-HIT, DIAMOND or MMseqs2',
                        required=True)
    required.add_argument('-f', action='store', dest='format', choices=['CD-HIT', 'DIAMOND', 'MMseqs2'],
                        help='Default - "CD-HIT": Which clustering algorithm used (CD-HIT, DIAMOND or MMseqs2)', required=True)


    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rc', action='store', dest='reclustered', help='Clustering output file from secondary round of clustering',
                        required=False)
    optional.add_argument('-st', action='store', dest='sequence_tag', help='Default - "StORF": Unique identifier to be used to distinguish the second of two rounds of clustered sequences',
                        required=False)
    optional.add_argument('-groups', action="store", dest='core_groups', default="99,80,15",
                        help='Default - (\'99,95,90,80,15\'): Gene family groups to use')
    optional.add_argument('-cf', action='store', dest='coinfinder_out', help='Default - No Coinfinder output: File name to be used for Coinfinder data output',
                        required=False)

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-verbose', action='store', dest='verbose', default=False, type=eval, choices=[True, False],
                        help='Default - False: Print out runtime messages')
    misc.add_argument('-v', action='store_true', dest='version',
                        help='Default - False: Print out version number and exit')


    options = parser.parse_args()
    if options.clusters == None or options.format == None:
        if options.version:
            sys.exit(PyamilySeq_Version)
        else:
            exit('PyamilySeq: error: the following arguments are required: -c, -f')

    if options.sequence_tag == None:
        options.sequence_tag = 'StORF'

    options.core_groups = options.core_groups + ',0'

    cluster(options)

    print("Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq\n"
          "Please report any issues to: https://github.com/NickJD/PyamilySeq/issues\n#####")





if __name__ == "__main__":
    main()
    print("Complete")

