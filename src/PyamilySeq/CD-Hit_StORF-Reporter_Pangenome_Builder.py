from collections import OrderedDict,defaultdict
import copy
import math
import sys
import argparse


try:
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *


def combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered):
    ####################
    ALL_PEPS = 0
    Com_PEPs = 0
    Com_PEPs_List = []
    total_combined_pep_clusters = []
    only_PEP = []
    num_clustered_PEP = defaultdict(list)
    recorded_PEP = []
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    pep_max_comb = [0, 0]
    for cluster, pep_genomes in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster
        # if rep not in recorded_PEP:
        recorded_PEP.append(rep)
        Com_PEP_Genomes = 0
        StORFs = 0
        seen_StORFs = []
        Added_StORF_Genomes = 0
        try:  # get the cluster from the storf clusters which contains this rep
            clustered_combined = combined_pangenome_clusters_First_Second_clustered[
                rep]  # Not true clusters - I put a PEP as key myself
            seen_clust_Genomes = []
            num_clustered_PEP[cluster].append(rep + '_' + str(len(pep_genomes)))

            for clust in clustered_combined:
                if options.reclustered not in clust:  # Not good enough at the moment
                    ### Need to get the number of pep genomes for each pep clustered into this
                    Com_PEPs += 1  #
                    recorded_PEP.append(clust)
                    if (rep + '_' + str(len(pep_genomes))) not in Com_PEPs_List:
                        Com_PEPs_List.append(rep + '_' + str(len(pep_genomes)))
                    clust_Genome = clust.split('|')[0]
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in pep_genomes:
                            Com_PEP_Genomes += 1
                    num_clustered_PEP[cluster].append(clust + '_' + str(reps[clust][1]))
                elif options.reclustered in clust:
                    StORFs += 1
                    clust_Genome = clust.split('|')[0]
                    if clust_Genome not in seen_StORFs:
                        seen_StORFs.append(clust_Genome)
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in pep_genomes:
                            Added_StORF_Genomes += 1
                else:
                    print("WHAT")

            size_of_pep_clusters = []
            peps = num_clustered_PEP[cluster]
            if len(peps) > 1:
                total_combined_pep_clusters.append(peps)
            for pep in peps:
                pep = pep.rsplit('_', 1)
                # if pep[0] not in recorded_PEP: # Do not record PEPs are combined and on their own.
                size_of_pep_clusters.append(int(pep[1]))
                ALL_PEPS += int(pep[1])
                recorded_PEP.append(pep[0])
                # else:
                #     print("W")
            pangenome_clusters_Type[cluster] = [len(num_clustered_PEP[cluster]), sum(size_of_pep_clusters),
                                                size_of_pep_clusters, Added_StORF_Genomes, StORFs, len(seen_StORFs)]
            if len(peps) > 1:
                print("")
                if Added_StORF_Genomes > 10:
                    print("")
            if len(peps) > 2:
                print("")
            if len(peps) > pep_max_comb[0]:
                pep_max_comb = [len(peps), peps]
            num_pep_combined = sum(size_of_pep_clusters)
            if len(peps) > 1 and all(i >= 5 for i in size_of_pep_clusters):
                print("")
        except KeyError:
            ###Singleton
            num_pep_genomes = [len(pep_genomes)]
            ALL_PEPS += len(pep_genomes)
            pangenome_clusters_Type[cluster] = [1, len(pep_genomes), num_pep_genomes, Added_StORF_Genomes, StORFs,
                                                len(seen_StORFs)]
            only_PEP.append(cluster)

    return pangenome_clusters_Type

def single_clustering_counting(options, pangenome_clusters_First, reps):
    ####################
    ALL_PEPS = 0
    Com_PEPs = 0
    Com_PEPs_List = []
    total_combined_pep_clusters = []
    only_PEP = []
    num_clustered_PEP = defaultdict(list)
    recorded_PEP = []
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    pep_max_comb = [0, 0]
    for cluster, pep_genomes in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster
        # if rep not in recorded_PEP:
        recorded_PEP.append(rep)

        try:  # get the cluster from the storf clusters which contains this rep
            # clustered_combined = combined_pangenome_clusters_First_Second_clustered[
            #     rep]  # Not true clusters - I put a PEP as key myself
            # seen_clust_Genomes = []
            num_clustered_PEP[cluster].append(rep + '_' + str(len(pep_genomes)))

            # for clust in clustered_combined:
            #     if options.reclustered not in clust:  # Not good enough at the moment
            #         ### Need to get the number of pep genomes for each pep clustered into this
            #         Com_PEPs += 1  #
            #         recorded_PEP.append(clust)
            #         if (rep + '_' + str(len(pep_genomes))) not in Com_PEPs_List:
            #             Com_PEPs_List.append(rep + '_' + str(len(pep_genomes)))
            #         clust_Genome = clust.split('|')[0]
            #         if clust_Genome not in seen_clust_Genomes:
            #             seen_clust_Genomes.append(clust_Genome)
            #             if clust_Genome not in pep_genomes:
            #                 Com_PEP_Genomes += 1
            #         num_clustered_PEP[cluster].append(clust + '_' + str(reps[clust][1]))
            #     elif options.reclustered in clust:
            #         StORFs += 1
            #         clust_Genome = clust.split('|')[0]
            #         if clust_Genome not in seen_StORFs:
            #             seen_StORFs.append(clust_Genome)
            #         if clust_Genome not in seen_clust_Genomes:
            #             seen_clust_Genomes.append(clust_Genome)
            #             if clust_Genome not in pep_genomes:
            #                 Added_StORF_Genomes += 1
            #     else:
            #         print("WHAT")

            size_of_pep_clusters = []
            peps = num_clustered_PEP[cluster]
            if len(peps) > 1:
                total_combined_pep_clusters.append(peps)
            for pep in peps:
                pep = pep.rsplit('_', 1)
                # if pep[0] not in recorded_PEP: # Do not record PEPs are combined and on their own.
                size_of_pep_clusters.append(int(pep[1]))
                ALL_PEPS += int(pep[1])
                recorded_PEP.append(pep[0])
                # else:
                #     print("W")
            pangenome_clusters_Type[cluster] = [len(num_clustered_PEP[cluster]), sum(size_of_pep_clusters),
                                                size_of_pep_clusters, 0, 0, 0]



        except KeyError:
            ###Singleton
            num_pep_genomes = [len(pep_genomes)]
            ALL_PEPS += len(pep_genomes)
            pangenome_clusters_Type[cluster] = [1, len(pep_genomes), num_pep_genomes, 0, 0,
                                                0]
            only_PEP.append(cluster)

    return pangenome_clusters_Type




def combined_clustering(options, genome_dict):
    unique_genomes = []
    Second_in = open(options.reclustered, 'r')
    combined_pangenome_clusters_First = OrderedDict()
    combined_pangenome_clusters_First_sequences = OrderedDict()
    combined_pangenome_clusters_Second = OrderedDict()
    combined_pangenome_clusters_Second_sequences = OrderedDict()
    combined_pangenome_clusters_First_Second_clustered = OrderedDict()

    clusters_with_Seconds = []

    not_StORF_Only_Cluster_IDs = []
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
                #  if len(Combined_pangenome_clusters_PEP[cluster_id]) > 1:
                #      print("Here")
                if len(combined_pangenome_clusters_Second_sequences[cluster_id]) > 0 and len(
                        combined_pangenome_clusters_First_sequences[cluster_id]) > 0:
                    if len(combined_pangenome_clusters_First_sequences[
                               cluster_id]) > 1:  # If we have clustered >1 PEP family, we need to record 1 as key and all others are val
                        all_but_first = combined_pangenome_clusters_First_sequences[cluster_id][1:]
                        storfs_clustered = combined_pangenome_clusters_Second_sequences[cluster_id]
                        VALUE = all_but_first + storfs_clustered
                    else:
                        VALUE = combined_pangenome_clusters_Second_sequences[cluster_id]
                    KEY = combined_pangenome_clusters_First_sequences[cluster_id][0]
                    combined_pangenome_clusters_First_Second_clustered.update({KEY: VALUE})
                if len(Combined_clusters[cluster_id]) == 1 : # Stop at end of file
                    # Combined_pangenome_clusters_PEP.popitem()
                    # Combined_pangenome_clusters_PEP_SEQS.popitem()
                    # Combined_pangenome_clusters_StORF_SEQS.popitem()
                    # Combined_pangenome_clusters_Con_StORF.popitem()
                    # Combined_reps.popitem()
                    break
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
            if genome not in unique_genomes:
                unique_genomes.append(genome)
            genome_dict[genome] += 1
            if '*' in line:
                rep = clustered
                Combined_reps.update({rep: 0})
            if first == False:
                Combined_clusters[cluster_id].append(clustered)
                clustered_genome = clustered.split('|')[0]
                if options.reclustered in line:
                    if cluster_id not in clusters_with_Seconds:  # For counting?
                        clusters_with_Seconds.append(cluster_id)
                    if clustered_genome not in combined_pangenome_clusters_Second[cluster_id]:
                        combined_pangenome_clusters_Second[cluster_id].append(clustered_genome)
                    combined_pangenome_clusters_Second_sequences[cluster_id].append(clustered)
                else:
                    if cluster_id not in not_StORF_Only_Cluster_IDs:
                        not_StORF_Only_Cluster_IDs.append(
                            cluster_id)  # Tell us which StORF_Reporter clustered are unmatched to a PEP
                    if clustered_genome not in combined_pangenome_clusters_First[cluster_id]:
                        combined_pangenome_clusters_First[cluster_id].append(clustered_genome)
                    combined_pangenome_clusters_First_sequences[cluster_id].append(clustered)


    return combined_pangenome_clusters_First_Second_clustered,not_StORF_Only_Cluster_IDs, combined_pangenome_clusters_Second, unique_genomes


def cluster(options):
    First_in = open(options.clusters, 'r')
    clusters = OrderedDict()
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()

    count = 0
    first = True
    genome_dict = defaultdict(int)
    reps = OrderedDict()
    county = 0

    singleton_cluster = []
    ## Load in all data for easier reuse later
    for line in First_in:
        if line.startswith('>'):
            if first == False:
                cluster_size = len(clusters[cluster_id])
                reps.update({rep:[cluster_size,len(pangenome_clusters_First[cluster_id])]})
                if len(clusters[cluster_id]) == 1:# and not singleton_cluster: # Stop at clusters smaller than 10
                    #singleton_cluster.append(cluster_id)
                #     pangenome_clusters_PEP.popitem()
                #     pangenome_clusters_PEP_SEQS.popitem()

                #     reps.popitem()
                    break
            Ensem_genomes, Con_genomes = [], []
            cluster_id = line.strip('>')
            cluster_id = cluster_id.strip('\n')
            cluster_id = cluster_id.split(' ')[1]
            clusters.update({cluster_id: []})
            pangenome_clusters_First.update({cluster_id:[]})
            pangenome_clusters_First_sequences.update({cluster_id:[]})

            first = False
        else:
            clustered = line.split('\t')[1]
            clustered = clustered.split('>')[1]
            clustered = clustered.split('...')[0]
            genome = clustered.split('|')[0]
            genome_dict[genome] +=1
            if '*' in line:
                rep = clustered
                reps.update({rep:[0,0]})
            if first == False:
                clusters[cluster_id].append(clustered)
                clustered_genome = clustered.split('|')[0]

                if clustered_genome not in pangenome_clusters_First[cluster_id]:
                    pangenome_clusters_First[cluster_id].append(clustered_genome)
                pangenome_clusters_First_sequences[cluster_id].append(clustered)

    ######################################

    if options.reclustered != None:
        combined_pangenome_clusters_First_Second_clustered,not_StORF_Only_Cluster_IDs,combined_pangenome_clusters_Second,\
            unique_genomes = combined_clustering(options, genome_dict)
        pangenome_clusters_Type = combined_clustering_counting()


    else:
        pangenome_clusters_Type = single_clustering_counting(options, pangenome_clusters_First, reps)



    #######################################

    core_99 = math.floor(9.9/10 * len(genome_dict))
    core_95 = math.floor(9.5/10 * len(genome_dict))
    core_90 = math.floor(9/10 * len(genome_dict))
    core_15 = math.floor(1.5/10 * len(genome_dict))

    cores = OrderedDict({'first_core_99':0,'first_core_95':0,'first_core_15':0,'extended_99':0,'extended_95':0
             ,'extended_15':0,'comb_extended_99':0,'comb_extended_95':0,'comb_extended_15':0,'second_core_99':0,'second_core_95':0,'second_core_15':0,
                                     'only_second_core_99':0,'only_second_core_95':0,'only_second_core_15':0})



    StORF_Seqs_Extended = []
    StORF_Genomes_Extended = []

    record_all_pep_15 = []

    core_list = []
    soft_core_list = []
    accessory_list = []

    second_core_only = []

    list_all_pep_numbers = defaultdict(int)


    ############################ Count PEP separately first to get TRUE Ensembl gene families
    def calc_pep_only_core(pep_num):
        list_all_pep_numbers[pep_num] +=1
        if pep_num >= math.floor(core_99):# and StORF_num == 0:
            cores['first_core_99'] += 1
        elif pep_num >= math.floor(core_95) and pep_num < math.floor(core_99):# and StORF_num == 0:
            cores['first_core_95'] += 1
        # elif pep_num >= math.floor(core_90) and pep_num < math.floor(core_95):# and StORF_num == 0:
        #     cores['first_core_90'] += 1
        if pep_num >= math.floor(core_15) and pep_num < math.floor(core_95):# and StORF_num == 0:  # this catch captures some from first_core_90
            cores['first_core_15'] += 1
            record_all_pep_15.append(pep_num)
        #####################
    def calc_single_pep_extended_StORF_only_core(cluster,pep_num,storf_num): # Count gene families extended with StORFs
        if pep_num < math.floor(core_99) and pep_num != 0 and pep_num+storf_num >= math.floor(core_99):
            cores['extended_99'] +=1
            core_list.append(cluster)
        elif pep_num < math.floor(core_95) and pep_num != 0 and pep_num+storf_num >= math.floor(core_95) and pep_num+storf_num < math.floor(core_99):
            cores['extended_95'] +=1
            soft_core_list.append(cluster)
        # elif pep_num < math.floor(core_90) and pep_num != 0 and pep_num+storf_num >= math.floor(core_90) and pep_num+storf_num < math.floor(core_95):
        #     cores['extended_90'] +=1
        if pep_num < math.floor(core_15) and pep_num != 0 and pep_num+storf_num >= math.floor(core_15) and pep_num+storf_num < math.floor(core_95):
            cores['extended_15'] +=1
            accessory_list.append(cluster)
    #####################################
    def calc_multi_pep_extended_StORF_only_core(pep_num,storf_num): # Count seperately those gene families extended with StORF_Reporter but combined >1 PEP
        if pep_num < math.floor(core_99) and pep_num != 0 and pep_num+storf_num >= math.floor(core_99):
            cores['comb_extended_99'] +=1
        elif pep_num < math.floor(core_95) and pep_num != 0 and pep_num+storf_num >= math.floor(core_95) and pep_num+storf_num < math.floor(core_99):
            cores['comb_extended_95'] +=1
        # elif pep_num < math.floor(core_90) and pep_num != 0 and pep_num+storf_num >= math.floor(core_90) and pep_num+storf_num < math.floor(core_95):
        #     cores['comb_extended_90'] +=1
        if pep_num < math.floor(core_15) and pep_num != 0 and pep_num+storf_num >= math.floor(core_15) and pep_num+storf_num < math.floor(core_95):
            cores['comb_extended_15'] +=1
    ######################### StORFs Only >>><<<
    def calc_StORF_only_core(storf_num):
        if storf_num >= math.floor(core_99):# and StORF_num == 0:
            cores['second_core_99'] += 1
        elif storf_num >= math.floor(core_95) and storf_num < math.floor(core_99):# and StORF_num == 0:
            cores['second_core_95'] += 1
        # elif storf_num >= math.floor(core_90) and storf_num < math.floor(core_95):# and StORF_num == 0:
        #     cores['second_core_90'] += 1
        if storf_num >= math.floor(core_15) and storf_num < math.floor(core_95):# and StORF_num == 0:  # this catch captures some from first_core_90
            cores['second_core_15'] += 1
    ###########################
    def calc_only_StORF_only_core(cluster,storf_num): # only count the true storf onlies
        if storf_num >= math.floor(core_99):# and StORF_num == 0:
            cores['only_second_core_99'] += 1
            second_core_only.append(cluster)
        elif storf_num >= math.floor(core_95) and storf_num < math.floor(core_99):# and StORF_num == 0:
            cores['only_second_core_95'] += 1
        # elif storf_num >= math.floor(core_90) and storf_num < math.floor(core_95):# and StORF_num == 0:
        #     cores['only_second_core_90'] += 1
        if storf_num >= math.floor(core_15) and storf_num < math.floor(core_95):# and StORF_num == 0:  # this catch captures some from first_core_90
            cores['only_second_core_15'] += 1

    record_all_pep = []

    counter = 0

    Number_Of_StORF_Extending_But_Same_Genomes = 0


    print("Running")
    for cluster, numbers in pangenome_clusters_Type.items(): # put limits here to make sure storf and enembl only are in more than one genome.
        if numbers[3] >= 1:
            StORF_Genomes_Extended.append(numbers[3])
        if numbers[4] >= 1:
            StORF_Seqs_Extended.append(numbers[4])
    ############################### Calc PEP only
        ######### TO fix - Only loop through the first 1's to get the baseline pep numbs?
        if numbers[0] == 1 and numbers[1] >=2: # If StORFs did not combine PEP reps
            calc_pep_only_core(numbers[1])#,numbers[3])
            counter +=1
        elif numbers[0] >1 and numbers[1] >=2: # IF StORFs combined multiple PEP
            calc_pep_only_core(numbers[2][0])
            counter += 1

    ############################# Calc PEP and StORF_Reporter - M
        if numbers[0] == 1 and numbers[3] >= 1: # If StORFs did not combine PEP reps
            calc_single_pep_extended_StORF_only_core(cluster,numbers[1],numbers[3])
        elif numbers[0] >1 and numbers[3] >= 1: # IF unique StORFs combined multiple PEP
            #grouped_pep = sum(numbers[2])
            #for num in numbers[2]:
            calc_multi_pep_extended_StORF_only_core(numbers[1],numbers[3])

        elif numbers[4] >= 1:
            Number_Of_StORF_Extending_But_Same_Genomes +=1
            # for num in numbers[2]:

            ############## Typing for the StORF_Reporter-Data
    if options.reclustered != None:
        combined_pangenome_clusters_ONLY_Second_Type = defaultdict(list)
        combined_pangenome_clusters_Second_Type = defaultdict(list)
        for cluster, genomes in combined_pangenome_clusters_Second.items():
            if cluster in not_StORF_Only_Cluster_IDs:
                combined_pangenome_clusters_Second_Type[cluster] = [cluster, len(genomes)]
            else:
                combined_pangenome_clusters_ONLY_Second_Type[cluster] = [cluster, len(genomes)]
                #     multi_PEP_Combined_By_StORFs_num_of_PEP_Clusters +=1
        for cluster, data in combined_pangenome_clusters_Second_Type.items():
            # if data[1] >= 2:
            calc_StORF_only_core(data[1])  # ,numbers[3])multi_PEP_Combined_By_StORFs
        for cluster, data in combined_pangenome_clusters_ONLY_Second_Type.items():
            if data[1] >= 2:

                calc_only_StORF_only_core(cluster, data[1])  # ,numbers[3])
    ###########################
    ###########################
    print("End")
    print(cores)



def main():

    parser = argparse.ArgumentParser(description='PyamilySeq ' + PyamilySeq_Version + ': PyamilySeq Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-c', action='store', dest='clusters', help='Clustering output file from CD-HIT, DIAMOND or MMseqs2',
                        required=True)
    required.add_argument('-f', action='store', dest='format', choices=['CD-HIT', 'DIAMOND', 'MMseqs2'],
                        help='Which clustering algorithm used (CD-HIT, DIAMOND or MMseqs2)', required=True)


    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rc', action='store', dest='reclustered', help='Clustering output file from secondary round of clustering',
                        required=False)
    # optional.add_argument('-min_score', action='store', dest='minscore', default='30', type=int,
    #                     help='Minimum BitScore to keep StORF: Default 30')
    #
    # output = parser.add_argument_group('Output')
    # output.add_argument('-oname', action="store", dest='o_name', required=False,
    #                     help='Default - Appends \'_UR\' to end of input GFF filename')
    # output.add_argument('-odir', action="store", dest='o_dir', required=False,
    #                     help='Default -  Same directory as input GFF')
    # output.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
    #                     help='Default - False: Output as .gz')


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

    if options.format != None:
        options.format = 'StORF'

    cluster(options)

    print("Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq\n"
          "Please report any issues to: https://github.com/NickJD/PyamilySeq/issues\n#####")



if __name__ == "__main__":
    main()
    print("Complete")

