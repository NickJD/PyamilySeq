import subprocess
import shutil
import os
import glob
import sys
import copy
from collections import OrderedDict
from collections import defaultdict

def cluster_CDHIT(options, splitter):
    First_in = open(options.clusters, 'r')
    clusters = OrderedDict()
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    first = True
    taxa_dict = defaultdict(int)
    reps = OrderedDict()
    ## Load in all data for easier reuse later
    for line in First_in:
        if '>Cluster 7575' in line:
            print()
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
            taxa = clustered.split(splitter)[0]
            taxa_dict[taxa] += 1
            if '*' in line:
                rep = clustered
                reps.update({rep: [0, 0]})
            if first == False:
                clusters[cluster_id].append(clustered)
                clustered_taxa = clustered.split(splitter)[0]
                if clustered_taxa not in pangenome_clusters_First[cluster_id]:
                    pangenome_clusters_First[cluster_id].append(clustered_taxa)
                pangenome_clusters_First_sequences[cluster_id].append(clustered)
    return taxa_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps



#@profile
def combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered, splitter):
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
                    clust_Genome = clust.split(splitter)[0]
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in pep_genomes:
                            Com_PEP_Genomes += 1
                    num_clustered_First[cluster].append(clust + '_' + str(reps[clust][1]))
                elif options.sequence_tag in clust:
                    Seconds += 1
                    clust_Genome = clust.split(splitter)[0]
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
    # pangenome_clusters_Type = [Number of First clustered genomes or genera, Size of the cluster, Ditto, Added Seconds,Number of Seconds,Unique Seconds ]
    return pangenome_clusters_Type


#@profile
def single_clustering_counting(pangenome_clusters_First, reps):
    num_clustered_First = defaultdict(list)
    recorded_First = []
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    for cluster, First_taxa in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster

        try:  # get the cluster from the storf clusters which contains this rep
            num_clustered_First[cluster].append(rep + '_' + str(len(First_taxa)))
            size_of_First_clusters = []
            Firsts = num_clustered_First[cluster]
            for First in Firsts:
                First = First.rsplit('_', 1)
                size_of_First_clusters.append(int(First[1]))
                recorded_First.append(First[0])
            pangenome_clusters_Type[cluster] = [len(num_clustered_First[cluster]), sum(size_of_First_clusters),
                                                size_of_First_clusters, 0, 0, 0]

        except KeyError:
            ###Singleton
            num_First_taxa = [len(First_taxa)]
            pangenome_clusters_Type[cluster] = [1, len(First_taxa), num_First_taxa, 0, 0, 0]

    # pangenome_clusters_Type = [Number of First clustered genomes or genera, Size of the cluster, Ditto, 0,0,0 ]
    return pangenome_clusters_Type



#@profile
def combined_clustering_CDHIT(options, taxa_dict, splitter):
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
            genome = clustered.split(splitter)[0]
            taxa_dict[genome] += 1
            if '*' in line:
                rep = clustered
                Combined_reps.update({rep: 0})
            if first == False:
                Combined_clusters[cluster_id].append(clustered)
                clustered_taxa = clustered.split(splitter)[0]
                if options.sequence_tag in line:
                    if clustered_taxa not in combined_pangenome_clusters_Second[cluster_id]:
                        combined_pangenome_clusters_Second[cluster_id].append(clustered_taxa)
                    combined_pangenome_clusters_Second_sequences[cluster_id].append(clustered)
                else:
                    if cluster_id not in not_Second_only_cluster_ids:
                        not_Second_only_cluster_ids.append(cluster_id)  # Tell us which StORF_Reporter clustered are unmatched to a PEP
                    if clustered_taxa not in combined_pangenome_clusters_First[cluster_id]:
                        combined_pangenome_clusters_First[cluster_id].append(clustered_taxa)
                    combined_pangenome_clusters_First_sequences[cluster_id].append(clustered)


    return combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids, combined_pangenome_clusters_Second


def cluster_EdgeList(options,splitter):
    if options.cluster_format == 'TSV':
        separator = '\t'
    elif options.cluster_format == 'CSV':
        separator = ','
    cluster_id = 0
    last_rep = ''
    first = True
    First_in = open(options.clusters, 'r')
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    taxa_dict = defaultdict(int)
    reps = OrderedDict()
    for line in First_in:
        rep, child = line.strip().split(separator)
        child_taxa = child.split(splitter)[0]  # Extracting the genome identifier from the child sequence
        # Counting occurrences of genomes
        taxa_dict[child_taxa] += 1
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
        if child_taxa not in pangenome_clusters_First[cluster_id]:
            pangenome_clusters_First[cluster_id].append(child_taxa)

        pangenome_clusters_First_sequences[cluster_id].append(child)
        last_rep = rep
        cluster_size = len(pangenome_clusters_First_sequences[cluster_id])
        reps.update({rep: [cluster_size, len(pangenome_clusters_First[cluster_id])]})


    return taxa_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps


def combined_clustering_Edge_List(options, splitter):
    if options.cluster_format == 'TSV':
        separator = '\t'
    elif options.cluster_format == 'CSV':
        separator = ','

    cluster_id = 0
    last_rep = ''
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
        rep, child = line.strip().split(separator)
        child_taxa = child.split(splitter)[0]  # Extracting the genome identifier from the child sequence

        if first == True:
            Combined_clusters.update({cluster_id: []})
            combined_pangenome_clusters_First.update({cluster_id: []})
            combined_pangenome_clusters_First_sequences.update({cluster_id: []})
            combined_pangenome_clusters_Second.update({cluster_id: []})
            combined_pangenome_clusters_Second_sequences.update({cluster_id: []})
            Combined_reps.update({rep: 0})
            first = False

        if first == False:
            if rep != last_rep and last_rep != '':
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

                cluster_id += 1
                Combined_clusters.update({cluster_id: []})
                combined_pangenome_clusters_First.update({cluster_id: []})
                combined_pangenome_clusters_First_sequences.update({cluster_id: []})
                combined_pangenome_clusters_Second.update({cluster_id: []})
                combined_pangenome_clusters_Second_sequences.update({cluster_id: []})
                Combined_reps.update({rep: 0})

        Combined_clusters[cluster_id].append(child)
        if options.sequence_tag in line:
            if child_taxa not in combined_pangenome_clusters_Second[cluster_id]:
                combined_pangenome_clusters_Second[cluster_id].append(child_taxa)
            combined_pangenome_clusters_Second_sequences[cluster_id].append(child)
        else:
            if cluster_id not in not_Second_only_cluster_ids:
                not_Second_only_cluster_ids.append(cluster_id)  # Tell us which StORF_Reporter clustered are unmatched to a PEP
            if child_taxa not in combined_pangenome_clusters_First[cluster_id]:
                combined_pangenome_clusters_First[cluster_id].append(child_taxa)
            combined_pangenome_clusters_First_sequences[cluster_id].append(child)

        last_rep = rep

    return combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids, combined_pangenome_clusters_Second
