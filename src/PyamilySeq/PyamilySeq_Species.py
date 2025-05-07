
try:
    from .constants import *
    from .clusterings import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from clusterings import *
    from utils import *


def gene_presence_absence_output(options, genome_dict, pangenome_clusters_First_sorted, pangenome_clusters_First_sequences_sorted):
    print("Outputting gene_presence_absence file")
    output_dir = os.path.abspath(options.output_dir)
    #in_name = options.clusters.split('.')[0].split('/')[-1]
    gpa_outfile = os.path.join(output_dir, 'gene_presence_absence.csv')
    gpa_outfile = open(gpa_outfile, 'w')
    genome_dict = OrderedDict(sorted(genome_dict.items()))
    gpa_outfile.write('"Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment",'
                     '"Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc","')
    gpa_outfile.write('","'.join(genome_dict.keys()))
    gpa_outfile.write('"\n')
    for cluster, sequences in pangenome_clusters_First_sequences_sorted.items():
        average_sequences_per_genome = len(sequences) / len(pangenome_clusters_First_sorted[cluster])
        gpa_outfile.write('"group_'+str(cluster)+'","","","'+str(len(pangenome_clusters_First_sorted[cluster]))+'","'+str(len(sequences))+'","'+str(average_sequences_per_genome)+
                         '","","","","","","","",""')


        for genome in genome_dict.keys():
            full_out = ''
            tmp_list = []
            for value in sequences:
                if value.split('|')[0] == genome:
                    tmp_list.append(value.split('|')[1])
            if tmp_list:
                full_out += ',"'+'  '.join(tmp_list)+'"'
            else:
                full_out = ',""'
            gpa_outfile.write(full_out)
        gpa_outfile.write('\n')

### Below is some unfinished code
    # edge_list_outfile = open(in_name+'_edge_list.csv','w')
    # for cluster, sequences in pangenome_clusters_First_sequences_sorted.items():
    #     output = []
    #     for entry in sequences:
    #         # Split each entry at '|'
    #         genome, gene = entry.split('|')
    #         # Format the result as "gene  genome"
    #         output.append(f"{gene}\t{genome}")
    #     for line in output:
    #         edge_list_outfile.write(line + '\n')




def get_cores(options,genome_dict):
    ##Calculate core groups
    groups = OrderedDict()
    cores = OrderedDict()
    prev_top = len(genome_dict)
    first = True
    for group in options.species_groups.split(','):
        calculated_floor = math.floor(int(group) / 100 * len(genome_dict))
        if first == False:
            groups[group] = (calculated_floor,prev_top)
        else:
            groups[group] = (calculated_floor, prev_top)
            first = False
        prev_top = calculated_floor
        first_core_group = 'First_core_' + group
        cores[first_core_group] = []
        if options.reclustered != None:
            extended_core_group = 'extended_core_' + group
            cores[extended_core_group] = []
            combined_core_group = 'combined_core_' + group
            cores[combined_core_group] = []
            second_core_group = 'Second_core_' + group
            cores[second_core_group] = []
            only_second_core_group = 'only_Second_core_' + group
            cores[only_second_core_group] = []
    return cores, groups

#@profile
def calc_First_only_core(cluster, First_num, groups, cores):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= int(First_num) <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['First_core_'+family_group].append(cluster)

#@profile
def calc_single_First_extended_Second_only_core(cluster, First_num, groups, cores, Second_num): # Count gene families extended with StORFs
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= First_num+Second_num <= fir):
        res = idx

    family_group = list(groups)[res]
    cores['extended_core_' + family_group].append(cluster)


#@profile
def calc_multi_First_extended_Second_only_core(cluster, First_num, groups, cores, Second_num): # Count seperately those gene families extended with StORF_Reporter but combined >1 PEP
    groups_as_list = list(groups.values())
    # Looping through the list to find the matching condition
    for idx, (sec, fir) in enumerate(groups_as_list):
        if sec <= First_num + Second_num <= fir:
            res = idx
            break
    family_group = list(groups)[res]
    cores['combined_core_' + family_group].append(cluster)


#@profile
def calc_Second_only_core(cluster, Second_num, groups, cores):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= Second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['Second_core_' + family_group].append(cluster)

#@profile
def calc_only_Second_only_core(cluster, Second_num, groups, cores): # only count the true storf onlies
    try:
        groups_as_list = list(groups.values())
        for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= Second_num <= fir):
            res = idx
        family_group = list(groups)[res]
        cores['only_Second_core_' + family_group].append(cluster)
    except UnboundLocalError:
        sys.exit("Error in calc_only_Second_only_core")


#@profile
def cluster(options):

    if options.cluster_format == 'CD-HIT':
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps = cluster_CDHIT(options, '|')
    elif 'BLAST' in options.cluster_format:
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps = cluster_BLAST(options, '|')
    elif 'MMseqs' in options.cluster_format:
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps = cluster_MMseqs(options, '|')

    ###
    cores, groups = get_cores(options, genome_dict)
    ###

    if options.reclustered != None: #FIX
        if options.cluster_format == 'CD-HIT':
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second, combined_pangenome_clusters_Second_sequences = combined_clustering_CDHIT(options, genome_dict, '|')
        elif 'TSV' in options.cluster_format or 'CSV' in options.cluster_format:
            #Fix
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second,combined_pangenome_clusters_Second_sequences  = combined_clustering_Edge_List(options, '|')
        pangenome_clusters_Type = combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered, pangenome_clusters_First_genomes, '|')
    else:
        pangenome_clusters_Type = single_clustering_counting(pangenome_clusters_First, reps)



    Number_Of_Second_Extending_But_Same_Genomes = 0

    sorted_first_keys = sort_keys_by_values(pangenome_clusters_First, pangenome_clusters_First_sequences)
    pangenome_clusters_First_sorted = reorder_dict_by_keys(pangenome_clusters_First, sorted_first_keys)
    pangenome_clusters_First_sequences_sorted = reorder_dict_by_keys(pangenome_clusters_First_sequences, sorted_first_keys)
    pangenome_clusters_Type_sorted = reorder_dict_by_keys(pangenome_clusters_Type, sorted_first_keys)

    print("Calculating Groups")
    seen_groupings = []
    for cluster, numbers in pangenome_clusters_Type_sorted.items():
    ############################### Calculate First only
        cluster = str(cluster)
        for grouping in numbers[2]: #!!# Could do with a more elegant solution
            current_cluster = grouping[0].split(':')[0]
            if current_cluster not in seen_groupings:
                seen_groupings.append(current_cluster)
                current_cluster_size = grouping[0].split(':')[1]
                calc_First_only_core(current_cluster, current_cluster_size,groups,cores)
            ############################# Calculate First and Reclustered-Second
                if numbers[0] == 1 and numbers[3] >= 1:  # If Seconds did not combine First reps
                    calc_single_First_extended_Second_only_core(cluster, numbers[1], groups, cores, numbers[3])
                elif numbers[0] > 1 and numbers[3] >= 1:  # If unique Seconds combined multiple Firsts
                    calc_multi_First_extended_Second_only_core(cluster, numbers[1], groups, cores, numbers[3])
                elif numbers[4] >= 1:
                    Number_Of_Second_Extending_But_Same_Genomes += 1
            else:
                if options.verbose == True:
                    print("First cluster " + current_cluster + " already processed - This is likely because it was clustered with another First representative.")

    if options.reclustered != None:
        combined_pangenome_clusters_ONLY_Second_Type = defaultdict(list)
        combined_pangenome_clusters_Second_Type = defaultdict(list)
        for cluster, genomes in combined_pangenome_clusters_Second.items():
            if cluster in not_Second_only_cluster_ids:
                combined_pangenome_clusters_Second_Type[cluster] = [cluster, len(genomes)]
            else:
                combined_pangenome_clusters_ONLY_Second_Type[cluster] = [cluster, len(genomes)]
        for cluster, data in combined_pangenome_clusters_Second_Type.items():
            if data[1] >= 1:
                calc_Second_only_core(cluster, data[1], groups, cores)
        for cluster, data in combined_pangenome_clusters_ONLY_Second_Type.items():
            if data[1] >= 1:
                calc_only_Second_only_core(cluster, data[1], groups, cores)
    ###########################
    ### Output
    output_path = os.path.abspath(options.output_dir)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    stats_out = os.path.join(output_path,'summary_statistics.txt')
    key_order = ['First_core_', 'extended_core_', 'combined_core_', 'Second_core_','only_Second_core_']
    with open(stats_out, 'w') as outfile:
        print("Number of Genomes: " + str(len(genome_dict)))
        outfile.write("Number of Genomes: " + str(len(genome_dict)) + "\n")
        print("Gene Groups:")
        outfile.write("Gene Groups\n")
        for key_prefix in key_order:
            for key, value in cores.items():
                if key.startswith(key_prefix):
                    print(f"{key}: {len(value)}")
                    outfile.write(f"{key}: {len(value)}\n")
        print("Total Number of First Gene Groups (Including Singletons): " + str(len(pangenome_clusters_First_sequences_sorted)))
        outfile.write("Total Number of First Gene Groups (Including Singletons): " + str(len(pangenome_clusters_First_sequences_sorted)))
        if options.reclustered!= None:
            print("Total Number of Second Gene Groups (Including Singletons): " + str(
                len(combined_pangenome_clusters_Second_sequences)))
            print("Total Number of First Gene Groups That Had Additional Second Sequences But Not New Genomes: " + str(
                Number_Of_Second_Extending_But_Same_Genomes))
            outfile.write("\nTotal Number of Second Gene Groups (Including Singletons): " + str(
                len(combined_pangenome_clusters_Second_sequences)))
            outfile.write("\nTotal Number of First Gene Groups That Had Additional Second Sequences But Not New Genomes: " + str(
                Number_Of_Second_Extending_But_Same_Genomes))
        #Report number of first and second clusters and do the ame for genus
    if options.gene_presence_absence_out != False:
        gene_presence_absence_output(options,genome_dict, pangenome_clusters_First_sorted, pangenome_clusters_First_sequences_sorted)


    ###Need to fix this below. If full/partial the ifs need to be different. If full we first need to output the gfs then align. if -wruite-groups not presented then it needs
    # to be done for alignment full anyway...

    genome_list = list(genome_dict.keys())
    if options.run_mode == 'Full':
        sequences = read_fasta(options.fasta)
        if options.reclustered == None:
            combined_pangenome_clusters_Second_sequences = None
        ## Output representative sequences
        representatives_out = os.path.join(output_path,'pan_genome_reference.fa')
        with open(representatives_out, 'w') as outfile:
            for cluster, ids in pangenome_clusters_First_sequences.items():
                outfile.write('>group_'+str(cluster)+'\n')
                wrapped_aa_seq = wrap_sequence(sequences[ids[0]], 60)
                outfile.write(wrapped_aa_seq+'\n')
        if options.write_groups != None:
            print("Outputting gene group FASTA files")
            #output_dir = os.path.dirname(os.path.abspath(options.output_dir))
            output_dir = os.path.join(options.output_dir, 'Gene_Groups_Output')
            write_groups_func(options,output_dir, key_order, cores, sequences,
                         pangenome_clusters_First_sequences_sorted, combined_pangenome_clusters_Second_sequences)

            if options.align_core != None:
                print("Processing gene group alignment")
                process_gene_groups(options, output_dir, None, None, genome_list, 'core_gene_alignment.aln')

    elif options.run_mode == 'Partial':
        sequences = read_fasta(options.fasta)
        if options.reclustered == None:
            combined_pangenome_clusters_Second_sequences = None
        # else: ## Output representative sequences - Under development
        #     representatives_out = os.path.join(output_path, 'pan_genome_reference_reclustered.fa')
        #     with open(representatives_out, 'w') as outfile:
        #         for cluster, ids in combined_pangenome_clusters_Second_sequences.items():
        #             outfile.write('>group_' + str(cluster) + '\n')
        #             try:
        #                 wrapped_aa_seq = wrap_sequence(sequences[ids[0]], 60)
        #             except:
        #                 print(2)
        #             outfile.write(wrapped_aa_seq + '\n')
        ## Output representative sequences
        representatives_out = os.path.join(output_path,'pan_genome_reference.fa')
        with open(representatives_out, 'w') as outfile:
            for cluster, ids in pangenome_clusters_First_sequences.items():
                outfile.write('>group_'+str(cluster)+'\n')
                wrapped_aa_seq = wrap_sequence(sequences[ids[0]], 60)
                outfile.write(wrapped_aa_seq+'\n')
        if options.write_groups != None:
            print("Outputting gene group FASTA files")
            output_dir = os.path.join(options.output_dir, 'Gene_Groups_Output')
            write_groups_func(options,output_dir, key_order, cores, sequences,
                         pangenome_clusters_First_sequences_sorted, combined_pangenome_clusters_Second_sequences)

            if options.align_core != None:
                print("Processing gene group alignment")
                process_gene_groups(options, output_dir, None, None, genome_list, 'core_gene_alignment.aln')



        #
        # if options.align_core != None:
        #     #output_dir = os.path.dirname(os.path.abspath(options.output_dir))
        #     output_dir = os.path.join(options.output_dir, 'Gene_Families_Output')
        #     if not os.path.exists(output_dir):
        #         os.makedirs(output_dir)
        #     process_gene_families(options, output_dir, 'concatenated_genes_aligned.fasta')

    #
    # elif options.run_mode == 'Partial':
    #     if options.align_core != None and options.fasta != None and options.write_groups != None:
    #         process_gene_families(options, os.path.join(output_dir, 'Gene_Families_Output'), 'concatenated_genes_aligned.fasta')
    #
    #
    #
    #
    #

