
import sys
import copy
import gzip
from collections import OrderedDict
from collections import defaultdict
from collections import Counter

try:
    from .families import (
        FamilyInputError,
        iter_blast_clusters,
        iter_cdhit_clusters,
        iter_mmseqs_clusters,
        parse_blast,
        parse_cdhit,
        parse_mmseqs,
    )
except ImportError:
    from families import (
        FamilyInputError,
        iter_blast_clusters,
        iter_cdhit_clusters,
        iter_mmseqs_clusters,
        parse_blast,
        parse_cdhit,
        parse_mmseqs,
    )


def _open_text(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    return opener(path, 'rt')


def _iter_fasta_headers(fasta_file):
    if not fasta_file:
        return
    seen = set()
    with _open_text(fasta_file) as handle:
        for line in handle:
            if line.startswith('>'):
                identifier = line[1:].strip().split()[0]
                if identifier in seen:
                    raise FamilyInputError(f"Duplicate FASTA identifier {identifier!r} in {fasta_file}")
                seen.add(identifier)
                yield identifier


def _next_numeric_cluster_id(cluster_dict):
    if not cluster_dict:
        return 0
    numeric_ids = []
    for cluster_id in cluster_dict.keys():
        try:
            numeric_ids.append(int(cluster_id))
        except (TypeError, ValueError):
            continue
    if not numeric_ids:
        return len(cluster_dict)
    return max(numeric_ids) + 1


def _index_member_genomes(pangenome_clusters_First, pangenome_clusters_First_genomes,
                          pangenome_clusters_First_sequences):
    for cluster_id, sequences in pangenome_clusters_First_sequences.items():
        genomes = list(pangenome_clusters_First.get(cluster_id, []))
        for seq in sequences:
            pangenome_clusters_First_genomes[seq] = genomes


def backfill_missing_fasta_singletons(options, splitter, taxa_dict,
                                      pangenome_clusters_First,
                                      pangenome_clusters_First_genomes,
                                      pangenome_clusters_First_sequences,
                                      reps):
    """Ensure every sequence in the original FASTA belongs to a family.

    Edge-list cluster formats often omit nodes with no above-threshold edge.
    Partial mode requires the original FASTA, so use it as the authoritative
    set of input genes and synthesize singleton clusters for missing headers.
    """
    fasta_file = getattr(options, 'fasta', None)
    if not fasta_file:
        return 0

    assigned = set()
    for sequences in pangenome_clusters_First_sequences.values():
        assigned.update(sequences)

    missing_count = 0
    next_cluster_id = _next_numeric_cluster_id(pangenome_clusters_First)
    for header in _iter_fasta_headers(fasta_file):
        if header in assigned:
            continue
        cluster_id = str(next_cluster_id)
        next_cluster_id += 1
        taxa = header.split(splitter)[0]
        taxa_dict[taxa] += 1
        pangenome_clusters_First[cluster_id] = [taxa]
        pangenome_clusters_First_sequences[cluster_id] = [header]
        pangenome_clusters_First_genomes[header] = [taxa]
        reps[header] = [1, 1]
        assigned.add(header)
        missing_count += 1

    if missing_count:
        _index_member_genomes(pangenome_clusters_First, pangenome_clusters_First_genomes,
                              pangenome_clusters_First_sequences)
    return missing_count


def _legacy_cluster_result(parsed_clusters, options, splitter):
    """Adapt canonical parsed clusters to the historical five-dictionary API."""

    taxa_dict = defaultdict(int)
    cluster_taxa = OrderedDict()
    member_taxa = OrderedDict()
    cluster_sequences = OrderedDict()
    reps = OrderedDict()
    assigned = set()
    for index, parsed in enumerate(parsed_clusters):
        cluster_id = str(index)
        members = list(parsed.members)
        duplicate = assigned.intersection(members)
        if duplicate:
            sample = ', '.join(sorted(duplicate)[:10])
            raise FamilyInputError(f"Cluster members occur in multiple families: {sample}")
        assigned.update(members)
        taxa = list(OrderedDict.fromkeys(member.split(splitter, 1)[0] for member in members))
        cluster_taxa[cluster_id] = taxa
        cluster_sequences[cluster_id] = members
        reps[parsed.representative] = [len(members), len(taxa)]
        for member in members:
            taxon = member.split(splitter, 1)[0]
            taxa_dict[taxon] += 1
            member_taxa[member] = list(taxa)
    backfill_missing_fasta_singletons(
        options, splitter, taxa_dict, cluster_taxa, member_taxa, cluster_sequences, reps,
    )
    _index_member_genomes(cluster_taxa, member_taxa, cluster_sequences)
    return taxa_dict, cluster_taxa, member_taxa, cluster_sequences, reps

def cluster_CDHIT(options, splitter):
    return _legacy_cluster_result(iter_cdhit_clusters(options.clusters), options, splitter)

    First_in = open(options.clusters, 'r')
    # Read file contents so we can validate format before parsing.
    lines = First_in.read().splitlines()
    First_in.close()
    # Validate that this is a CD-HIT .clstr file. CD-HIT cluster files contain headers like:
    #   >Cluster 0
    # If the first '>' header looks like a FASTA header (e.g. '>genome_id...'), fail early with a helpful message.
    header = None
    for l in lines:
        if not l:
            continue
        if l.startswith('>'):
            header = l
            break
    if header is not None:
        parts = header.split(' ')
        if len(parts) < 2 or not parts[0].lower().startswith('>cluster'):
            sys.exit(f"Error: expected a CD-HIT .clstr formatted cluster file for '-clustering_format CD-HIT'.\n"
                     f"Found header: {header!r}\nPlease provide the CD-HIT .clstr file (usually ends with .clstr) not a FASTA file.")
    clusters = OrderedDict()
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_genomes = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    first = True
    taxa_dict = defaultdict(int)
    reps = OrderedDict()
    tmp_genomes = None
    ## Load in all data for easier reuse later
    for line in lines:
        if line.startswith('>'):
            if tmp_genomes != None:
                pangenome_clusters_First_genomes[rep] = tmp_genomes
            tmp_genomes = []
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
                    tmp_genomes.append(clustered_taxa)
                pangenome_clusters_First_sequences[cluster_id].append(clustered)

    if not first:
        pangenome_clusters_First_genomes[rep] = tmp_genomes
    backfill_missing_fasta_singletons(options, splitter, taxa_dict, pangenome_clusters_First,
                                      pangenome_clusters_First_genomes,
                                      pangenome_clusters_First_sequences, reps)
    _index_member_genomes(pangenome_clusters_First, pangenome_clusters_First_genomes,
                          pangenome_clusters_First_sequences)

    return taxa_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps

def cluster_BLAST(options, splitter):
    return _legacy_cluster_result(iter_blast_clusters(options.clusters, options=options), options, splitter)

    separator = '\t'
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_genomes = defaultdict(list)
    pangenome_clusters_First_sequences = defaultdict(list)
    taxa_dict = defaultdict(int)
    reps = OrderedDict()
    edges = defaultdict(set)
    seen_nodes = set()

    def register_node(node):
        edges[node]
        if node not in seen_nodes:
            taxa_dict[node.split(splitter)[0]] += 1
            seen_nodes.add(node)

    with open(options.clusters, 'r') as First_in:
        for line_number, line in enumerate(First_in, start=1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            elements = line.split(separator)
            if len(elements) < 2:
                sys.exit(f"Error parsing BLAST cluster file at line {line_number}: expected at least two tab-separated columns.")
            rep, child = elements[0], elements[1]
            register_node(rep)
            register_node(child)
            if rep == child:
                # Self-hit: still register the node so singletons are captured.
                continue
            edges[rep].add(child)
            edges[child].add(rep)

    visited = set()
    cluster_id = 0
    cluster_genomes_seen = {}  # cluster_id -> set of genomes, for O(1) membership

    def dfs(node, cluster_id):
        stack = [node]
        tmp_genomes = []
        genomes_seen = cluster_genomes_seen[cluster_id]
        while stack:
            current = stack.pop()
            if current not in visited:
                visited.add(current)
                clustered_taxa = current.split(splitter)[0]
                pangenome_clusters_First_sequences[cluster_id].append(current)
                if clustered_taxa not in genomes_seen:
                    genomes_seen.add(clustered_taxa)
                    pangenome_clusters_First[cluster_id].append(clustered_taxa)
                    tmp_genomes.append(clustered_taxa)
                for neighbor in edges[current]:
                    if neighbor not in visited:
                        stack.append(neighbor)

        pangenome_clusters_First_genomes[node] = tmp_genomes

    for node in edges:
        if node not in visited:
            pangenome_clusters_First[cluster_id] = []
            pangenome_clusters_First_sequences[cluster_id] = []
            pangenome_clusters_First_genomes[node] = []
            cluster_genomes_seen[cluster_id] = set()
            dfs(node, cluster_id)
            cluster_id += 1

    for cluster in pangenome_clusters_First:
        cluster_size = len(pangenome_clusters_First_sequences[cluster])
        reps[cluster] = [cluster_size, len(pangenome_clusters_First[cluster])]

    backfill_missing_fasta_singletons(options, splitter, taxa_dict, pangenome_clusters_First,
                                      pangenome_clusters_First_genomes,
                                      pangenome_clusters_First_sequences, reps)
    _index_member_genomes(pangenome_clusters_First, pangenome_clusters_First_genomes,
                          pangenome_clusters_First_sequences)

    return taxa_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps

def cluster_MMseqs(options,splitter):
    return _legacy_cluster_result(iter_mmseqs_clusters(options.clusters), options, splitter)

    separator = '\t'
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_genomes = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    taxa_dict = defaultdict(int)
    reps = OrderedDict()
    clusters_by_rep = OrderedDict()
    seen_nodes = set()

    def register_node(node):
        if node not in seen_nodes:
            taxa_dict[node.split(splitter)[0]] += 1
            seen_nodes.add(node)

    with open(options.clusters, 'r') as First_in:
        for line_number, line in enumerate(First_in, start=1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            elements = line.split(separator)
            if len(elements) < 2:
                sys.exit(f"Error parsing MMseqs cluster file at line {line_number}: expected at least two tab-separated columns.")
            rep, child = elements[0], elements[1]
            register_node(rep)
            register_node(child)
            clusters_by_rep.setdefault(rep, OrderedDict())
            clusters_by_rep[rep][rep] = None
            clusters_by_rep[rep][child] = None

    for cluster_id, (rep, members) in enumerate(clusters_by_rep.items()):
        cluster_key = str(cluster_id)
        pangenome_clusters_First[cluster_key] = []
        pangenome_clusters_First_sequences[cluster_key] = []
        genomes_seen = set()
        for child in members.keys():
            child_taxa = child.split(splitter)[0]
            if child_taxa not in genomes_seen:
                genomes_seen.add(child_taxa)
                pangenome_clusters_First[cluster_key].append(child_taxa)
            pangenome_clusters_First_sequences[cluster_key].append(child)
        reps[rep] = [len(pangenome_clusters_First_sequences[cluster_key]),
                     len(pangenome_clusters_First[cluster_key])]

    backfill_missing_fasta_singletons(options, splitter, taxa_dict, pangenome_clusters_First,
                                      pangenome_clusters_First_genomes,
                                      pangenome_clusters_First_sequences, reps)
    _index_member_genomes(pangenome_clusters_First, pangenome_clusters_First_genomes,
                          pangenome_clusters_First_sequences)

    return taxa_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps


def _canonical_combined_counting(options, pangenome_clusters_First, reps,
                                 combined_clustered, pangenome_clusters_First_genomes,
                                 splitter):
    cluster_ids = list(pangenome_clusters_First)
    rep_ids = list(reps)
    if len(cluster_ids) != len(rep_ids):
        raise FamilyInputError(
            "Initial cluster representatives do not align with initial families; cannot reconcile reclustering"
        )
    rep_to_cluster = dict(zip(rep_ids, cluster_ids))
    adjacency = {rep: set() for rep in rep_ids}
    second_by_rep = defaultdict(list)
    sequence_tag = getattr(options, 'sequence_tag', None) or ''
    for key, values in combined_clustered.items():
        if key not in adjacency:
            continue
        for value in values:
            if value in adjacency and (not sequence_tag or sequence_tag not in value):
                adjacency[key].add(value)
                adjacency[value].add(key)
            else:
                second_by_rep[key].append(value)

    component_for = {}
    for seed in rep_ids:
        if seed in component_for:
            continue
        stack = [seed]
        component = []
        while stack:
            rep = stack.pop()
            if rep in component_for:
                continue
            component_for[rep] = component
            component.append(rep)
            stack.extend(adjacency[rep])
        for rep in component:
            component_for[rep] = component

    result = OrderedDict()
    for cluster_id, rep in zip(cluster_ids, rep_ids):
        component = sorted(component_for[rep], key=lambda value: cluster_ids.index(rep_to_cluster[value]))
        source_clusters = [rep_to_cluster[value] for value in component]
        first_taxa = set()
        groupings = []
        for source_cluster in source_clusters:
            taxa = set(pangenome_clusters_First[source_cluster])
            first_taxa.update(taxa)
            groupings.append([f'{source_cluster}:{len(taxa)}'])
        seconds = []
        for component_rep in component:
            seconds.extend(second_by_rep.get(component_rep, []))
        seconds = list(OrderedDict.fromkeys(seconds))
        second_taxa = {member.split(splitter, 1)[0] for member in seconds}
        result[cluster_id] = [
            len(source_clusters),
            len(first_taxa),
            groupings,
            len(second_taxa - first_taxa),
            len(seconds),
            len(second_taxa),
        ]
    return result


#@profile
def combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered, pangenome_clusters_First_genomes, combined_pangenome_clusters_Second, combined_pangenome_clusters_Second_sequences, splitter):
    return _canonical_combined_counting(
        options,
        pangenome_clusters_First,
        reps,
        combined_pangenome_clusters_First_Second_clustered,
        pangenome_clusters_First_genomes,
        splitter,
    )

    num_clustered_First = defaultdict(list)
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    for cluster, First_genomes in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster
        Com_PEP_Genomes = 0
        Seconds = 0
        seen_Seconds = []
        added_Second_genomes = 0
        temp_pep_genomes = copy.deepcopy(First_genomes)
        try:  # get the cluster from the Second clusters which contains this rep
            clustered_combined = combined_pangenome_clusters_First_Second_clustered[rep]
            #We have to do this to correctly account for Seconds grouping multiple original First clusters
            for clust in clustered_combined:
                ### Get the original clustered genomes first:
                if options.sequence_tag not in clust:
                    original_clustered_genomes = pangenome_clusters_First_genomes[clust]
                    for genome in original_clustered_genomes:
                        if genome not in temp_pep_genomes:
                            temp_pep_genomes.append(genome)

            seen_clust_Genomes = []
            num_clustered_First[cluster].append(rep + '_' + str(len(First_genomes)))
            for clust in clustered_combined:
                if options.sequence_tag not in clust:  # Not good enough at the moment
                    clust_Genome = clust.split(splitter)[0]
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in First_genomes:
                            Com_PEP_Genomes += 1
                    num_clustered_First[cluster].append(clust + '_' + str(reps[clust][1]))
                elif options.sequence_tag in clust:
                    Seconds += 1
                    clust_Genome = clust.split(splitter)[0]
                    if clust_Genome not in seen_Seconds:
                        seen_Seconds.append(clust_Genome)
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in temp_pep_genomes:
                            added_Second_genomes += 1
                            temp_pep_genomes.append(clust_Genome)
                else:
                    sys.exit("Error: looking for sequence_tag")

            size_of_pep_clusters = []
            genomes = num_clustered_First[cluster]


            if len(genomes) > 1: #!!# So that we don't double count - This still needs to account for whether the same genome/genus is present however. Probably need to unique ti
               collecting_genomes = []
               for genome in genomes:
                   genome = genome.rsplit('_', 1)
                   collecting_genomes.append(pangenome_clusters_First[str(list_of_reps.index(genome[0]))])
                   size_of_pep_clusters.append([str(list_of_reps.index(genome[0])) + ':' + genome[1]])
               flattened_list = [item for sublist in collecting_genomes for item in sublist]
               element_counts = Counter(flattened_list)
               unique_elements = [element for element, count in element_counts.items() if count == 1]
               sum_size_of_pep_clusters = len(unique_elements)
            else:
                genome = genomes[0].rsplit('_', 1)
                size_of_pep_clusters.append([str(list_of_reps.index(genome[0]))+':'+genome[1]])
                sum_size_of_pep_clusters = int(genome[1])

            pangenome_clusters_Type[cluster] = [len(num_clustered_First[cluster]), sum_size_of_pep_clusters,
                                                size_of_pep_clusters, added_Second_genomes, Seconds, len(seen_Seconds)]

        except KeyError:
            ###Singleton
            num_First_genomes = [[str(cluster)+':'+str(len(First_genomes))]]
            pangenome_clusters_Type[cluster] = [1, len(First_genomes), num_First_genomes, added_Second_genomes, Seconds,
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
        rep = list_of_reps[int(cluster)]  # get the rep of the current cluster

        try:  # get the cluster from the storf clusters which contains this rep
            num_clustered_First[str(cluster)].append(str(rep) + '_' + str(len(First_taxa)))
            size_of_First_clusters = []
            Firsts = num_clustered_First[str(cluster)]
            for First in Firsts:
                First = First.rsplit('_', 1)
                size_of_First_clusters.append(int(First[1]))
                recorded_First.append(First[0])
            num_First_genomes = [[str(cluster) + ':' + str(len(First_taxa))]]
            pangenome_clusters_Type[cluster] = [len(num_clustered_First[cluster]), sum(size_of_First_clusters),
                                                num_First_genomes, 0, 0, 0]

        except KeyError:
            ###Singleton
            num_First_genomes = [[str(cluster)+':'+str(len(First_taxa))]]
            pangenome_clusters_Type[cluster] = [1, len(First_taxa), num_First_genomes, 0, 0, 0]

    # pangenome_clusters_Type = [Number of First clustered genomes or genera, Size of the cluster, Ditto, 0,0,0 ]
    return pangenome_clusters_Type


def _legacy_combined_result(parsed_clusters, options, splitter, taxa_dict=None):
    """Adapt second-stage clusters without row-order or self-hit assumptions."""

    combined_clustered = OrderedDict()
    not_second_only = set()
    second_taxa = OrderedDict()
    second_sequences = OrderedDict()
    sequence_tag = getattr(options, 'sequence_tag', None) or ''
    for index, parsed in enumerate(parsed_clusters):
        cluster_id = str(index)
        first = []
        second = []
        for member in parsed.members:
            taxon = member.split(splitter, 1)[0]
            if taxa_dict is not None:
                taxa_dict[taxon] += 1
            if sequence_tag and sequence_tag in member:
                second.append(member)
            else:
                first.append(member)
        if first:
            not_second_only.add(cluster_id)
            if len(first) + len(second) > 1:
                combined_clustered[first[0]] = first[1:] + second
        second_sequences[cluster_id] = second
        second_taxa[cluster_id] = list(OrderedDict.fromkeys(
            member.split(splitter, 1)[0] for member in second
        ))
    return combined_clustered, not_second_only, second_taxa, second_sequences



#@profile
def combined_clustering_CDHIT(options, taxa_dict, splitter):
    return _legacy_combined_result(iter_cdhit_clusters(options.reclustered), options, splitter, taxa_dict)

    Second_in = open(options.reclustered, 'r')
    combined_pangenome_clusters_First = OrderedDict()
    combined_pangenome_clusters_First_sequences = OrderedDict()
    combined_pangenome_clusters_Second = OrderedDict()
    combined_pangenome_clusters_Second_sequences = OrderedDict()
    combined_pangenome_clusters_First_Second_clustered = OrderedDict()

    not_Second_only_cluster_ids = set()
    already_seen_PEP = set()
    Combined_clusters = OrderedDict()
    Combined_reps = OrderedDict()
    first = True

    def finalize_cluster(current_cluster_id):
        if current_cluster_id is None:
            return
        cluster_size = len(Combined_clusters[current_cluster_id])
        Combined_reps.update({rep: cluster_size})
        for pep in combined_pangenome_clusters_First_sequences[current_cluster_id]:
            if pep != []:
                if pep in already_seen_PEP:
                    continue
                else:
                    already_seen_PEP.add(pep)
        if len(combined_pangenome_clusters_Second_sequences[current_cluster_id]) > 0 and len(combined_pangenome_clusters_First_sequences[current_cluster_id]) > 0:
            if len(combined_pangenome_clusters_First_sequences[current_cluster_id]) > 1:  # If we have clustered >1 First family, we need to record 1 as key and all others are val
                all_but_first = combined_pangenome_clusters_First_sequences[current_cluster_id][1:]
                storfs_clustered = combined_pangenome_clusters_Second_sequences[current_cluster_id]
                VALUE = all_but_first + storfs_clustered
            else:
                VALUE = combined_pangenome_clusters_Second_sequences[current_cluster_id]
            KEY = combined_pangenome_clusters_First_sequences[current_cluster_id][0]
            combined_pangenome_clusters_First_Second_clustered.update({KEY: VALUE})

    cluster_id = None
    rep = None
    for line in Second_in:
        if line.startswith('>'):
            if first == False:
                finalize_cluster(cluster_id)

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
                    not_Second_only_cluster_ids.add(cluster_id)
                    if clustered_taxa not in combined_pangenome_clusters_First[cluster_id]:
                        combined_pangenome_clusters_First[cluster_id].append(clustered_taxa)
                    combined_pangenome_clusters_First_sequences[cluster_id].append(clustered)

    if first == False:
        finalize_cluster(cluster_id)
    Second_in.close()

    return combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids, combined_pangenome_clusters_Second, combined_pangenome_clusters_Second_sequences




# def cluster_BLAST(options, splitter):
#     separator = '\t'
#     First_in = open(options.clusters, 'r')
#     pangenome_clusters_First = OrderedDict()
#     pangenome_clusters_First_genomes = defaultdict(list)
#     pangenome_clusters_First_sequences = defaultdict(list)
#     taxa_dict = defaultdict(int)
#     reps = OrderedDict()
#
#     for line in First_in:
#         elements = line.strip().split(separator)
#         rep, child = elements[0], elements[1]
#         child_taxa = child.split(splitter)[0]  # Extracting the genome identifier from the child sequence
#         # Counting occurrences of genomes
#         taxa_dict[child_taxa] += 1
#
#         if rep not in pangenome_clusters_First:
#             pangenome_clusters_First[rep] = []
#             pangenome_clusters_First_sequences[rep] = []
#
#         if child_taxa not in pangenome_clusters_First[rep]:
#             pangenome_clusters_First[rep].append(child_taxa)
#             pangenome_clusters_First_genomes[rep].append(child_taxa)
#
#         pangenome_clusters_First_sequences[rep].append(child)
#
#     for rep in pangenome_clusters_First:
#         cluster_size = len(pangenome_clusters_First_sequences[rep])
#         reps[rep] = [cluster_size, len(pangenome_clusters_First[rep])]
#
#     return taxa_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps




def combined_clustering_Edge_List(options, splitter):
    cluster_format = getattr(options, 'cluster_format', 'MMseqs')
    if cluster_format == 'CSV':
        parsed = iter_mmseqs_clusters(options.reclustered, delimiter=',')
    elif 'BLAST' in cluster_format:
        parsed = iter_blast_clusters(options.reclustered, options=options)
    else:
        parsed = iter_mmseqs_clusters(options.reclustered)
    return _legacy_combined_result(parsed, options, splitter)

    if options.cluster_format == 'CSV':
        separator = ','
    else:
        separator = '\t'

    cluster_id = 0
    last_rep = ''
    Second_in = open(options.reclustered, 'r')
    combined_pangenome_clusters_First = OrderedDict()
    combined_pangenome_clusters_First_sequences = OrderedDict()
    combined_pangenome_clusters_Second = OrderedDict()
    combined_pangenome_clusters_Second_sequences = OrderedDict()
    combined_pangenome_clusters_First_Second_clustered = OrderedDict()

    not_Second_only_cluster_ids = set()
    already_seen_PEP = set()
    Combined_clusters = OrderedDict()
    Combined_reps = OrderedDict()
    first = True

    def finalize_cluster(current_cluster_id, current_rep):
        if current_cluster_id is None:
            return
        cluster_size = len(Combined_clusters[str(current_cluster_id)])
        Combined_reps.update({current_rep: cluster_size})
        for pep in combined_pangenome_clusters_First_sequences[str(current_cluster_id)]:
            if pep != []:
                if pep in already_seen_PEP:
                    continue
                else:
                    already_seen_PEP.add(pep)
        if len(combined_pangenome_clusters_Second_sequences[str(current_cluster_id)]) > 0 and len(combined_pangenome_clusters_First_sequences[str(current_cluster_id)]) > 0:
            if len(combined_pangenome_clusters_First_sequences[str(current_cluster_id)]) > 1:  # If we have clustered >1 PEP family, we need to record 1 as key and all others are val
                all_but_first = combined_pangenome_clusters_First_sequences[str(current_cluster_id)][1:]
                storfs_clustered = combined_pangenome_clusters_Second_sequences[str(current_cluster_id)]
                VALUE = all_but_first + storfs_clustered
            else:
                VALUE = combined_pangenome_clusters_Second_sequences[str(current_cluster_id)]
            KEY = combined_pangenome_clusters_First_sequences[str(current_cluster_id)][0]
            combined_pangenome_clusters_First_Second_clustered.update({KEY: VALUE})

    for line in Second_in:
        elements = line.strip().split(separator)
        if len(elements) < 2:
            continue
        rep, child = elements[0], elements[1]
        child_taxa = child.split(splitter)[0]  # Extracting the genome identifier from the child sequence

        if first == True:
            Combined_clusters.update({str(cluster_id): []})
            combined_pangenome_clusters_First.update({str(cluster_id): []})
            combined_pangenome_clusters_First_sequences.update({str(cluster_id): []})
            combined_pangenome_clusters_Second.update({str(cluster_id): []})
            combined_pangenome_clusters_Second_sequences.update({str(cluster_id): []})
            Combined_reps.update({rep: 0})
            first = False

        if first == False:
            if rep != last_rep and last_rep != '':
                finalize_cluster(cluster_id, last_rep)

                cluster_id += 1
                Combined_clusters.update({str(cluster_id): []})
                combined_pangenome_clusters_First.update({str(cluster_id): []})
                combined_pangenome_clusters_First_sequences.update({str(cluster_id): []})
                combined_pangenome_clusters_Second.update({str(cluster_id): []})
                combined_pangenome_clusters_Second_sequences.update({str(cluster_id): []})
                Combined_reps.update({rep: 0})

        Combined_clusters[str(cluster_id)].append(child)
        if options.sequence_tag in line:
            if child_taxa not in combined_pangenome_clusters_Second[str(cluster_id)]:
                combined_pangenome_clusters_Second[str(cluster_id)].append(child_taxa)
            combined_pangenome_clusters_Second_sequences[str(cluster_id)].append(child)
        else:
            not_Second_only_cluster_ids.add(str(cluster_id))
            if child_taxa not in combined_pangenome_clusters_First[str(cluster_id)]:
                combined_pangenome_clusters_First[str(cluster_id)].append(child_taxa)
            combined_pangenome_clusters_First_sequences[str(cluster_id)].append(child)

        last_rep = rep

    if first == False:
        finalize_cluster(cluster_id, last_rep)
    Second_in.close()

    return combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids, combined_pangenome_clusters_Second, combined_pangenome_clusters_Second_sequences
