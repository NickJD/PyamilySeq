import argparse
from collections import OrderedDict
from collections import defaultdict

try:
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from utils import *


def categorise_percentage(percent):
    """Categorise the percentage of genomes with multicopy genes."""
    if 20 <= percent < 40:
        return "20-40%"
    elif 40 <= percent < 60:
        return "40-60%"
    elif 60 <= percent < 80:
        return "60-80%"
    elif 80 <= percent < 95:
        return "80-95%"
    elif 95 <= percent < 99:
        return "95-99%"
    elif 99 <= percent <= 100:
        return "99-100%"
    return None

# Read cd-hit .clstr file and extract information
def read_cd_hit_output(clustering_output):
    clusters = OrderedDict()

    with open(clustering_output, 'r') as f:
        current_cluster_id = None

        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster_id = line.split(' ')[1]
                clusters[current_cluster_id] = []
            elif line and current_cluster_id is not None:
                parts = line.split('\t')
                if len(parts) > 1:
                    clustered_info = parts[1]
                    length = clustered_info.split(',')[0]
                    length = int(''.join(c for c in length if c.isdigit()))
                    clustered_header = clustered_info.split('>')[1].split('...')[0]
                    clustered_header = '>' + clustered_header

                    if 'at ' in clustered_info and '%' in clustered_info.split('at ')[-1]:
                        percent_identity = extract_identity(clustered_info)
                    elif line.endswith('*'):
                        percent_identity = 100.0
                    else:
                        raise ValueError("Percent identity not found in the string.")

                    clusters[current_cluster_id].append({
                        'header': clustered_header,
                        'length': length,
                        'percent_identity': percent_identity
                    })

    return clusters


# Summarise the information for each cluster
def summarise_clusters(options,clusters, output):
    multicopy_groups = defaultdict(int)  # Counter for groups with multicopy genes

    with open(output, 'w') as out_f:
        out_f.write("Cluster_ID\tNum_Sequences\tAvg_Length\tLength_Range\tAvg_Identity\tIdentity_Range\n")

        for cluster_id, seqs in clusters.items():
            num_seqs = len(seqs)
            lengths = [seq['length'] for seq in seqs]
            identities = [seq['percent_identity'] for seq in seqs]

            avg_length = sum(lengths) / num_seqs if num_seqs > 0 else 0
            length_range = f"{min(lengths)}-{max(lengths)}" if num_seqs > 0 else "N/A"

            avg_identity = sum(identities) / num_seqs if num_seqs > 0 else 0
            identity_range = f"{min(identities):.2f}-{max(identities):.2f}" if num_seqs > 0 else "N/A"

            out_f.write(
                f"{cluster_id}\t{num_seqs}\t{avg_length:.2f}\t{length_range}\t{avg_identity:.2f}\t{identity_range}\n")

            # Count genomes with more than one gene
            genome_to_gene_count = defaultdict(int)
            for seq in seqs:
                genome = seq['header'].split('|')[0].replace('>','')
                genome_to_gene_count[genome] += 1

            num_genomes_with_multiple_genes = sum(1 for count in genome_to_gene_count.values() if count > 1)

            # Calculate the percentage of genomes with multicopy genes

            multicopy_percentage = (num_genomes_with_multiple_genes / options.genome_num) * 100
            category = categorise_percentage(multicopy_percentage)
            if category:
                multicopy_groups[category] += 1

        # Define the order of categories for printout
        category_order = ["20-40%", "40-60%", "60-80%", "80-95%", "95-99%", "99-100%"]

        # Print the number of clusters with multicopy genes in each percentage range, in the correct order
        for category in category_order:
            print(f"Number of clusters with multicopy genes in {category} range: {multicopy_groups[category]}")


# Main function to parse arguments and run the analysis
def main():
    parser = argparse.ArgumentParser(description='PyamilySeq ' + PyamilySeq_Version + ': Cluster-Summary - A tool to summarise CD-HIT clustering files.')
    ### Required Arguments
    required = parser.add_argument_group('Required Parameters')
    required.add_argument('-input_clstr', action="store", dest="input_clstr",
                          help='Input CD-HIT .clstr file',
                          required=True)
    required.add_argument('-output', action="store", dest="output",
                          help="Output TSV file to store cluster summaries - Will add '.tsv' if not provided by user",
                          required=True)
    required.add_argument('-genome_num', action='store', dest='genome_num', type=int,
                          help='The total number of genomes must be provide',
                          required=True)
    #required.add_argument("-clustering_format", action="store", dest="clustering_format", choices=['CD-HIT','TSV','CSV'],
    #                      help="Clustering format to use: CD-HIT or TSV (MMseqs2, BLAST, DIAMOND) / CSV edge-list file (Node1\tNode2).",
    #                      required=True)

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-output_dir', action="store", dest="output_dir",
                          help='Default: Same as input file',
                          required=False)

    misc = parser.add_argument_group("Misc Parameters")
    misc.add_argument("-verbose", action="store_true", dest="verbose",
                      help="Print verbose output.",
                      required=False)
    misc.add_argument("-v", "--version", action="version",
                      version=f"PyamilySeq: Group-Summary version {PyamilySeq_Version} - Exiting",
                      help="Print out version number and exit")


    options = parser.parse_args()
    print("Running PyamilySeq " + PyamilySeq_Version+ ": Group-Summary ")

    ### File handling
    options.input_clstr = fix_path(options.input_clstr)
    if options.output_dir is None:
        options.output_dir = os.path.dirname(os.path.abspath(options.input_clstr))
    output_path = os.path.abspath(options.output_dir)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_name = options.output
    if not output_name.endswith('.tsv'):
        output_name += '.tsv'
    output_file_path = os.path.join(output_path, output_name)
    ###

    clusters = read_cd_hit_output(options.input_clstr)
    summarise_clusters(options,clusters, output_file_path)


if __name__ == "__main__":
    main()
