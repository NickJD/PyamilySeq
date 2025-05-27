import subprocess
import shutil
import os
import glob
import collections
from tempfile import NamedTemporaryFile
import sys
import re
import math

####
# Placeholder for the distance function
levenshtein_distance_cal = None
# Check for Levenshtein library once
try:
    import Levenshtein as LV
    # Assign the optimised function
    def levenshtein_distance_calc(seq1, seq2):
        return LV.distance(seq1, seq2)
except (ModuleNotFoundError, ImportError):
    print("Levenshtein package not installed - Will fallback to slower Python implementation.")
    # Fallback implementation
    def levenshtein_distance_calc(seq1, seq2):
        # Slower Python implementation of Levenshtein distance
        len1, len2 = len(seq1), len(seq2)
        dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]

        for i in range(len1 + 1):
            dp[i][0] = i
        for j in range(len2 + 1):
            dp[0][j] = j

        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                if seq1[i - 1] == seq2[j - 1]:
                    cost = 0
                else:
                    cost = 1
                dp[i][j] = min(dp[i - 1][j] + 1,  # Deletion
                               dp[i][j - 1] + 1,  # Insertion
                               dp[i - 1][j - 1] + cost)  # Substitution

        return dp[len1][len2]
#####

################### We are currently fixed using Table 11
codon_table = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def translate_frame(sequence):
    translate = ''.join([codon_table.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

def translate_dna_to_aa(dna_fasta, aa_fasta):
    def translate_dna_sequence(dna_seq):
        aa_seq = ""
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3]
            aa_seq += codon_table.get(codon, 'X')  # 'X' for unknown codons
        return aa_seq

    with open(dna_fasta, 'r') as infile, open(aa_fasta, 'w') as outfile:
        dna_seq = ""
        header = ""
        for line in infile:
            if line.startswith('>'):
                if dna_seq:
                    aa_seq = translate_dna_sequence(dna_seq)
                    wrapped_aa_seq = wrap_sequence(aa_seq, 60)
                    outfile.write(f"{header}\n{wrapped_aa_seq}\n")
                header = line.strip()
                dna_seq = ""
            else:
                dna_seq += line.strip()
        if dna_seq:
            aa_seq = translate_dna_sequence(dna_seq)
            wrapped_aa_seq = wrap_sequence(aa_seq, 60)
            outfile.write(f"{header}\n{wrapped_aa_seq}\n")


def detect_sequence_type(fasta_file):
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            if any(base in line for base in 'EFILPQZ'):
                return False  # Contains amino acids
    return True  # Contains DNA


def is_tool_installed(tool_name):
    """Check if a tool is installed and available in PATH."""
    # Check if the tool is in the system PATH
    if shutil.which(tool_name) is None:
        return False

    # Try running the tool to ensure it's executable
    try:
        subprocess.run([tool_name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except subprocess.CalledProcessError:
        return True  # The tool is installed and ran, even if it returns an error code
    except FileNotFoundError:
        return False  # This shouldn't happen due to the earlier check

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))


def fix_path(path):
    fixed_path = os.path.normpath(path)
    fixed_path = os.path.realpath(fixed_path)
    return fixed_path


def extract_identity(clustered_info):
    # Use regex to capture percentage, including optional '-' or '+' before it
    match = re.search(r'at [+-/]*(\d+\.\d+)%', clustered_info)

    if match:
        percent_identity = float(match.group(1))  # Extract the percentage value
        return percent_identity
    else:
        raise ValueError("Percent identity not found in the string.")


def wrap_sequence(sequence, width=60):
    wrapped_sequence = []
    for i in range(0, len(sequence), width):
        wrapped_sequence.append(sequence[i:i + width])
    return "\n".join(wrapped_sequence)


def read_genomes_from_fasta(fasta_file):
    genomes = set()
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                genome = line.split('|')[1]
                genomes.add(genome)
    return list(genomes)

def read_fasta(fasta_file):
    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            if line.startswith('>'):
                current_sequence = line[1:]  # Remove '>' character
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line
    return sequences


def reorder_dict_by_keys(original_dict, sorted_keys):
    return {k: original_dict[k] for k in sorted_keys}
def custom_sort_key(k, dict1, dict2):
    return (len(dict1[k]), len(dict2[k]))

def sort_keys_by_values(dict1, dict2):
    sorted_keys = sorted(dict1.keys(), key=lambda k: custom_sort_key(k, dict1, dict2), reverse=True)
    return sorted_keys

def select_longest_gene(sequences, subgrouped):
    """Select the longest sequence for each genome."""
    longest_sequences = {}
    for seq_id, sequence in sequences.items():
        if subgrouped == False:
            genome = seq_id.split('|')[0]  # Assuming genome name can be derived from the sequence ID
        elif subgrouped == True:
            genome = seq_id.split('|')[1]
        if genome not in longest_sequences or len(sequence) > len(longest_sequences[genome][1]):
            longest_sequences[genome] = (seq_id, sequence)
    return longest_sequences


def run_mafft_on_sequences(options, sequences, output_file):
    #print("Conducting MAFFT alignment.")
    """Run mafft on the given sequences and write to output file."""
    # Create a temporary input file for mafft
    with NamedTemporaryFile('w', delete=False) as temp_input_file:
        for header, sequence in sequences.items():
            temp_input_file.write(f">{header}\n{sequence}\n")
        temp_input_file_path = temp_input_file.name

    # Run mafft
    try:
        with open(output_file, 'w') as output_f:
            if options.verbose == True:
                subprocess.run(
                    ['mafft', '--auto', '--thread', str(options.threads), temp_input_file_path],
                    stdout=output_f,
                    stderr=sys.stderr,
                    check=True
                )

            else:
                subprocess.run(
                    ['mafft', '--auto', '--thread', str(options.threads), temp_input_file_path],
                    stdout=output_f,
                    stderr=subprocess.DEVNULL,  # Suppress stderr
                    check=True
                )
    finally:
        os.remove(temp_input_file_path)  # Clean up the temporary file




def read_separate_files(input_dir, name_split_gff, name_split_fasta, gene_ident, combined_out, translate, run_as_combiner):
    if run_as_combiner == True:
        combined_out_file_aa = None
    else:
        combined_out_file_aa = combined_out.replace('_dna.fasta','_aa.fasta')

    with open(combined_out, 'w') as combined_out_file, (open(combined_out_file_aa, 'w') if combined_out_file_aa else open(os.devnull, 'w')) as combined_out_file_aa:
        paired_files_found = None
    #with open(combined_out, 'w') as combined_out_file, open(combined_out.replace('_dna.fasta','_aa.fasta'), 'w') as combined_out_file_aa:
        gff_files = glob.glob(os.path.join(input_dir, '*' + name_split_gff))
        if not gff_files:
            sys.exit("Error: No GFF files found.")
        for gff_file in gff_files:
            genome_name = os.path.basename(gff_file).split(name_split_gff)[0]
            if name_split_fasta == None:
                possible_extensions = ['.fa', '.fasta', '.fna']
                corresponding_fasta_file = None
                for ext in possible_extensions:
                    temp_file = os.path.splitext(gff_file)[0] + ext
                    if os.path.exists(temp_file):
                        corresponding_fasta_file = temp_file
                        break
                if corresponding_fasta_file is None:
                    print("Corresponding FASTA file for GFF file '" + gff_file + "' not found. Skipping. - Try using the -name_split_fasta option.")
                    continue
            else:
                corresponding_fasta_file = os.path.join(input_dir, genome_name + name_split_fasta)
                if not os.path.exists(corresponding_fasta_file):
                    print("Corresponding FASTA file for GFF file '" + gff_file + "' not found. Skipping. - Try using the -name_split_fasta option.")
                    continue

            gff_features = []
            paired_files_found = True
            with open(gff_file, 'r') as file:
                seen_seq_ids = collections.defaultdict(int)
                lines = file.readlines()
                for line in lines:
                    line_data = line.split('\t')
                    if len(line_data) == 9:
                        if any(gene_type in line_data[2] for gene_type in gene_ident):
                            seq_id = line_data[8].split('ID=')[1].split(';')[0]
                            contig = line_data[0]
                            feature = line_data[2]
                            strand = line_data[6]
                            start, end = int(line_data[3]), int(line_data[4])
                            if seq_id in seen_seq_ids:
                                seq_id += '_' + str(seen_seq_ids[seq_id])
                                seen_seq_ids[seq_id] + 1
                            else:
                                seen_seq_ids[seq_id] = 1
                            gff_features.append((contig, start, end, strand, feature,  seq_id))
            fasta_dict = collections.defaultdict(str)
            with open(corresponding_fasta_file, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith('>'):
                        current_contig = line[1:].split()[0]
                        fasta_dict[current_contig] = ['', '']
                    else:
                        fasta_dict[current_contig][0] += line.strip()

            for contig, fasta in fasta_dict.items():
                reverse_sequence = reverse_complement(fasta[0])
                fasta_dict[contig][1] = reverse_sequence

            if fasta_dict and gff_features:
                for contig, start, end, strand, feature, seq_id in gff_features:
                    if contig in fasta_dict:
                        if strand == '+':
                            full_sequence = fasta_dict[contig][0]
                            seq = full_sequence[start - 1:end]
                        elif strand == '-':
                            corrected_start = max(len(fasta_dict[contig][0]) - int(end), 1)
                            corrected_stop = max(len(fasta_dict[contig][0]) - int(start - 1), 1)
                            full_sequence = fasta_dict[contig][1]
                            seq = full_sequence[corrected_start:corrected_stop]

                        if run_as_combiner == True:
                            if translate == True:
                                seq_aa = translate_frame(seq)
                                wrapped_sequence_aa = '\n'.join([seq_aa[i:i + 60] for i in range(0, len(seq_aa), 60)])
                                combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence_aa}\n")
                            else:
                                wrapped_sequence = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                                combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")
                        else:
                            if translate == True:
                                seq_aa = translate_frame(seq)
                                wrapped_sequence_aa = '\n'.join([seq_aa[i:i + 60] for i in range(0, len(seq_aa), 60)])
                                combined_out_file_aa.write(f">{genome_name}|{seq_id}\n{wrapped_sequence_aa}\n")
                            wrapped_sequence = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                            combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")

    if not paired_files_found:
        sys.exit("Could not find matching GFF/FASTA files - Please check input directory and -name_split_gff and -name_split_fasta parameters.")
    if translate == False or translate == None:
        #Clean up unused file
        try: # Catches is combined_out_file_aa is None
            if combined_out_file.name != combined_out_file_aa.name:
                os.remove(combined_out_file_aa.name)
        except AttributeError:
            pass


def read_combined_files(input_dir, name_split, gene_ident, combined_out, translate, run_as_combiner):
    if run_as_combiner == True:
        combined_out_file_aa = None
    else:
        combined_out_file_aa = combined_out.replace('_dna.fasta','_aa.fasta')
    #with open(combined_out, 'w') as combined_out_file, open(combined_out_file_aa, 'w') if combined_out_file_aa else open(os.devnull, 'w'):
    with open(combined_out, 'w') as combined_out_file, (open(combined_out_file_aa, 'w') if combined_out_file_aa else open(os.devnull, 'w')) as combined_out_file_aa:
        gff_files = glob.glob(os.path.join(input_dir, '*' + name_split))
        if not gff_files:
            sys.exit("Error: No GFF files found - check input directory and -name_split_gff parameter.")
        for gff_file in gff_files:
            genome_name = os.path.basename(gff_file).split(name_split)[0]
            fasta_dict = collections.defaultdict(str)
            gff_features = []
            with open(gff_file, 'r') as file:
                seen_seq_ids = collections.defaultdict(int)
                lines = file.readlines()
                fasta_section = False
                for line in lines:
                    if line.startswith('##FASTA'):
                        fasta_section = True
                        continue
                    if fasta_section:
                        if line.startswith('>'):
                            current_contig = line[1:].split()[0]
                            fasta_dict[current_contig] = ['','']
                        else:
                            fasta_dict[current_contig][0] +=line.strip()
                    else:
                        line_data = line.split('\t')
                        if len(line_data) == 9:
                            if any(gene_type in line_data[2] for gene_type in gene_ident):
                                contig = line_data[0]
                                feature = line_data[2]
                                strand = line_data[6]
                                start, end = int(line_data[3]), int(line_data[4])
                                seq_id = line_data[8].split('ID=')[1].split(';')[0]
                                if seq_id in seen_seq_ids:
                                    seq_id += '_' + str(seen_seq_ids[seq_id])
                                    seen_seq_ids[seq_id] + 1
                                else:
                                    seen_seq_ids[seq_id] = 1
                                gff_features.append((contig, start, end, strand, feature,  seq_id))

                for contig, fasta in fasta_dict.items():
                    reverse_sequence = reverse_complement(fasta[0])
                    fasta_dict[contig][1] = reverse_sequence

                if fasta_dict and gff_features:
                    for contig, start, end, strand, feature, seq_id in gff_features:
                        if contig in fasta_dict:
                            if strand == '+':
                                full_sequence = fasta_dict[contig][0]
                                seq = full_sequence[start - 1:end]
                            elif strand == '-':
                                corrected_start = max(len(fasta_dict[contig][0]) - int(end), 1)
                                corrected_stop = max(len(fasta_dict[contig][0]) - int(start - 1), 1)
                                full_sequence = fasta_dict[contig][1]
                                seq = full_sequence[corrected_start:corrected_stop]

                            if run_as_combiner == True:
                                if translate == True:
                                    seq_aa = translate_frame(seq)
                                    wrapped_sequence_aa = '\n'.join([seq_aa[i:i + 60] for i in range(0, len(seq_aa), 60)])
                                    combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence_aa}\n")
                                else:
                                    wrapped_sequence = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                                    combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")
                            else:
                                if translate == True:
                                    seq_aa = translate_frame(seq)
                                    wrapped_sequence_aa = '\n'.join([seq_aa[i:i + 60] for i in range(0, len(seq_aa), 60)])
                                    combined_out_file_aa.write(f">{genome_name}|{seq_id}\n{wrapped_sequence_aa}\n")
                                wrapped_sequence = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                                combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")

    if translate == False or translate == None:
        #Clean up unused file
        try: # Catches is combined_out_file_aa is None
            if combined_out_file.name != combined_out_file_aa.name:
                os.remove(combined_out_file_aa.name)
        except AttributeError:
            pass



def read_fasta_files(input_dir, name_split_fasta, combined_out, translate, run_as_combiner):
    if run_as_combiner == True:
        combined_out_file_aa = None
    else:
        combined_out_file_aa = combined_out.replace('_dna.fasta','_aa.fasta')
    with open(combined_out, 'w') as combined_out_file, (open(combined_out_file_aa, 'w') if combined_out_file_aa else open(os.devnull, 'w')) as combined_out_file_aa:
        fasta_files = glob.glob(os.path.join(input_dir, '*' + name_split_fasta))
        if not fasta_files:
            sys.exit("Error: No GFF files found.")
        for fasta_file in fasta_files:
            genome_name = os.path.basename(fasta_file).split(name_split_fasta)[0]
            fasta_dict = collections.defaultdict(str)
            with open(fasta_file, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith('>'):
                        current_seq = line[1:].split()[0]
                        fasta_dict[current_seq] = ''
                    else:
                        fasta_dict[current_seq] +=line.strip()
                for seq_id, seq in fasta_dict.items():
                    if run_as_combiner == True:
                        if translate == True:
                            seq_aa = translate_frame(seq)
                            wrapped_sequence_aa = '\n'.join([seq_aa[i:i + 60] for i in range(0, len(seq_aa), 60)])
                            combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence_aa}\n")
                        else:
                            wrapped_sequence = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                            combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")
                    else:
                        if translate == True:
                            seq_aa = translate_frame(seq)
                            wrapped_sequence_aa = '\n'.join([seq_aa[i:i + 60] for i in range(0, len(seq_aa), 60)])
                            combined_out_file_aa.write(f">{genome_name}|{seq_id}\n{wrapped_sequence_aa}\n")
                        wrapped_sequence = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                        combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")

    if translate == False or translate == None:
        #Clean up unused file
        try: # Catches is combined_out_file_aa is None
            if combined_out_file.name != combined_out_file_aa.name:
                os.remove(combined_out_file_aa.name)
        except AttributeError:
            pass


def write_groups_func(options, output_dir, key_order, cores, sequences,
                 pangenome_clusters_First_sequences_sorted, combined_pangenome_clusters_Second_sequences):
    """
    Writes individual FASTA files and a combined FASTA file for all sequences.

    Parameters:
    - options: Command-line options.
    - output_dir: Directory where output FASTA files will be saved.
    - key_order: The order in which to process keys.
    - cores: Dictionary of core genes.
    - sequences: Dictionary mapping headers to sequences.
    - pangenome_clusters_First_sequences_sorted: Dictionary of first sequence clusters.
    - combined_pangenome_clusters_Second_sequences: Dictionary of second sequence clusters.
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for group in options.write_groups.split(','):

        combined_fasta_filename = os.path.join(output_dir, "combined_group_sequences_" + group + "_dna.fasta")

        # Open combined FASTA file for writing all sequences
        with open(combined_fasta_filename, 'w') as combined_fasta, open(combined_fasta_filename.replace('_dna.fasta','_aa.fasta'), 'w') as combined_fasta_aa:
            for key_prefix in key_order:
                for key, values in cores.items():
                    if any(part in group for part in key.split('_')):
                        if key.startswith(key_prefix):
                            for value in values:
                                output_filename = f"{key}_{value}_dna.fasta"
                                if 'First' in key_prefix:
                                    sequences_to_write = pangenome_clusters_First_sequences_sorted[value]
                                else:
                                    sequences_to_write = combined_pangenome_clusters_Second_sequences[value]

                                # Write individual FASTA file
                                with open(os.path.join(output_dir,output_filename), 'w') as outfile, open(os.path.join(output_dir, output_filename.replace('_dna.fasta','_aa.fasta')), 'w') as outfile_aa:
                                    for header in sequences_to_write:
                                        if header in sequences:
                                            sequence = sequences[header]
                                            wrapped_sequence = wrap_sequence(sequence)
                                            # Handle Amino Acid Sequences (AA)
                                            if options.sequence_type == 'AA':
                                                seq_aa = translate_frame(sequence)
                                                wrapped_sequence_aa = wrap_sequence(seq_aa)
                                                # Write individual group file for AA, if option is enabled
                                                if options.write_individual_groups:
                                                    outfile_aa.write(f">{header}\n")
                                                    outfile_aa.write(f"{wrapped_sequence_aa}\n")
                                                else:
                                                    os.remove(outfile_aa.name)  # Delete individual file if option is disabled
                                                # Always write to the combined AA file
                                                combined_fasta_aa.write(f">Group_{value}|{header}\n")
                                                combined_fasta_aa.write(f"{wrapped_sequence_aa}\n")
                                            # Handle Nucleotide Sequences
                                            else:
                                                # If the option is disabled, delete individual AA file (if created)
                                                try:
                                                    os.remove(outfile_aa.name)  # Ensure outfile_aa is removed when sequence_type isn't 'AA'
                                                except FileNotFoundError:
                                                    pass
                                            # Write individual group file for nucleotide sequence, if option is enabled
                                            if options.write_individual_groups:
                                                outfile.write(f">{header}\n")
                                                outfile.write(f"{wrapped_sequence}\n")
                                            else:
                                                os.remove(outfile.name)  # Delete individual file if option is disabled
                                            # Always write to the combined nucleotide file
                                            combined_fasta.write(f">Group_{value}|{header}\n")
                                            combined_fasta.write(f"{wrapped_sequence}\n")

                                        else:
                                            if options.verbose == True:
                                                print(f"Sequence {header} not found in original_fasta file.")
        if options.sequence_type != 'AA':
            #Clean up unused file
            os.remove(combined_fasta_aa.name)
    print(f"Combined FASTA file saved to: {combined_fasta_filename}")


# def process_gene_groups(options, group_directory, sub_group_directory, paralog_groups, output_file):
#     """Process each gene family file to select the longest sequence per genome and concatenate aligned sequences."""
#     concatenated_sequences = {}
#     output_file = group_directory.replace('Gene_Groups_Output',output_file)
#
#     # Iterate over each gene family file
#     for gene_file in os.listdir(group_directory):
#         if gene_file.endswith('.fasta') and not gene_file.endswith('combined_group_sequences.fasta') :
#             gene_path = os.path.join(group_directory, gene_file)
#
#             # Read sequences from the gene family file
#             sequences = read_fasta(gene_path)
#
#             # Select the longest sequence for each genome
#             longest_sequences = select_longest_gene(sequences)
#
#             # Run mafft on the longest sequences
#             aligned_file = f"{group_directory}/{gene_file}_aligned.fasta.tmp"
#             run_mafft_on_sequences(options, {seq_id: seq for seq_id, seq in longest_sequences.values()}, aligned_file)
#
#             # Read aligned sequences and concatenate them
#             aligned_sequences = read_fasta(aligned_file)
#             for genome, aligned_seq in aligned_sequences.items():
#                 genome_name = genome.split('|')[0]
#                 if genome_name not in concatenated_sequences:
#                     concatenated_sequences[genome_name] = ""
#                 concatenated_sequences[genome_name] += aligned_seq
#
#             # Clean up aligned file
#             os.remove(aligned_file)
#
#     # Write the concatenated sequences to the output file
#     with open(output_file, 'w') as out:
#         for genome, sequence in concatenated_sequences.items():
#             out.write(f">{genome}\n")
#             wrapped_sequence = wrap_sequence(sequence, 60)
#             out.write(f"{wrapped_sequence}\n")

def perform_alignment(gene_path,group_directory, gene_file, options, concatenated_sequences, subgrouped):
    # Read sequences from the gene family file
    sequences = read_fasta(gene_path)
    if len(sequences) == 1: # We can't align a single sequence
        return concatenated_sequences
    # Select the longest sequence for each genome
    longest_sequences = select_longest_gene(sequences, subgrouped)

    # Run mafft on the longest sequences
    aligned_file = f"{group_directory}/{gene_file}_aligned.fasta.tmp"
    run_mafft_on_sequences(options, {seq_id: seq for seq_id, seq in longest_sequences.values()}, aligned_file)

    # Read aligned sequences and concatenate them
    aligned_sequences = read_fasta(aligned_file)
    # Find the length of the longest sequence in aligned_sequences
    max_length = max(len(seq) for seq in aligned_sequences.values())

    for genome, sequence in concatenated_sequences.items():
        if any(genome in key for key in aligned_sequences.keys()):
            genome_name_in_aligned = next(key for key in aligned_sequences.keys() if genome in key)#.split('|')[split_by]
            concatenated_sequences[genome] += aligned_sequences[genome_name_in_aligned]
        else:
            concatenated_sequences[genome] += "-" * max_length

    # Clean up aligned file
    os.remove(aligned_file)

    return concatenated_sequences

def process_gene_groups(options, group_directory, sub_group_directory, paralog_groups, genome_list, output_file):
    """Process each gene family file to select the longest sequence per genome and concatenate aligned sequences."""
    concatenated_sequences = {genome: "" for genome in genome_list}
    output_file = group_directory.replace('Gene_Groups_Output', output_file)
    if paralog_groups != None:
        threshold_size = math.floor(int(options.align_core) * int(options.genome_num) / 100)

    if options.align_aa == True:
        affix = '_aa.fasta'
    else:
        affix = '_dna.fasta'

    if options.align_core == True:
        # Iterate over each gene family file
        for gene_file in os.listdir(group_directory):
            if gene_file.endswith(affix) and not gene_file.startswith('combined_group_sequences'):
                current_group = int(gene_file.split('_')[3].split('.')[0])
                gene_path = os.path.join(group_directory, gene_file)
                # Could add more catches here to work with First and Secondary groups - This ensures only core '99/100' are aligned
                if 'First_core_99' in gene_file or 'First_core_100' in gene_file:
                    # Check for matching group in paralog_groups
                    if sub_group_directory and paralog_groups and '>Group_'+str(current_group) in paralog_groups:
                        for subgroup, size in enumerate(paralog_groups['>Group_' + str(current_group)]['sizes']):
                            if size >= threshold_size:
                                gene_path = os.path.join(sub_group_directory,f"Group_{current_group}_subgroup_{subgroup}{affix}")
                                concatenated_sequences = perform_alignment(gene_path, group_directory, gene_file, options, concatenated_sequences, True)
                    else:
                        concatenated_sequences = perform_alignment(gene_path, group_directory, gene_file, options, concatenated_sequences, False)

    # Write the concatenated sequences to the output file
    with open(output_file, 'w') as out:
        for genome, sequence in concatenated_sequences.items():
            out.write(f">{genome}\n")
            wrapped_sequence = wrap_sequence(sequence, 60)
            out.write(f"{wrapped_sequence}\n")

