import subprocess
import shutil
import os
import glob
import collections
from tempfile import NamedTemporaryFile
import sys
from line_profiler_pycharm import profile
import re


################### We are currently fixed using Table 11
gencode = {
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
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

@profile
def calculate_similarity(seq1, seq2):
    len1, len2 = len(seq1), len(seq2)

    # If lengths are the same, directly compare without alignment
    if len1 == len2:
        matches = sum(c1 == c2 for c1, c2 in zip(seq1, seq2))
        return (matches / len1) * 100  # Return similarity based on the length

    # For different lengths, proceed with global alignment
    # Initialize the scoring matrix
    score_matrix = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    # Fill the first row and first column with gap penalties
    for i in range(len1 + 1):
        score_matrix[i][0] = -i  # Gap penalty for seq1
    for j in range(len2 + 1):
        score_matrix[0][j] = -j  # Gap penalty for seq2

    # Fill the score matrix
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match = score_matrix[i - 1][j - 1] + (1 if seq1[i - 1] == seq2[j - 1] else -1)
            delete = score_matrix[i - 1][j] - 1  # Gap in seq2
            insert = score_matrix[i][j - 1] - 1  # Gap in seq1
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback to find the alignment (if needed for detailed output)
    aligned_seq1, aligned_seq2 = "", ""
    i, j = len1, len2

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        if i > 0 and j > 0 and current_score == score_matrix[i - 1][j - 1] + (1 if seq1[i - 1] == seq2[j - 1] else -1):
            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif i > 0 and current_score == score_matrix[i - 1][j] - 1:
            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += "-"
            i -= 1
        else:
            aligned_seq1 += "-"
            aligned_seq2 += seq2[j - 1]
            j -= 1

    # Reverse the aligned sequences if needed
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    # Calculate matches from aligned sequences
    matches = sum(c1 == c2 for c1, c2 in zip(aligned_seq1, aligned_seq2))

    # Calculate the similarity percentage based on the maximum length
    max_length = max(len(seq1), len(seq2))
    return (matches / max_length) * 100



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
    # Use regular expressions to capture the percentage value at the end of the line
    match = re.search(r'at ([-+]*)(\d+\.\d+)%', clustered_info)

    if match:
        percent_identity = float(match.group(2))  # Extract the percentage value
        return percent_identity
    else:
        raise ValueError("Percent identity not found in the string.")

def wrap_sequence(sequence, width=60):
    wrapped_sequence = []
    for i in range(0, len(sequence), width):
        wrapped_sequence.append(sequence[i:i + width])
    return "\n".join(wrapped_sequence)


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

def select_longest_gene(sequences):
    """Select the longest sequence for each genome."""
    longest_sequences = {}
    for seq_id, sequence in sequences.items():
        genome = seq_id.split('|')[0]  # Assuming genome name can be derived from the sequence ID
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




def read_separate_files(input_dir, name_split, gene_ident, combined_out, translate):
    with open(combined_out, 'w') as combined_out_file:
        for gff_file in glob.glob(os.path.join(input_dir, '*' + name_split)):
            genome_name = os.path.basename(gff_file).split(name_split)[0]
            corresponding_fasta_file = os.path.splitext(gff_file)[0] + '.fa'
            if not os.path.exists(corresponding_fasta_file):
                continue

            gff_features = []
            with open(gff_file, 'r') as file:
                seen_seq_ids = collections.defaultdict(int)
                lines = file.readlines()
                for line in lines:
                    line_data = line.split('\t')
                    if len(line_data) == 9:
                        if any(gene_type in line_data[2] for gene_type in gene_ident):
                            contig = line_data[0]
                            feature = line_data[2]
                            strand = line_data[6]
                            start, end = int(line_data[3]), int(line_data[4])
                            if seq_id in seen_seq_ids:
                                seq_id += '_' + str(seen_seq_ids[seq_id])
                                seen_seq_ids[seq_id] + 1
                            else:
                                seen_seq_ids[seq_id] = 1
                            seq_id = line_data[8].split('ID=')[1].split(';')[0]
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
                            cds_sequence = full_sequence[start - 1:end]
                        elif strand == '-':
                            corrected_start = max(len(fasta_dict[contig][0]) - int(end), 1)
                            corrected_stop = max(len(fasta_dict[contig][0]) - int(start - 1), 1)
                            full_sequence = fasta_dict[contig][1]
                            cds_sequence = full_sequence[corrected_start:corrected_stop]
                        if translate == True:
                            cds_sequence = translate_frame(cds_sequence)
                        wrapped_sequence = '\n'.join([cds_sequence[i:i + 60] for i in range(0, len(cds_sequence), 60)])
                        combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")


def read_combined_files(input_dir, name_split, gene_ident, combined_out, translate):
    with open(combined_out, 'w') as combined_out_file:
        for gff_file in glob.glob(os.path.join(input_dir, '*' + name_split)):
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
                    fasta_dict[contig][1]=reverse_sequence

                if fasta_dict and gff_features:
                    for contig, start, end, strand, feature, seq_id in gff_features:
                        if contig in fasta_dict:
                            if strand == '+':
                                full_sequence = fasta_dict[contig][0]
                                cds_sequence = full_sequence[start - 1:end]
                            elif strand == '-':
                                corrected_start = max(len(fasta_dict[contig][0]) - int(end), 1)
                                corrected_stop = max(len(fasta_dict[contig][0]) - int(start - 1), 1)
                                full_sequence = fasta_dict[contig][1]
                                cds_sequence = full_sequence[corrected_start:corrected_stop]

                            if translate == True:
                                cds_sequence = translate_frame(cds_sequence)
                            wrapped_sequence = '\n'.join([cds_sequence[i:i + 60] for i in range(0, len(cds_sequence), 60)])
                            combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")


def read_fasta_files(input_dir, name_split, combined_out, translate):
    with open(combined_out, 'w') as combined_out_file:
        for fasta_file in glob.glob(os.path.join(input_dir, '*' + name_split)):
            genome_name = os.path.basename(fasta_file).split(name_split)[0]
            fasta_dict = collections.defaultdict(str)
            with open(fasta_file, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith('>'):
                        current_seq = line[1:].split()[0]
                        fasta_dict[current_seq] = ''
                    else:
                        fasta_dict[current_seq] +=line.strip()
                for id, seq in fasta_dict.items():
                    if translate == True:
                        seq = translate_frame(seq)
                    wrapped_sequence = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                    combined_out_file.write(f">{genome_name}|{id}\n{wrapped_sequence}\n")


def write_groups(options, output_dir, key_order, cores, sequences,
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

    combined_fasta_filename = os.path.join(output_dir, "combined_group_sequences.fasta")

    # Open combined FASTA file for writing all sequences
    with open(combined_fasta_filename, 'w') as combined_fasta:
        for key_prefix in key_order:
            for key, values in cores.items():
                if any(part in options.write_groups.split(',') for part in key.split('_')):
                    if key.startswith(key_prefix):
                        for value in values:
                            output_filename = f"{key}_{value}.fasta"
                            if 'First' in key_prefix:
                                sequences_to_write = pangenome_clusters_First_sequences_sorted[value]
                            else:
                                sequences_to_write = combined_pangenome_clusters_Second_sequences[value]

                            # Write individual FASTA file
                            with open(os.path.join(output_dir, output_filename), 'w') as outfile:
                                for header in sequences_to_write:
                                    if header in sequences:
                                        sequence = sequences[header]
                                        outfile.write(f">{header}\n")
                                        wrapped_sequence = wrap_sequence(sequence)
                                        outfile.write(f"{wrapped_sequence}\n")

                                        # Also write to the combined FASTA file
                                        combined_fasta.write(f">Group_{value}|{header}\n")
                                        combined_fasta.write(f"{wrapped_sequence}\n")
                                    else:
                                        if options.verbose:
                                            print(f"Sequence {header} not found in original_fasta file.")

    print(f"Combined FASTA file saved to: {combined_fasta_filename}")


def process_gene_families(options, directory, output_file):
    """Process each gene family file to select the longest sequence per genome and concatenate aligned sequences."""
    concatenated_sequences = {}
    output_file = directory.replace('Gene_Families_Output',output_file)

    # Iterate over each gene family file
    for gene_file in os.listdir(directory):
        if gene_file.endswith('.fasta') and not gene_file.endswith('combined_group_sequences.fasta'):
            gene_path = os.path.join(directory, gene_file)

            # Read sequences from the gene family file
            sequences = read_fasta(gene_path)

            # Select the longest sequence for each genome
            longest_sequences = select_longest_gene(sequences)

            # Run mafft on the longest sequences
            aligned_file = f"{directory}/{gene_file}_aligned.fasta.tmp"
            run_mafft_on_sequences(options, {seq_id: seq for seq_id, seq in longest_sequences.values()}, aligned_file)

            # Read aligned sequences and concatenate them
            aligned_sequences = read_fasta(aligned_file)
            for genome, aligned_seq in aligned_sequences.items():
                genome_name = genome.split('|')[0]
                if 'Group' in genome_name:
                    print(2)
                if genome_name not in concatenated_sequences:
                    concatenated_sequences[genome_name] = ""
                concatenated_sequences[genome_name] += aligned_seq

            # Clean up aligned file
            os.remove(aligned_file)

    # Write the concatenated sequences to the output file
    with open(output_file, 'w') as out:
        for genome, sequence in concatenated_sequences.items():
            out.write(f">{genome}\n")
            wrapped_sequence = wrap_sequence(sequence, 60)
            out.write(f"{wrapped_sequence}\n")