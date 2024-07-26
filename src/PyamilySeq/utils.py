import subprocess
import shutil
import os
import glob
import collections
from tempfile import NamedTemporaryFile
import sys


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
    print("Conducting MAFFT alignment.")
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
                    ['mafft', '--auto', temp_input_file_path],
                    stdout=output_f,
                    stderr=sys.stderr,
                    check=True
                )
            else:
                subprocess.run(
                    ['mafft', '--auto', temp_input_file_path],
                    stdout=output_f,
                    stderr=subprocess.DEVNULL,  # Suppress stderr
                    check=True
                )
    finally:
        os.remove(temp_input_file_path)  # Clean up the temporary file




def read_separate_files(input_dir, name_split, combined_out):
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
                        if line_data[2] == 'CDS':
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

                        wrapped_sequence = '\n'.join([cds_sequence[i:i + 60] for i in range(0, len(cds_sequence), 60)])
                        combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")


def read_combined_files(input_dir, name_split, combined_out):
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
                            if line_data[2] == 'CDS':
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

                            wrapped_sequence = '\n'.join([cds_sequence[i:i + 60] for i in range(0, len(cds_sequence), 60)])
                            combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")