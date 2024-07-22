import subprocess
import shutil
import os
import glob
import collections


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







def read_separate_files(input_dir, name_split, combined_out):
    with open(combined_out, 'w') as combined_out_file:
        for gff_file in glob.glob(os.path.join(input_dir, '*' + name_split)):
            genome_name = os.path.basename(gff_file).split(name_split)[0]
            corresponding_fasta_file = os.path.splitext(gff_file)[0] + '.fa'
            if not os.path.exists(corresponding_fasta_file):
                continue

            gff_features = []
            with open(gff_file, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    line_data = line.split('\t')
                    if len(line_data) == 9:
                        if line_data[2] == 'CDS':
                            contig = line_data[0]
                            feature = line_data[2]
                            strand = line_data[6]
                            start, end = int(line_data[3]), int(line_data[4])
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