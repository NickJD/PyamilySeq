import argparse
import collections
import os
import glob
import subprocess
from PyamilySeq_Species import *


try:
    from .PyamilySeq_Species import cluster
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from PyamilySeq_Species import cluster
    from Constants import *

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))


def read_separate_files(input_dir, name_split, combined_out):
    with open(combined_out, 'w') as combined_out_file:
        for fasta_file in glob.glob(os.path.join(input_dir, '*' + name_split)):
            genome_name = os.path.basename(fasta_file).split(name_split)[0]
            corresponding_gff_file = fasta_file.replace('.fasta', '.gff')
            if not os.path.exists(corresponding_gff_file):
                continue
            cds_sequences = extract_cds_from_gff(fasta_file, corresponding_gff_file)
            for gene_name, seq in cds_sequences:
                header = f">{genome_name}_{gene_name}\n"
                combined_out_file.write(header)
                combined_out_file.write(seq + '\n')

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
                            fasta_dict[current_contig] = []
                        else:
                            fasta_dict[current_contig].append(line.strip())
                    else:
                        line_data = line.split('\t')
                        if len(line_data) == 9:
                            if line_data[2] == 'CDS':
                                contig = line_data[0]
                                feature = line_data[2]
                                start, end = int(line_data[3]), int(line_data[4])
                                seq_id = line_data[8].split('ID=')[1].split(';')[0]
                                gff_features.append((contig, start, end, seq_id))

                if fasta_dict and gff_features:
                    for contig, start, end, seq_id in gff_features:
                        if contig in fasta_dict:
                            full_sequence = ''.join(fasta_dict[contig])
                            cds_sequence = full_sequence[start - 1:end]
                            wrapped_sequence = '\n'.join([cds_sequence[i:i + 60] for i in range(0, len(cds_sequence), 60)])
                            combined_out_file.write(f">{genome_name}|{seq_id}\n{wrapped_sequence}\n")


def run_cd_hit(input_file, clustering_output, options):
    cdhit_command = [
        'cd-hit-est',
        '-i', input_file,
        '-o', clustering_output,
        '-c', str(options.pident),
        '-s', str(options.len_diff),
        '-T', "20",
        '-d', "0",
        '-sc', "1",
        '-sf', "1"
    ]
    subprocess.run(cdhit_command)








def main():
    parser = argparse.ArgumentParser(
        description='PyamilySeq ' + PyamilySeq_Version + ': PyamilySeq Run Parameters.')
    required = parser.add_argument_group('Required Arguments')
    required.add_argument("-id", action="store", dest="input_dir",
                          help="Directory containing GFF/FASTA files.",
                          required=True)
    required.add_argument("-od", action="store", dest="output_dir",
                          help="Directory for all output files.",
                          required=True)
    required.add_argument("-it", action="store", dest="input_type", choices=['separate', 'combined'],
                        help="Type of input files: 'separate' for separate FASTA and GFF files,"
                             " 'combined' for GFF files with embedded FASTA sequences.",
                          required=True)
    required.add_argument("-ns", action="store", dest="name_split",
                          help="Character used to split the filename and extract the genome name.",
                          required=True)
    required.add_argument("-pid", action="store", dest="pident", type=float,
                          help="Pident threshold for CD-HIT clustering.",
                          required=True)
    required.add_argument("-ld", action="store", dest="len_diff", type=float,
                          help="Length difference (-s) threshold for CD-HIT clustering.",
                          required=True)
    required.add_argument("-co", action="store", dest="clustering_out",
                          help="Output file for initial clustering.",
                          required=True)
    required.add_argument("-ct", action="store", dest="clustering_type", choices=['CD-HIT', 'BLAST', 'DIAMOND', "MMseqs2"],
                        help="Clustering format for PyamilySeq.",
                          required=True)

    output_args = parser.add_argument_group('Output Parameters')
    output_args.add_argument('-w', action="store", dest='write_families', default=None,
                          help='Default - No output: Output sequences of identified families (provide levels at which to output "-w 99,95"'
                               ' - Must provide FASTA file with -fasta',
                          required=False)
    output_args.add_argument('-con', action="store", dest='con_core', default=None,
                          help='Default - No output: Output aligned and concatinated sequences of identified families - used for MSA (provide levels at which to output "-w 99,95"'
                               ' - Must provide FASTA file with -fasta',
                          required=False)
    output_args.add_argument('-fasta', action='store', dest='fasta',
                          help='FASTA file to use in conjunction with "-w" or "-con"',
                          required=False)

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rc', action='store', dest='reclustered', help='Clustering output file from secondary round of clustering',
                        required=False)
    optional.add_argument('-st', action='store', dest='sequence_tag', help='Default - "StORF": Unique identifier to be used to distinguish the second of two rounds of clustered sequences',
                        required=False)
    optional.add_argument('-groups', action="store", dest='core_groups', default="99,95,15",
                        help='Default - (\'99,95,15\'): Gene family groups to use')
    optional.add_argument('-gpa', action='store', dest='gene_presence_absence_out', help='Default - False: If selected, a Roary formatted gene_presence_absence.csv will be created - Required for Coinfinder and other downstream tools',
                        required=False)

    parser.add_argument("pyamilyseq_args", nargs=argparse.REMAINDER, help="Additional arguments for PyamilySeq.")
    options = parser.parse_args()



    output_path = os.path.abspath(options.output_dir)
    combined_out_file = os.path.join(output_path,"end_to_end_combined_sequences.fasta")
    clustering_output = os.path.join(output_path,'clustering_'+options.clustering_type)



    # Step 1: Read and rename sequences from files based on input type
    if options.input_type == 'separate':
        read_separate_files(options.input_dir, options.name_split, combined_out_file)
    else:
        read_combined_files(options.input_dir, options.name_split, combined_out_file)

    # Step 2: Run CD-HIT on the renamed sequences
    run_cd_hit(combined_out_file, clustering_output, options)


    class clustering_options:
        def __init__(self):
            self.format = 'CD-HIT'
            self.reclustered = options.reclustered
            self.sequence_tag = 'StORF'
            self.core_groups = '99,95,15,0'
            self.clusters = clustering_output+'.clstr'
            self.gene_presence_absence_out = options.gene_presence_absence_out
            self.write_families = options.write_families
            self.con_core = options.con_core

    clustering_options = clustering_options()

    # Step 3: Run PyamilySeq with the CD-HIT output
    cluster(clustering_options)
    #run_pyamilyseq(options.clustering_out, options.clustering_type, combined_out_file, options.pyamilyseq_args)


if __name__ == "__main__":
    main()
