# PyamilySeq - !BETA!
**PyamilySeq** (Family Seek) is a Python tool for clustering gene sequences into families based on sequence similarity identified by tools such as CD-HIT, BLAST, DIAMOND or MMseqs2.
This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter publication in NAR (https://doi.org/10.1093/nar/gkad814).

## Features
- **End-to-End**: PyamilySeq can take a directory of GFF+FASTA files, run CD-HIT for clustering and process the results.
- **Clustering**: Supports input from CD-HIT formatted files as well as CSV and TSV edge lists (-outfmt 6 from BLAST/DIAMOND).
- **Reclustering**: Allows for the addition of new sequences post-initial clustering.
- **Output**: Generates a gene 'Roary/Panaroo' formatted presence-absence CSV formatted file for downstream analysis.
  - Align representative sequences using MAFFT.
  - Output concatenated aligned sequences for downstream analysis.
  - Optionally output sequences of identified families.


### Installation
PyamilySeq requires Python 3.6 or higher. Install using pip:

```bash
pip install PyamilySeq
```

### Examples: Below are two examples of running PyamilySeq in its two main modes.
#### 'Full Mode': Will conduct clustering of sequences as part of PyamilySeq run
```bash 
PyamilySeq -run_mode Full -group_mode Species -output_dir ../../test_data/testing -input_type combined -input_dir .../test_data/genomes -name_split _combined.gff3 -pid 0.99 -len_diff 0.99 -clust_tool CD-HIT -gpa True -con True -w 99 -verbose True
```
#### 'Partial Mode': Will take the output of a sequence clustering
```bash
PyamilySeq -run_mode Partial -group_mode Species -output_dir .../test_data/testing -cluster_file .../test_data/CD-HIT/combined_Ensmbl_pep_CD_90_60.clstr -clust_tool CD-HIT -original_fasta .../test_data/combined_Ensmbl_cds.fasta -gpa True -con True -w 99 -verbose True
```

```bash
Calculating Groups
Gene Groups:
first_core_99: 3103
first_core_95: 0
first_core_15: 3217
first_core_0: 4808
Total Number of Gene Groups (Including Singletons): 11128
```


## Usage - Menu
```
usage: PyamilySeq.py [-h] -run_mode {Full,Partial} -group_mode {Species,Genus} -clust_tool {CD-HIT} -output_dir OUTPUT_DIR [-input_type {separate,combined}] [-input_dir INPUT_DIR] [-name_split NAME_SPLIT]
                     [-pid PIDENT] [-len_diff LEN_DIFF] [-mem CLUSTERING_MEMORY] [-t CLUSTERING_THREADS] [-cluster_file CLUSTER_FILE] [-reclustered RECLUSTERED] [-seq_tag SEQUENCE_TAG]
                     [-core_groups CORE_GROUPS] [-genus_groups GENUS_GROUPS] [-w WRITE_FAMILIES] [-con CON_CORE] [-original_fasta ORIGINAL_FASTA] [-gpa GENE_PRESENCE_ABSENCE_OUT] [-verbose {True,False}] [-v]

PyamilySeq v0.6.0: PyamilySeq Run Parameters.

options:
  -h, --help            show this help message and exit

Required Arguments:
  -run_mode {Full,Partial}
                        Run Mode: Should PyamilySeq be run in "Full" or "Partial" mode?
  -group_mode {Species,Genus}
                        Group Mode: Should PyamilySeq be run in "Species" or "Genus" mode?
  -clust_tool {CD-HIT}  Clustering tool to use: CD-HIT, DIAMOND, BLAST or MMseqs2.
  -output_dir OUTPUT_DIR
                        Directory for all output files.

Full-Mode Arguments - Required when "-run_mode Full" is used:
  -input_type {separate,combined}
                        Type of input files: 'separate' for separate FASTA and GFF files, 'combined' for GFF files with embedded FASTA sequences.
  -input_dir INPUT_DIR  Directory containing GFF/FASTA files.
  -name_split NAME_SPLIT
                        substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff').
  -pid PIDENT           Default 0.95: Pident threshold for clustering.
  -len_diff LEN_DIFF    Default 0.80: Minimum length difference between clustered sequences - (-s) threshold for CD-HIT clustering.

Clustering Runtime Arguments - Optional when "-run_mode Full" is used:
  -mem CLUSTERING_MEMORY
                        Default 4000: Memory to be allocated for clustering (in MBs).
  -t CLUSTERING_THREADS
                        Default 4: Threads to be allocated for clustering.

Partial-Mode Arguments - Required when "-run_mode Partial" is used:
  -cluster_file CLUSTER_FILE
                        Clustering output file containing CD-HIT, TSV or CSV Edge List

Grouping Arguments - Use to fine-tune grouping of genes after clustering:
  -reclustered RECLUSTERED
                        Currently only works on Partial Mode: Clustering output file from secondary round of clustering.
  -seq_tag SEQUENCE_TAG
                        Default - "StORF": Unique identifier to be used to distinguish the second of two rounds of clustered sequences
  -core_groups CORE_GROUPS
                        Default - ('99,95,15'): Gene family groups to use for "Species" mode
  -genus_groups GENUS_GROUPS
                        Default - ('1,2,3,4,5,6'): Gene family groups to use for "Genus" mode

Output Parameters:
  -w WRITE_FAMILIES     Default - No output: Output sequences of identified families (provide levels at which to output "-w 99,95" - Must provide FASTA file with -fasta
  -con CON_CORE         Default - No output: Output aligned and concatinated sequences of identified families - used for MSA (provide levels at which to output "-w 99,95" - Must provide FASTA file with -fasta
  -original_fasta ORIGINAL_FASTA
                        FASTA file to use in conjunction with "-w" or "-con" when running in Partial Mode.
  -gpa GENE_PRESENCE_ABSENCE_OUT
                        Default - False: If selected, a Roary formatted gene_presence_absence.csv will be created - Required for Coinfinder and other downstream tools

Misc:
  -verbose {True,False}
                        Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit

```


## Seq-Combiner: This tool is provided to enable the pre-processing of multiple GFF/FASTA files together ready to be clustered by the user
### Example:
```bash
Seq-Combiner -input_dir .../test_data/genomes -name_split _combined.gff3 -output_dir.../test_data -output_name combine_fasta_seqs.fa -input_type combined
```
## Seq-Combiner Menu:
```bash
usage: Seq_Combiner.py [-h] -input_dir INPUT_DIR -input_type {separate,combined} -name_split NAME_SPLIT -output_dir OUTPUT_DIR -output_name OUTPUT_FILE

Seq-Combiner v0.6.0: Seq-Combiner Run Parameters.

options:
  -h, --help            show this help message and exit

Required Arguments:
  -input_dir INPUT_DIR  Directory location where the files are located.
  -input_type {separate,combined}
                        Type of input files: 'separate' for separate FASTA and GFF files, 'combined' for GFF files with embedded FASTA sequences.
  -name_split NAME_SPLIT
                        substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff').
  -output_dir OUTPUT_DIR
                        Directory for all output files.
  -output_name OUTPUT_FILE
                        Output file name.
```
