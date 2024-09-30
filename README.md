# PyamilySeq - !BETA!
**PyamilySeq** is a Python tool for clustering gene sequences into groups based on sequence similarity identified by tools such as CD-HIT, BLAST, DIAMOND or MMseqs2.
This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter publication in NAR (https://doi.org/10.1093/nar/gkad814).

## Features
- **End-to-End**: PyamilySeq can take a directory of GFF+FASTA files, run CD-HIT for clustering and process the results.
- **Clustering input**: Supports input from CD-HIT formatted files as well as CSV and TSV edge lists (MMseqs2 and -outfmt 6 from BLAST/DIAMOND).
- **Reclustering**: Allows for the addition of new sequences post-initial clustering - Ensures continuity of contemporary clustering results and highlights impact of novel gene predictions.
- **'Genus Mode'**: Unlike other 'pangenome' tools, PyamilySeq can identify gene groups found across multiple genera as unique entities (see below).  
- **Output**: Generates a 'Roary/Panaroo' formatted presence-absence CSV formatted file for downstream analysis.
  - User-define species-/genus-wide gene groups - User has control over grouping parameters (core = 99/95% or min 6 genera etc).
  - Aligns representative sequences using MAFFT.
  - Output concatenated aligned sequences for tree building.
  - Optionally output sequences of each separate identified gene group.

## Installation
PyamilySeq probably requires Python 3.6 or higher. Install using pip:

```bash
pip install PyamilySeq
```

## Example usage: Below are two examples of running PyamilySeq in its two main modes.
### 'Full Mode': Will conduct clustering of sequences with CD-HIT as part of PyamilySeq run
```
PyamilySeq -run_mode Full -group_mode Species -clustering_format CD-HIT  -output_dir .../test_data/testing/Full 
-input_type combined -input_dir .../test_data/genomes -name_split _combined.gff3 -pid 0.95 -len_diff 0.80 
-gpa -a -w 99
```
### 'Partial Mode': Will take the output of a sequence clustering.
```
PyamilySeq -run_mode Partial -group_mode Species -clustering_format TSV -output_dir .../test_data/Species/testing/Partial 
-cluster_file .../test_data/Species/MMseqs2/combined_Ensmbl_pep_cluster.tsv 
-original_fasta .../test_data/species/combined_Ensmbl_cds.fasta -gpa -a -w 99 -verbose 

```
#### Note: '-clustering_format TSV/CSV' requires input to be two in two columns as below (Same format as MMseqs2 tsv) - Genome name and sequence name are separated by '|'.
```
Escherichia_coli_110957|ENSB:lL-zIKt-gh0oSno	Escherichia_coli_110957|ENSB:lL-zIKt-gh0oSno
Escherichia_coli_110957|ENSB:lL-zIKt-gh0oSno	Escherichia_coli_113290|ENSB:2fj4rJ8e8Z9PNdX
Escherichia_coli_110957|ENSB:lL-zIKt-gh0oSno	Escherichia_coli_b185|ENSB:G_PVe28-ej8q-3S
Escherichia_coli_110957|ENSB:TIZS9kbTvShDvyX	Escherichia_coli_110957|ENSB:TIZS9kbTvShDvyX
```
### Example output:
```
Running PyamilySeq v0.8.1
Calculating Groups
Gene Groups:
First_core_99: 2682
First_core_95: 0
First_core_15: 3789
First_core_0: 6469
Total Number of First Gene Groups (Including Singletons): 12940
Outputting gene_presence_absence file
Outputting gene group FASTA files
Processing gene group alignment
Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq
Please report any issues to: https://github.com/NickJD/PyamilySeq/issues
```
## Genus mode: 
### In addition to "Species mode" (see above) which reports gene groups the same as pangenome tools such as Roary and Panaroo, Genus mode reports gene groups identified across multiple genera.
#### Example:
```
PyamilySeq -run_mode Partial -group_mode Genus -clustering_format CD-HIT -output_dir .../test_data/genus/testing/
 -cluster_file .../test_data/genus/CD-HIT/combined_cds_cd-hit_80_60.clstr -gpa 
```
```commandline
Running PyamilySeq v0.8.1
Calculating Groups
Genus Groups:
First_genera_1:	28549
First_genera_2:	12
First_genera_3:	0
First_genera_>:	0
Total Number of First Gene Groups (Including Singletons): 28561
Outputting gene_presence_absence file
Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq
Please report any issues to: https://github.com/NickJD/PyamilySeq/issues
#####
```

## Reclustering:
### Reclustering can be used to see where additional sequences/genes lay in relation to a contemporary pangenome/gene grouping.
```
PyamilySeq -run_mode Partial -group_mode Species -clustering_format CD-HIT -output_dir .../test_data/species/CD-HIT/testing 
-cluster_file .../test_data/species/CD-HIT/E-coli_extracted_cds_cd-hit_80_60.clstr -gpa 
-reclustered .../test_data/species/CD-HIT/E-coli_extracted_cds_cd-hit_80_60_And_StORFs_cds_80_60.clstr
```
#### As can be seen below, the additional sequences recovered by the StORF-Reporter annotation tool have 'extended' contemporary or created entirely new gene groups. 'First' corresponds to the groups identified from the first clustering round and 'Second' for the second. In 'reclustering' mode, First_core_# groups are unaffected thus retaining the initial grouping information. 
```commandline
Calculating Groups
Gene Groups:
First_core_99: 587
First_core_95: 1529
First_core_15: 3708
First_core_0: 29992
extended_core_99: 29
extended_core_95: 67
extended_core_15: 431
extended_core_0: 1331
combined_core_99: 2
combined_core_95: 4
combined_core_15: 5
combined_core_0: 4
Second_core_99: 0
Second_core_95: 6
Second_core_15: 172
Second_core_0: 1825
only_Second_core_99: 53
only_Second_core_95: 493
only_Second_core_15: 3806
only_Second_core_0: 27569
Total Number of First Gene Groups (Including Singletons): 35816
Total Number of Second Gene Groups (Including Singletons): 67728
Total Number of First Gene Groups That Had Additional Second Sequences But Not New Genomes: 136
Outputting gene_presence_absence file
Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq
Please report any issues to: https://github.com/NickJD/PyamilySeq/issues
#####
```

## PyamilySeq - Menu: 
### PyamilySeq is separated into two main 'run modes', Full and Partial. They each have their own set of required and optional arguments. 
```
Running PyamilySeq v0.8.1
usage: PyamilySeq.py [-h] -run_mode {Full,Partial} -group_mode {Species,Genus} -clustering_format {CD-HIT,TSV,CSV} -output_dir OUTPUT_DIR
                     [-input_type {separate,combined}] [-input_dir INPUT_DIR] [-name_split NAME_SPLIT] [-sequence_type {AA,DNA}] [-gene_ident GENE_IDENT]
                     [-pid PIDENT] [-len_diff LEN_DIFF] [-mem CLUSTERING_MEMORY] [-t CLUSTERING_THREADS] [-cluster_file CLUSTER_FILE]
                     [-reclustered RECLUSTERED] [-seq_tag SEQUENCE_TAG] [-core_groups CORE_GROUPS] [-genus_groups GENUS_GROUPS] [-w WRITE_GROUPS] [-a]
                     [-original_fasta ORIGINAL_FASTA] [-gpa] [-verbose] [-v]

PyamilySeq v0.8.1: A tool that groups genes into unique clusters.

options:
  -h, --help            show this help message and exit

Required Arguments:
  -run_mode {Full,Partial}
                        Run Mode: Should PyamilySeq be run in "Full" or "Partial" mode?
  -group_mode {Species,Genus}
                        Group Mode: Should PyamilySeq be run in "Species" or "Genus" mode?
  -clustering_format {CD-HIT,TSV,CSV}
                        Clustering format to use: CD-HIT or TSV (MMseqs2, BLAST, DIAMOND) / CSV edge-list file (Node1 Node2).
  -output_dir OUTPUT_DIR
                        Directory for all output files.

Full-Mode Arguments - Required when "-run_mode Full" is used:
  -input_type {separate,combined}
                        Type of input files: 'separate' for separate FASTA and GFF files, 'combined' for GFF files with embedded FASTA sequences.
  -input_dir INPUT_DIR  Directory containing GFF/FASTA files.
  -name_split NAME_SPLIT
                        substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff').
  -sequence_type {AA,DNA}
                        Default - DNA: Should clustering be performed in "DNA" or "AA" mode?
  -gene_ident GENE_IDENT
                        Identifier used for extraction of sequences such as
                        "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"
  -pid PIDENT           Default 0.95: Pident threshold for clustering.
  -len_diff LEN_DIFF    Default 0.80: Minimum length difference between clustered sequences - (-s) threshold for CD-HIT clustering.

Clustering Runtime Arguments - Optional when "-run_mode Full" is used:
  -mem CLUSTERING_MEMORY
                        Default 4000: Memory to be allocated for clustering (in MBs).
  -t THREADS            Default 8: Threads to be allocated for clustering
                        and/or alignment.


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
  -w WRITE_GROUPS       Default - No output: Output sequences of identified groups (provide levels at which to output - Species "-w 99,95" Genus "-w 2,3" -
                        Must provide FASTA file with -original_fasta if in Partial run mode.
  -a                    Default - No output: SLOW! (Only works for Species mode) Output aligned and concatinated sequences of identified groups -provide
                        group levels at which to output "-w 99,95" - Must provide FASTA file with -original_fasta in Partial run mode.
  -original_fasta ORIGINAL_FASTA
                        FASTA file to use in conjunction with "-w" or "-a" when running in Partial Mode.
  -gpa                  Default - False: If selected, a Roary/Panaroo formatted gene_presence_absence.csv will be created - Required for Coinfinder and
                        other downstream tools

Misc:
  -verbose              Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit
```





## Seq-Combiner: This tool is provided to enable the pre-processing of multiple GFF/FASTA files together ready to be clustered by the user.
### Example:
```bash
Seq-Combiner -input_dir .../test_data/genomes -name_split _combined.gff3 -output_dir.../test_data -output_name combine_fasta_seqs.fa -input_type combined
```
### Seq-Combiner Menu:
```
usage: Seq_Combiner.py [-h] -input_dir INPUT_DIR -input_type {separate,combined,fasta} -name_split NAME_SPLIT -output_dir OUTPUT_DIR -output_name OUTPUT_FILE [-gene_ident GENE_IDENT] [-translate] [-v]

Seq-Combiner v0.8.1: A tool to extract sequences from GFF/FASTA files.

options:
  -h, --help            show this help message and exit

Required Arguments:
  -input_dir INPUT_DIR  Directory location where the files are located.
  -input_type {separate,combined,fasta}
                        Type of input files: "separate" for separate FASTA and GFF files, "combined" for GFF files with embedded FASTA sequences and "fasta" for combining multiple FASTA files together.
  -name_split NAME_SPLIT
                        substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff').
  -output_dir OUTPUT_DIR
                        Directory for all output files.
  -output_name OUTPUT_FILE
                        Output file name.

Optional Arguments:
  -gene_ident GENE_IDENT
                        Default - "CDS": Identifier used for extraction of sequences such as "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo" - Not compatible with "fasta" input mode.
  -translate            Default - False: Translate extracted sequences to their AA counterpart?

Misc Arguments:
  -v                    Print out version number and exit


```

### Group-Splitter menu: 

```
usage: Group_Splitter.py [-h] -input_fasta INPUT_FASTA -output_dir OUTPUT_DIR [-pident PIDENT] [-len_diff LEN_DIFF] [-clustering_threads CLUSTERING_THREADS]
                         [-clustering_memory CLUSTERING_MEMORY] [-percent_threshold PERCENT_THRESHOLD] [-verbose] [-delete_temp_files] [-v]

Group-Splitter: v0.8.1: A tool to split "paralogous" groups identified by PyamilySeq.

options:
  -h, --help            show this help message and exit

Required Arguments:
  -input_fasta INPUT_FASTA
                        Input FASTA file containing gene groups.
  -sequence_type {AA,DNA}
                        Default - DNA: Are groups "DNA" or "AA" sequences?
  -output_dir OUTPUT_DIR
                        Output directory.

Optional Arguments:
  -pident PIDENT        Sequence identity threshold (default: 0.9)
  -len_diff LEN_DIFF    Length difference threshold (default: 0.05)
  -clustering_threads CLUSTERING_THREADS
                        Number of threads for clustering (default: 4)
  -clustering_memory CLUSTERING_MEMORY
                        Memory limit in MB for clustering (default: 2000)
  -percent_threshold PERCENT_THRESHOLD
                        Minimum percentage of genomes with paralogs (default: 80.0)
  -verbose              Print verbose output.
  -delete_temp_files    Delete all temporary files after processing.

Misc Arguments:
  -v                    Print out version number and exit
```

### All example input and output data can be found  in the 'test_data' directory.