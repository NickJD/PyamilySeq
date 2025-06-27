# PyamilySeq 
**PyamilySeq** is a Python tool for clustering gene sequences into groups based on sequence similarity identified by tools such as CD-HIT, BLAST, DIAMOND or MMseqs2.
This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter publication in NAR (https://doi.org/10.1093/nar/gkad814).

## Features
- **End-to-End**: PyamilySeq can take a directory of GFF+FASTA files, run CD-HIT for clustering and process the results.
- **Clustering input**: Supports input from CD-HIT formatted files as well as CSV and TSV node-edge lists (MMseqs2 and -outfmt 6 from BLAST/DIAMOND).
- **Reclustering**: Allows for the addition of new sequences post-initial clustering - Ensures continuity of contemporary clustering results and highlights impact of novel gene predictions.
- **'Genus Mode'**: Unlike other 'pangenome' tools, PyamilySeq can identify gene groups found across multiple genera as unique entities (see below).  
- **Output**: Generates a 'Roary/Panaroo' formatted presence-absence CSV formatted file for downstream analysis.
  - User-define species-/genus-wide gene groups - User has control over grouping parameters (core = 99/95% or min 6 genera etc).
  - Aligns representative sequences using MAFFT.
  - Output concatenated aligned sequences for tree building.
  - Optionally output sequences of each separate identified gene group.
  - Group-Splitter tool to split multi-copy gene groups.
  - Numerous additional tools to assist in the pre- and post-processing of data.

## Installation
PyamilySeq probably requires Python 3.6 or higher and the levenshtein library (https://pypi.org/project/Levenshtein/) - \
If levenshtein is not available, a Python implementation is utilised which is significantly slower. 
#### Install using pip:

```bash
pip install PyamilySeq [optionally add -U]
```
PyamilySeq is currently still under active development so expect 'regular' updates with bugfixes and new features. \
To update to the newest version add '-U' to end of the pip install command.
## Example usage: Below are two examples of running PyamilySeq in its two main modes.
```commandline
usage: PyamilySeq.py [-h] {Full,Partial} ...

PyamilySeq v1.2.0: A tool for gene clustering and analysis.

positional arguments:
  {Full,Partial}  Choose a mode: 'Full' or 'Partial'.
    Full          Full mode: PyamilySeq to cluster with CD-HIT and process output.
    Partial       Partial mode: PyamilySeq to process pre-clustered data.

options:
  -h, --help      show this help message and exit

```
### 'Full Mode': Will conduct clustering of sequences with CD-HIT as part of PyamilySeq run
```
PyamilySeq Full -output_dir .../PyamilySeq_10_AA_90_80_Full_GFFs -input_type combined -input_dir .../genomes/ -name_split_gff _combined.gff3
```
### 'Partial Mode': Will process the output of a sequence clustering from MMseqs, BLAST, DIAMOND etc.
```
PyamilySeq Partial -clustering_format CD-HIT -cluster_file .../all_10_combined_pep_CD-HIT_90_80.clstr  -original_fasta .../all_10_combined_pep.fasta -output_dir .../PyamilySeq_10_AA_90_80_Partial -write_groups 99 -align 
```


#### Note: using a '-clustering_format' other than the default CD-HIT, requires input to be two in two columns as below (Same format as MMseqs2 tsv and BLAST outfmt 6) - Genome name and sequence name are separated by '|'.
```
Escherichia_coli_110957|ENSB_lL-zIKt-gh0oSno	Escherichia_coli_110957|ENSB_lL-zIKt-gh0oSno
Escherichia_coli_110957|ENSB_lL-zIKt-gh0oSno	Escherichia_coli_113290|ENSB_2fj4rJ8e8Z9PNdX
Escherichia_coli_110957|ENSB_lL-zIKt-gh0oSno	Escherichia_coli_b185|ENSB_G_PVe28-ej8q-3S
Escherichia_coli_110957|ENSB_TIZS9kbTvShDvyX	Escherichia_coli_110957|ENSB_TIZS9kbTvShDvyX
```
### Example output:
```
Running PyamilySeq v1.2.0
Calculating Groups
Number of Genomes: 10
Gene Groups
First_core_99: 2994
First_core_95: 0
First_core_15: 3266
First_core_0: 5466
Total Number of First Gene Groups (Including Singletons): 11726
Outputting gene_presence_absence file
Outputting gene group FASTA files
Combined FASTA file saved to: ../combined_group_sequences_dna.fasta
Processing gene group alignment
Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq
Please report any issues to: https://github.com/NickJD/PyamilySeq/issues
```


## Reclustering:
### Reclustering can be used to see where additional sequences/genes lay in relation to a contemporary pangenome/gene grouping.
```
PyamilySeq Partial -clustering_format CD-HIT -cluster_file .../all_10_combined_pep_CD-HIT_90_80.clstr -reclustered .../all_10_combined_pep_CD-HIT_90_80_AND_StORFs_CD-HIT_90_80.clstr -original_fasta .../all_10_combined_pep_AND_StORFs.fasta -output_dir .../PyamilySeq_10_AA_90_80_Partial_Reclustered_StORFs -write_groups 99 -align 
```
#### As can be seen below, the additional sequences recovered by the StORF-Reporter annotation tool have 'extended' contemporary or created entirely new gene groups. 'First' corresponds to the groups identified from the first clustering round and 'Second' for the second. In 'reclustering' mode, First_core_# groups are unaffected thus retaining the initial grouping information. 
```commandline
Number of Genomes: 10
Gene Groups
First_core_99: 2994
First_core_95: 0
First_core_15: 3266
First_core_0: 5466
extended_core_99: 3
extended_core_95: 0
extended_core_15: 49
extended_core_0: 0
combined_core_99: 0
combined_core_95: 0
combined_core_15: 3
combined_core_0: 0
Second_core_99: 0
Second_core_95: 0
Second_core_15: 20
Second_core_0: 39
only_Second_core_99: 768
only_Second_core_95: 0
only_Second_core_15: 4472
only_Second_core_0: 8395
Total Number of First Gene Groups (Including Singletons): 11726
Total Number of Second Gene Groups (Including Singletons): 25359
Total Number of First Gene Groups That Had Additional Second Sequences But Not New Genomes: 5
```

## PyamilySeq is separated into two main 'run modes', Full and Partial. They each have their own set of required and optional arguments.
### PyamilySeq - Full Menu: 
```
usage: PyamilySeq.py Full [-h] -output_dir OUTPUT_DIR -input_type {separate,combined,fasta} [-input_dir INPUT_DIR] [-input_fasta INPUT_FASTA] [-name_split_gff NAME_SPLIT_GFF] [-name_split_fasta NAME_SPLIT_FASTA] [-sequence_type {AA,DNA}] [-gene_ident GENE_IDENT] [-c PIDENT] [-s LEN_DIFF] [-fast_mode]
                          [-group_mode {Species,Genus}] [-species_groups SPECIES_GROUPS] [-genus_groups GENUS_GROUPS] [-write_groups WRITE_GROUPS] [-write_individual_groups] [-align] [-align_aa] [-no_gpa] [-M MEM] [-T THREADS] [-verbose] [-v]

options:
  -h, --help            show this help message and exit
  -output_dir OUTPUT_DIR
                        Directory for all output files.
  -input_type {separate,combined,fasta}
                        Type of input files: 'separate' for matching FASTA and GFF files, 'combined' for GFF+FASTA, or 'fasta' for a prepared FASTA file.
  -input_dir INPUT_DIR  Directory containing GFF/FASTA files - Use with -input_type separate/combined.
  -input_fasta INPUT_FASTA
                        Input FASTA file - Use with - input_type fasta.
  -name_split_gff NAME_SPLIT_GFF
                        Substring to split filenames and extract genome names for gff files (e.g., '_combined.gff3') - Use with -input_type separate/combined.
  -name_split_fasta NAME_SPLIT_FASTA
                        Substring to split filenames and extract genome names for fasta files if named differently to paired gff files (e.g., '_dna.fasta') - Use with -input_type separate/combined.
  -sequence_type {AA,DNA}
                        Clustering mode: 'DNA' or 'AA'.
  -gene_ident GENE_IDENT
                        Gene identifiers to extract sequences (e.g., 'CDS, tRNA').
  -c PIDENT             Sequence identity threshold for clustering (default: 0.90) - CD-HIT parameter '-c'.
  -s LEN_DIFF           Length difference threshold for clustering (default: 0.80) - CD-HIT parameter '-s'.
  -fast_mode            Enable fast mode for CD-HIT (not recommended) - CD-HIT parameter '-g'.
  -group_mode {Species,Genus}
                        Grouping mode: 'Species' or 'Genus'.
  -species_groups SPECIES_GROUPS
                        Gene groupings for 'Species' mode (default: '99,95,15').
  -genus_groups GENUS_GROUPS
                        Gene groupings for 'Genus' mode (default: '1-10').
  -write_groups WRITE_GROUPS
                        Output gene groups as a single FASTA file (specify levels: e.g., '-w 99,95'). - triggers '-wig'.
  -write_individual_groups
                        Output individual FASTA files for each group.
  -align                Align and concatenate sequences for 'core' groups specified with '-w'.
  -align_aa             Align sequences as amino acids.
  -no_gpa               Skip creation of gene_presence_absence.csv.
  -M MEM                Memory allocation for clustering (MB) - CD-HIT parameter '-M'.
  -T THREADS            Number of threads for clustering/alignment - CD-HIT parameter '-T' | MAFFT parameter '--thread'.
  -verbose              Print verbose output.
  -v, --version         Print version number and exit.

```
### PyamilySeq - Partial Menu: 
```commandline
usage: PyamilySeq.py Partial [-h] -clustering_format {CD-HIT,MMseqs,BLAST} -cluster_file CLUSTER_FILE -original_fasta ORIGINAL_FASTA -output_dir OUTPUT_DIR [-reclustered RECLUSTERED] [-seq_tag SEQUENCE_TAG] [-group_mode {Species,Genus}] [-species_groups SPECIES_GROUPS] [-genus_groups GENUS_GROUPS]
                             [-write_groups WRITE_GROUPS] [-write_individual_groups] [-align] [-align_aa] [-no_gpa] [-M MEM] [-T THREADS] [-verbose] [-v]

options:
  -h, --help            show this help message and exit
  -clustering_format {CD-HIT,MMseqs,BLAST}
                        Clustering format used: CD-HIT, MMseqs2, or BLAST.
  -cluster_file CLUSTER_FILE
                        Cluster file containing pre-clustered groups from CD-HIT, MMseqs, BLAST etc.
  -original_fasta ORIGINAL_FASTA
                        FASTA file used in pre-clustering (Provide sequences in DNA form).
  -output_dir OUTPUT_DIR
                        Directory for all output files.
  -reclustered RECLUSTERED
                        Clustering output file from a second round of clustering.
  -seq_tag SEQUENCE_TAG
                        Tag for distinguishing reclustered sequences.
  -group_mode {Species,Genus}
                        Grouping mode: 'Species' or 'Genus'.
  -species_groups SPECIES_GROUPS
                        Gene groupings for 'Species' mode (default: '99,95,15').
  -genus_groups GENUS_GROUPS
                        Gene groupings for 'Genus' mode (default: '1-10').
  -write_groups WRITE_GROUPS
                        Output gene groups as a single FASTA file (specify levels: e.g., '-w 99,95'). - triggers '-wig'.
  -write_individual_groups
                        Output individual FASTA files for each group.
  -align                Align and concatenate sequences for 'core' groups specified with '-w'.
  -align_aa             Align sequences as amino acids.
  -no_gpa               Skip creation of gene_presence_absence.csv.
  -M MEM                Memory allocation for clustering (MB) - CD-HIT parameter '-M'.
  -T THREADS            Number of threads for clustering/alignment - CD-HIT parameter '-T' | MAFFT parameter '--thread'.
  -verbose              Print verbose output.
  -v, --version         Print version number and exit.

```

## Seq-Combiner: This tool is provided to enable the pre-processing of multiple GFF/FASTA files together ready to be clustered by the user.
### Example:
```bash
Seq-Combiner -input_dir .../test_data/genomes -name_split_gff .gff3 -output_dir .../test_data/genomes -output_name combine_fasta_seqs.fa -input_type combined
```
### Seq-Combiner Menu:
```
usage: Seq_Combiner.py [-h] -input_dir INPUT_DIR -input_type {separate,combined,fasta} [-name_split_gff NAME_SPLIT_GFF] [-name_split_fasta NAME_SPLIT_FASTA] -output_dir OUTPUT_DIR -output_name OUTPUT_FILE [-gene_ident GENE_IDENT] [-translate] [-v]

PyamilySeq v1.2.0: Seq-Combiner - A tool to extract sequences from GFF/FASTA files and prepare them for PyamilySeq.

options:
  -h, --help            show this help message and exit

Required Arguments:
  -input_dir INPUT_DIR  Directory location where the files are located.
  -input_type {separate,combined,fasta}
                        Type of input files: "separate" for separate FASTA and GFF files, "combined" for GFF files with embedded FASTA sequences and "fasta" for combining multiple FASTA files together.
  -name_split_gff NAME_SPLIT_GFF
                        Substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff'). - Not needed with -input_type fasta
  -name_split_fasta NAME_SPLIT_FASTA
                        Substring used to split filenames and extract genome names for fasta files if named differently to paired gff files (e.g., '_dna.fasta').
  -output_dir OUTPUT_DIR
                        Directory for all output files.
  -output_name OUTPUT_FILE
                        Output file name.

Optional Arguments:
  -gene_ident GENE_IDENT
                        Default - "CDS": Identifier used for extraction of sequences such as "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo" - Not compatible with "fasta" input mode.
  -translate            Default - False: Translate extracted sequences to their AA counterpart? - appends _aa.fasta to given output_name

Misc Arguments:
  -v, --version         Print out version number and exit


```

## Group-Splitter: This tool can split multi-copy gene groups using CD-HIT after initial PyamilySeq analysis.
### Example:
```bash
Group-Splitter -genome_num 10 -input_fasta .../test/species/ -output_dir .../test/species/ -sequence_type AA
```
### Group-Splitter Menu:
```
usage: Group_Splitter.py [-h] -input_fasta INPUT_FASTA -sequence_type {AA,DNA}
                         -genome_num GENOME_NUM -output_dir OUTPUT_DIR
                         [-groups GROUPS] [-group_threshold GROUP_THRESHOLD]
                         [-c PIDENT] [-s LEN_DIFF] [-T CLUSTERING_THREADS]
                         [-M CLUSTERING_MEMORY] [-no_delete_temp_files]
                         [-verbose] [-v]

PyamilySeq v1.2.0: Group-Splitter - A tool to split multi-copy gene groups
identified by PyamilySeq.

options:
  -h, --help            show this help message and exit

Required Parameters:
  -input_fasta INPUT_FASTA
                        Input FASTA file containing gene groups.
  -sequence_type {AA,DNA}
                        Default - DNA: Are groups "DNA" or "AA" sequences?
  -genome_num GENOME_NUM
                        The total number of genomes must be provide
  -output_dir OUTPUT_DIR
                        Output directory.

Regrouping Parameters:
  -groups GROUPS        Default - auto: Detect groups to be split (see
                        -group_threshold). Provide "-groups 1,2,3,4" with
                        group IDs to split specific groups.
  -group_threshold GROUP_THRESHOLD
                        Minimum percentage of genomes with multi-copy
                        (default: 80.0) - Does not work with "-groups"

CD-HIT Reclustering Parameters:
  -c PIDENT             Sequence identity threshold (default: 0.8) - Probably
                        should be higher than what was used in initial
                        clustering.
  -s LEN_DIFF           Length difference cutoff (default: 0.20) - Often the
                        most impactful parameter to split 'multi-copy' gene
                        groups.
  -T CLUSTERING_THREADS
                        Number of threads for clustering (default: 4)
  -M CLUSTERING_MEMORY  Memory limit in MB for clustering (default: 2000)

Misc Parameters:
  -no_delete_temp_files
                        Default: Delete all temporary files after processing.
  -verbose              Print verbose output.
  -v, --version         Print out version number and exit

```

## Cluster-Summary menu: This tool can be used to summarise CD-HIT .clstr files:
### Example:
```bash
Cluster-Summary -genome_num 10 -input_clstr .../test_data/species/E-coli/E-coli_extracted_pep_cd-hit_80.clstr -output_tsv .../test_data/species/E-coli/E-coli_extracted_pep_cd-hit_80_Summary.tsv
```
### Cluster-Summary Menu: 
```
usage: Cluster_Summary.py [-h] -input_clstr INPUT_CLSTR -output OUTPUT -genome_num GENOME_NUM
                          [-output_dir OUTPUT_DIR] [-verbose] [-v]

PyamilySeq v1.2.0: Cluster-Summary - A tool to summarise CD-HIT clustering files.

options:
  -h, --help            show this help message and exit

Required Parameters:
  -input_clstr INPUT_CLSTR
                        Input CD-HIT .clstr file
  -output OUTPUT        Output TSV file to store cluster summaries - Will add '.tsv' if not
                        provided by user
  -genome_num GENOME_NUM
                        The total number of genomes must be provide

Optional Arguments:
  -output_dir OUTPUT_DIR
                        Default: Same as input file

Misc Parameters:
  -verbose              Print verbose output.
  -v, --version         Print out version number and exit

```

### All example input and output data can be found  in the 'test_data' directory.