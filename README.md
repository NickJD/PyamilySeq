# PyamilySeq
PyamilySeq (Family Seek) is a Python tool for clustering gene sequences into families based on sequence similarity identified by tools such as CD-HIT, DIAMOND or MMseqs2.
This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter publication in NAR (https://doi.org/10.1093/nar/gkad814).

## Features

- **Clustering**: Supports input from CD-HIT formatted files as well as TSV and CSV Edge List formats.
- **Reclustering**: Allows for the addition of new sequences post-initial clustering.
- **Output**: Generates a gene 'Roary' presence-absence CSV formatted file for downstream analysis.

## Installation

PyamilySeq requires Python 3.6 or higher. Install dependencies using pip:

```bash
pip install PyamilySeq
```

## Usage

### Clustering Analysis

To perform clustering analysis:

```bash
python pyamilyseq.py -c clusters_file -f format
```

Replace `clusters_file` with the path to your clustering output file and `format` with one of: `CD-HIT`, `CSV`, or `TSV`.

### Reclustering

To add new sequences and recluster:

```bash
PyamilySeq -c clusters_file -f format --reclustered reclustered_file
```

Replace `reclustered_file` with the path to the file containing additional sequences.

## Output

PyamilySeq generates various outputs, including:

- **Gene Presence-Absence File**: This CSV file details the presence and absence of genes across genomes.

## Gene Family Groups

After analysis, PyamilySeq categorizes gene families into several groups:

- **First Core**: Gene families present in all analysed genomes initially.
- **Extended Core**: Gene families extended with additional sequences.
- **Combined Core**: Gene families combined with both initial and additional sequences.
- **Second Core**: Gene families identified only in the additional sequences.
- **Only Second Core**: Gene families exclusively found in the additional sequences.
