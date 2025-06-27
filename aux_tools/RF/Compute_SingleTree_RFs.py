import argparse

import os
import subprocess


import numpy as np
import pandas as pd
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

from ete3 import Tree

def parse_gene_presence_absence(csv_path,refound,paralogs):
    recorded_max = None  # Placeholder for max genomes, can be adjusted based on data
    df = pd.read_csv(csv_path)
    gene_id_col = df.columns[0]
    genomes = df.columns[14:]  # all genome columns after column 14
    genomes_clean = [g.replace('.combined', '').replace('_combined', '') for g in genomes] # May need to adjust based on  data
    core_genes = []

    for _, row in df.iterrows():
        gene_id = row[gene_id_col]

        if refound is False:
            # Filter out 'refound' entries from the gene columns
            filtered_genes = [g for g in row[genomes] if isinstance(g, str) and 'refound' not in g]
        else:
            filtered_genes = [g for g in row[genomes] if isinstance(g, str)]

        # Check that all genomes are present and each only has one gene (no paralogs)
        all_single = (
                len(filtered_genes) == len(genomes) and
                all(not any(sep in g for sep in [';', ',', ' ', '\t']) for g in filtered_genes)
        )
        if all_single:
            if 'refound' in filtered_genes:
                print("s")
            genes = [f"{genomes_clean[i]}|{gene}" for i, gene in enumerate(filtered_genes)]
            core_genes.append((gene_id, genes))
    return core_genes

def parse_fasta(path):
    sequences = {}
    with open(path) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = ''.join(seq_lines)
                header = line[1:].split()[0]#.split("|")[-1]  # Use last part after pipe - Maybe change this?
                header = header.replace('.combined', '').replace('_combined', '') # May need to adjust based on  data

                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            sequences[header] = ''.join(seq_lines)
    return sequences

def extract_core_sequences(core_genes, fasta_path, output_dir):
    seq_dict = parse_fasta(fasta_path)
    gene_to_file = {}

    os.makedirs(output_dir, exist_ok=True)
    for gene_id, ids in core_genes:
        found = [(sid, seq_dict[sid]) for sid in ids if sid in seq_dict]
        output_file = os.path.join(output_dir, f"{gene_id}.fa")
        with open(output_file, 'w') as out:
            for sid, seq in found:
                sid = sid.split('|')[0]  # Use only the first part before pipe
                out.write(f">{sid}\n{seq}\n")
        gene_to_file[gene_id] = output_file

    return gene_to_file

def align_sequences(input_fasta):
    aligned_fasta = input_fasta + ".aln"
    subprocess.run(["mafft", "--auto","--thread",  "4", input_fasta], stdout=open(aligned_fasta, "w"), stderr=subprocess.DEVNULL)
    return aligned_fasta

def build_tree(aligned_fasta):
    tree_file = aligned_fasta + ".tree"
    subprocess.run(["FastTree", "-quiet", "-nt", aligned_fasta], stdout=open(tree_file, "w"))
    return tree_file


def compute_rf_matrix(tree_paths, output_csv=None):
    trees = [(os.path.basename(p).split('.fa')[0], Tree(open(p).read(), format=1)) for p in tree_paths]
    n = len(trees)
    rf_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                rf_matrix[i, j] = 0
            else:
                rf, max_rf = trees[i][1].robinson_foulds(trees[j][1], unrooted_trees=True)[:2]
                rf_matrix[i, j] = rf / max_rf if max_rf > 0 else 0

    labels = [t[0] for t in trees]

    # Optional CSV output
    if output_csv:
        df = pd.DataFrame(rf_matrix, index=labels, columns=labels)
        df.to_csv(output_csv)
        print(f"\u2713 Normalised RF matrix written to: {output_csv}")

    return rf_matrix, labels


def plot_mds(matrix, labels, output=None):
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(matrix)

    plt.figure(figsize=(8, 6))
    plt.scatter(coords[:, 0], coords[:, 1], s=100, edgecolors='k')
    # for label, (x, y) in zip(labels, coords):
    #     plt.text(x, y, label, fontsize=8, ha='center', va='center')
    plt.title("MDS plot of core gene trees (normalised RF distances)")
    plt.xlabel("MDS1")
    plt.ylabel("MDS2")
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.tight_layout()
    if output:
        plt.savefig(output+'_mds_plot.png')
    else:
        plt.show()

    print("\u2713 MDS plot generated")

def plot_rf_distribution(matrix, output=None):
    # Extract upper triangle (excluding diagonal) for unique pairwise RF scores
    rf_scores = matrix[np.triu_indices_from(matrix, k=1)]
    plt.figure(figsize=(8, 6))
    plt.hist(rf_scores, bins=20, color='skyblue', edgecolor='k')
    plt.title("Distribution of pairwise normalised RF scores")
    plt.xlabel("Normalised RF score")
    plt.ylabel("Frequency")
    plt.tight_layout()
    if output:
        plt.savefig(output+'_rf_distribution.png')
    else:
        plt.show()
    print("\u2713 RF score distribution plot generated")

# Explanation:
# Multidimensional Scaling (MDS) is a dimensionality reduction technique used to visualize similarities or dissimilarities
# (e.g. distances) between data points. In this case, it is applied to a distance matrix of Robinson–Foulds (RF) values
# between phylogenetic trees constructed from alignments of core gene sequences. MDS projects the high-dimensional
# pairwise distance matrix into two dimensions, allowing you to visually assess how similar or different the gene trees are.

def main():
    parser = argparse.ArgumentParser(description="Build MSA and RF-based MDS from core gene families")
    parser.add_argument("-csv", help="Roary/Panaroo-style CSV file")
    parser.add_argument("-fasta", help="FASTA file with all gene sequences")
    parser.add_argument("-refound", choices=[True, False], default = False, type=eval, help="Allow refound genes? True/False (required)")
    parser.add_argument("-paralogs", help="Allow paralogs? True/False")
    parser.add_argument("-output", help="Output image file (e.g. mds_plot.png)")
    parser.add_argument("-tmp", default="core_temp", help="Temporary output directory")
    parser.add_argument("-rf_csv", default="rf_matrix.csv", help="CSV file to save normalized RF matrix")
    args = parser.parse_args()

    print("\u2713 Parsing CSV for core genes...")
    core_genes = parse_gene_presence_absence(args.csv,args.refound,args.paralogs)
    # need to record genomes + multicopies?

    print(f"\u2713 Extracting core gene sequences to {args.tmp}...")
    gene_to_file = extract_core_sequences(core_genes, args.fasta, args.tmp)

    tree_paths = []
    for gene_id, fasta in gene_to_file.items():
        aln = align_sequences(fasta)
        tree = build_tree(aln)
        tree_paths.append(tree)

    print("\u2713 Computing RF distance matrix...")
    matrix, labels = compute_rf_matrix(tree_paths, output_csv=args.rf_csv)
    plot_mds(matrix, labels, output=args.output)
    plot_rf_distribution(matrix,output=args.output)

if __name__ == "__main__":
    main()
