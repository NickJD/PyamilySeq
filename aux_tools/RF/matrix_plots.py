# === plot_rf_distributions.py ===

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def extract_upper_triangular(matrix):
    return matrix.where(np.triu(np.ones(matrix.shape), k=1).astype(bool)).stack().values


def plot_kde(data, output):
    plt.figure(figsize=(6, 4))
    sns.kdeplot(data, fill=True, linewidth=1.5)
    plt.xlabel("Normalized RF Distance", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title("Distribution of Normalized RF Distances", fontsize=14)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output or "rf_kde_plot.pdf")
    plt.close()


def plot_violin(data, output):
    plt.figure(figsize=(4, 6))
    sns.violinplot(y=data, inner="box", linewidth=1.2)
    plt.ylabel("Normalized RF Distance", fontsize=12)
    plt.title("Violin Plot of RF Distances", fontsize=14)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output or "rf_violin_plot.pdf")
    plt.close()


def plot_histogram(data, output):
    plt.figure(figsize=(6, 4))
    plt.hist(data, bins=30, color="steelblue", edgecolor="black")
    plt.xlabel("Normalized RF Distance", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.title("Histogram of Normalized RF Distances", fontsize=14)
    plt.tight_layout()
    plt.savefig(output or "rf_histogram_plot.pdf")
    plt.close()


def plot_heatmap(matrix, labels, output):
    df = pd.DataFrame(matrix, index=labels, columns=labels)
    plt.figure(figsize=(10, 8))
    sns.heatmap(df, cmap="viridis", square=True, xticklabels=False, yticklabels=False, cbar_kws={"label": "Normalized RF"})
    plt.title("Heatmap of Normalized RF Distances", fontsize=14)
    plt.tight_layout()
    plt.savefig(output or "rf_heatmap_plot.pdf")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Generate distribution plots from RF matrix")
    parser.add_argument("matrix", help="CSV file with normalized RF matrix")
    parser.add_argument("--prefix", help="Prefix for output plots (default: derived from matrix filename)")
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, index_col=0)
    values = extract_upper_triangular(df)
    labels = df.index.tolist()
    prefix = args.prefix or args.matrix.rsplit(".", 1)[0]

    plot_kde(values, f"{prefix}_kde.pdf")
    plot_violin(values, f"{prefix}_violin.pdf")
    plot_histogram(values, f"{prefix}_histogram.pdf")
    plot_heatmap(df.values, labels, f"{prefix}_heatmap.pdf")
    print("\u2713 All distribution plots saved")


if __name__ == "__main__":
    main()
