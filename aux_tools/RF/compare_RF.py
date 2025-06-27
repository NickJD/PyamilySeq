# === compare_rf.py ===

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.stats.distance import DistanceMatrix, mantel
from scipy.stats import ks_2samp, mannwhitneyu


def mean_pairwise_distance(matrix):
    upper = matrix[np.triu_indices_from(matrix, k=1)]
    return np.mean(upper), np.std(upper), upper


def load_distance_matrix(path):
    df = pd.read_csv(path, index_col=0)
    mat = df.values
    return df.columns.tolist(), mat


def main():
    parser = argparse.ArgumentParser(description="Compare RF distance matrices from multiple tools")
    parser.add_argument("matrices", nargs='+', help="CSV files of RF matrices from different tools")
    parser.add_argument("--labels", nargs='+', help="Labels for tools (must match number of matrices)")
    parser.add_argument("--output", default="rf_summary.txt", help="Output summary file")
    parser.add_argument("--plot", default="rf_distributions.png", help="Output comparison plot")
    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.matrices):
        raise ValueError("Number of labels must match number of matrices")

    tool_labels = args.labels if args.labels else [f"Tool{i+1}" for i in range(len(args.matrices))]

    stats = []
    combined_df = pd.DataFrame()

    for label, path in zip(tool_labels, args.matrices):
        labels, mat = load_distance_matrix(path)
        mean, std, flat = mean_pairwise_distance(mat)
        stats.append((label, mean, std))
        df = pd.DataFrame({"RF": flat, "Tool": label})
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    with open(args.output, "w") as out:
        out.write("Mean Pairwise RF Distances:\n")
        for label, mean, std in stats:
            out.write(f"{label}: {mean:.4f} (std: {std:.4f})\n")

        out.write("\nPairwise Comparisons:\n")
        for i in range(len(tool_labels)):
            for j in range(i+1, len(tool_labels)):
                rf1 = combined_df[combined_df['Tool'] == tool_labels[i]]['RF']
                rf2 = combined_df[combined_df['Tool'] == tool_labels[j]]['RF']
                ks_stat, ks_p = ks_2samp(rf1, rf2)
                u_stat, u_p = mannwhitneyu(rf1, rf2, alternative='two-sided')
                out.write(f"{tool_labels[i]} vs {tool_labels[j]}:\n")
                out.write(f"  KS statistic = {ks_stat:.4f}, p = {ks_p:.4g}\n")
                out.write(f"  Mann–Whitney U = {u_stat:.4f}, p = {u_p:.4g}\n")

    plt.figure(figsize=(10, 6))
    sns.violinplot(data=combined_df, x="Tool", y="RF", inner=None, palette="Set2")
    sns.boxplot(data=combined_df, x="Tool", y="RF", width=0.2, showcaps=True, boxprops={'facecolor':'none'}, showfliers=False)
    plt.title("Distribution of Pairwise RF Distances by Tool")
    plt.ylabel("Normalised RF Distance")
    plt.tight_layout()
    plt.savefig(args.plot)
    print(f"✓ Summary saved to {args.output}\n✓ Plot saved to {args.plot}")


if __name__ == "__main__":
    main()
