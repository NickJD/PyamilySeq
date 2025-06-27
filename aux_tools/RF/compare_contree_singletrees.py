import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ete3 import Tree


def compare_single_to_concat(concat_path, single_tree_dir):
    concat_tree = Tree(concat_path, format=1, quoted_node_names=True)
    rf_scores = []

    for filename in os.listdir(single_tree_dir):
        if filename.endswith(".nwk") or filename.endswith(".tree"):
            tree_path = os.path.join(single_tree_dir, filename)
            single_tree = Tree(tree_path, format=1, quoted_node_names=True)
            rf, max_rf, *_ = single_tree.robinson_foulds(concat_tree, unrooted_trees=True)
            norm_rf = rf / max_rf if max_rf else 0
            rf_scores.append((filename, norm_rf))

    return pd.DataFrame(rf_scores, columns=["Gene", "Normalised_RF"])


def plot_violin(data, output=None):
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=data, x="Tool", y="Normalised_RF", inner="box", palette="muted")
    plt.title("Distribution of Single Gene Tree vs Concatenated Tree RF Distances")
    plt.tight_layout()
    if output:
        plt.savefig(output)
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Compare single-gene trees to a concatenated tree using RF distance.")
    parser.add_argument("-i", "--inputs", nargs=3, action="append", metavar=('CONCAT_TREE', 'SINGLE_TREE_DIR', 'LABEL'),
                        help="Provide triplets: <concat_tree.nwk> <directory_with_single_trees> <label>")
    parser.add_argument("-o", "--output_plot", help="Output image filename (e.g., plot.png)")
    parser.add_argument("-c", "--output_csv", help="CSV to write RF values to")
    args = parser.parse_args()


    all_data = []
    for concat_tree, tree_dir, label in args.inputs:
        df = compare_single_to_concat(concat_tree, tree_dir)
        df["Tool"] = label
        all_data.append(df)

    full_df = pd.concat(all_data, ignore_index=True)

    if args.output_csv:
        full_df.to_csv(args.output_csv, index=False)
        print(f"✓ RF data saved to {args.output_csv}")

    plot_violin(full_df, output=args.output_plot)
    print(f"✓ Plot saved to {args.output_plot}")




if __name__ == "__main__":
    main()
