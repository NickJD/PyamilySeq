import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from matplotlib.cm import ScalarMappable
import argparse

def load_rf_matrix(path):
    df = pd.read_csv(path, index_col=0)
    labels = df.index.tolist()
    matrix = df.values.astype(float)
    return matrix, labels

def plot_mds(matrix, labels, title="RF Distance MDS Projection", output=None):
    # MDS projection to 2D
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(matrix)

    # Color by rank index
    rf_mean = matrix.mean(axis=1)  # mean RF distance of each tree to all others
    colors = rf_mean  # continuous value between 0 and 1
    norm = plt.Normalize(vmin=0, vmax=1)
    cmap = plt.cm.inferno  # or any perceptually uniform colormap

    #color_vals = np.linspace(1, len(labels), len(labels))
    #norm = plt.Normalize(vmin=1, vmax=len(labels))
    #cmap = plt.cm.viridis
    #colors = cmap(norm(color_vals))

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=30)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])

    # Color bar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label="Index")

    if output:
        plt.savefig(output, bbox_inches="tight", dpi=300)
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot 2D MDS from RF distance matrix")
    parser.add_argument("matrix", help="CSV file of pairwise RF distances")
    parser.add_argument("-o", "--output", help="Output image file (e.g. plot.png)")
    args = parser.parse_args()

    matrix, labels = load_rf_matrix(args.matrix)
    plot_mds(matrix, labels, output=args.output)

if __name__ == "__main__":
    main()

