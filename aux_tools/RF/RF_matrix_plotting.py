import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
from scipy.spatial import ConvexHull

def plot_mds_with_density(csv_path, output=None, title=None):
    # Load the RF matrix
    df = pd.read_csv(csv_path, index_col=0)
    matrix = df.values
    labels = df.index.tolist()

    # MDS
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
    coords = mds.fit_transform(matrix)

    # Build a DataFrame for plotting
    data = pd.DataFrame(coords, columns=["MDS1", "MDS2"])
    data["label"] = labels

    # Optional: KDE or density plot
    plt.figure(figsize=(8, 6))
    sns.kdeplot(x=data["MDS1"], y=data["MDS2"], fill=True, cmap="Blues", alpha=0.3, levels=10, thresh=0.05)

    # Scatter
    plt.scatter(data["MDS1"], data["MDS2"], s=80, edgecolor='k')
    for i, row in data.iterrows():
        plt.text(row["MDS1"], row["MDS2"], row["label"], fontsize=7, ha='center', va='center')

    # Optional: Convex hull to outline the shape
    if len(data) >= 3:
        hull = ConvexHull(coords)
        for simplex in hull.simplices:
            plt.plot(coords[simplex, 0], coords[simplex, 1], 'k--', lw=1)

    plt.title(title or "MDS plot with density")
    plt.xlabel("MDS1")
    plt.ylabel("MDS2")
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.tight_layout()
    if output:
        plt.savefig(output)
    else:
        plt.show()

plot_mds_with_density("./Roary_74_rf.csv","./Roary_74_rf_mds_density.png", title="MDS plot of Roary 74 RF distances with density")