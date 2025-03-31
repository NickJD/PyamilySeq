import pandas as pd
import glob
import networkx as nx


def process_files(file_paths):
    G = nx.Graph()

    for file_path in file_paths:
        df = pd.read_csv(file_path)

        # Extract relevant columns (Group ID + Genes)
        group_column = df.columns[0]
        genes_start_index = 14  # 15th column (0-based index)

        for _, row in df.iterrows():
            group_id = row[group_column]
            genes = row.iloc[genes_start_index:].dropna().values.tolist()

            # Create edges between genes within the same group
            for i, gene1 in enumerate(genes):
                for j in range(i + 1, len(genes)):
                    gene2 = genes[j]
                    G.add_edge(gene1, gene2, group=group_id)

    return G


# Example usage: Load all CSV files in a directory
file_paths = glob.glob("/Users/macbookair/Git/PyamilySeq/test_data/CSVs/*.csv")
network = process_files(file_paths)

# Export to Cytoscape-compatible format
nx.write_graphml(network, "cytoscape_network.graphml")
print("Network saved as cytoscape_network.graphml")
