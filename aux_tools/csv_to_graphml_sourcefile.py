import pandas as pd
import glob
import networkx as nx


def process_files(file_paths):
    G = nx.Graph()

    for file_path in file_paths:
        file_name = file_path.split("/")[-1]  # Extract filename to use as attribute
        df = pd.read_csv(file_path)

        group_column = df.columns[0]
        genes_start_index = 14  # 15th column (0-based index)

        for _, row in df.iterrows():
            group_id = row[group_column]
            genes = row.iloc[genes_start_index:].dropna().values.tolist()

            # Add node attribute for file origin
            for gene in genes:
                if gene not in G.nodes:
                    G.add_node(gene, source_file=file_name)

            # Create edges between genes within the same group
            for i, gene1 in enumerate(genes):
                for j in range(i + 1, len(genes)):
                    gene2 = genes[j]
                    G.add_edge(gene1, gene2, group=group_id, source_file=file_name)

    return G


# Example usage: Load all CSV files in a directory
file_paths = glob.glob("/Users/macbookair/Git/PyamilySeq/test_data/CSVs/*.csv")
network = process_files(file_paths)

# Export to Cytoscape-compatible format
nx.write_graphml(network, "source_file_cytoscape_network.graphml")
print("Network saved as cytoscape_network.graphml")
