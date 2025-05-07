import pandas as pd
import networkx as nx


def process_files(file_path):
    G = nx.Graph()

    df = pd.read_csv(file_path)

    group_column = df.columns[0]
    genes_start_index = 14  # 15th column (0-based index)

    for _, row in df.iterrows():
        group_id = row[group_column]
        genes = row.iloc[genes_start_index:].dropna().values.tolist()

        # Only create edges if the group contains more than one gene
        if len(genes) > 1:
            for i, gene1 in enumerate(genes):
                for j in range(i + 1, len(genes)):
                    gene2 = genes[j]
                    G.add_edge(gene1, gene2, group=group_id)

    print(f"Nodes before removing isolates: {len(G.nodes)}")
    G.remove_nodes_from(list(nx.isolates(G)))  # remove isolated nodes (degree 0)

    # Now remove nodes with only one connection (degree 1)
    degree_one_nodes = [node for node, degree in G.degree() if degree == 5]
    print(f"Removing {len(degree_one_nodes)} nodes with degree 5")
    G.remove_nodes_from(degree_one_nodes)

    print(f"Nodes after pruning: {len(G.nodes)}")
    return G


# Example usage
file_path = "../test_data/CSVs/panaroo.csv"
network = process_files(file_path)

# Export to Cytoscape-compatible format
nx.write_graphml(network, "cytoscape_network_panaroo_pruned_5.graphml")
print("Network saved as cytoscape_network_panaroo_pruned.graphml")
