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
            genes = row.iloc[genes_start_index:].dropna().astype(str).tolist()

            # Skip groups with fewer than 10 genes
            if len(genes) < 74:
                continue

            # Add nodes with metadata: track all files the gene appears in
            for gene in genes:
                if gene not in G.nodes:
                    G.add_node(gene, source_files={file_name}, groups={group_id})
                else:
                    # Update existing node's metadata
                    G.nodes[gene]['source_files'].add(file_name)
                    G.nodes[gene]['groups'].add(group_id)

            # Create edges between genes within the same group
            for i, gene1 in enumerate(genes):
                for j in range(i + 1, len(genes)):
                    gene2 = genes[j]
                    if G.has_edge(gene1, gene2):
                        # Update existing edge's metadata
                        G[gene1][gene2]['source_files'].add(file_name)
                        G[gene1][gene2]['groups'].add(group_id)
                    else:
                        # Create a new edge
                        G.add_edge(gene1, gene2,
                                   source_files={file_name},
                                   groups={group_id})

    print(f"Nodes before final cleanup: {len(G.nodes)}")
    print(f"Edges before final cleanup: {len(G.edges)}")

    # OPTIONAL: Clean up node and edge attributes to be Cytoscape-friendly (convert sets to strings)
    for node, data in G.nodes(data=True):
        data['source_files'] = ";".join(sorted(data['source_files']))
        data['groups'] = ";".join(sorted(data['groups']))

    for u, v, data in G.edges(data=True):
        data['source_files'] = ";".join(sorted(data['source_files']))
        data['groups'] = ";".join(sorted(data['groups']))

    print(f"Final node count: {len(G.nodes)}")
    print(f"Final edge count: {len(G.edges)}")
    return G

# Example usage: Load all CSV files in a directory
file_paths = glob.glob("../test_data/CSVs/*.csv")
network = process_files(file_paths)

# Export to Cytoscape-compatible format
nx.write_graphml(network, "74_cytoscape_network_roary_panaroo_filtered_74.graphml")
print("Network saved as 74_cytoscape_network_roary_panaroo_filtered.graphml")
