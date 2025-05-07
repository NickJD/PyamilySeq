import pandas as pd
import glob
from collections import defaultdict

def collect_74_groups_distinct(file_paths):
    # Key: (Group_ID, sorted_gene_tuple)
    group_data = defaultdict(lambda: {'source_files': set()})

    for file_path in file_paths:
        file_name = file_path.split("/")[-1]
        df = pd.read_csv(file_path)

        group_column = df.columns[0]
        genes_start_index = 14  # 15th column (0-based index)

        for _, row in df.iterrows():
            group_id = row[group_column]
            genes = row.iloc[genes_start_index:].dropna().astype(str).tolist()

            if len(genes) == 74:
                # Sort genes to make the gene set order-independent
                genes_sorted = tuple(sorted(genes))

                # Add entry using (Group_ID, genes_sorted) as the unique key
                key = (group_id, genes_sorted)
                group_data[key]['source_files'].add(file_name)

    print(f"Total unique (Group_ID + gene set) combinations with 74 genes: {len(group_data)}")

    # Convert to DataFrame for output
    output_rows = []
    for (group_id, genes_sorted), data in group_data.items():
        source_files = ";".join(sorted(data['source_files']))
        output_rows.append({
            'Group_ID': group_id,
            'Num_Genes': len(genes_sorted),
            'Source_Files': source_files,
            'Genes': ";".join(genes_sorted)
        })

    output_df = pd.DataFrame(output_rows)
    output_df = output_df.sort_values(by='Group_ID')

    # Save to CSV
    output_df.to_csv("distinct_groups_with_74_genes.csv", index=False)
    print("Saved output to distinct_groups_with_74_genes.csv")

# Example usage: Load all CSV files in a directory
file_paths = glob.glob("../test_data/CSVs/*.csv")
collect_74_groups_distinct(file_paths)
