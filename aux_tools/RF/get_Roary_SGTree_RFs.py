import os
import glob
import csv
from ete3 import Tree

tree_dir = "/mnt/Internal/Nextcloud/Current_Work/PyamilySeq/E-coli/74_genome/roary_default_aligned_8cpu_keepfiles/trees"


tree_files = sorted(glob.glob(os.path.join(tree_dir, "*.tree")))
tree_names = [os.path.splitext(os.path.basename(f))[0] for f in tree_files]
trees = [Tree(open(f).read(), format=1) for f in tree_files]

# Prepare matrices
raw_matrix = []
header = [""] + tree_names
raw_matrix.append(header)

rf_values = []

# Step 1: Compute RF distances and collect all values
rf_data = []
for i, t1 in enumerate(trees):
    row = [tree_names[i]]
    row_vals = []
    for j, t2 in enumerate(trees):
        if i == j:
            rf = 0
        else:
            rf, *_ = t1.robinson_foulds(t2, unrooted_trees=True)
        row_vals.append(rf)
        rf_values.append(rf)
    rf_data.append((row[0], row_vals))

# Step 2: Normalize based on max RF distance (excluding diagonal zeros)
max_rf = max(rf_values) if rf_values else 1

# Step 3: Build both raw and normalised matrices
norm_matrix = [header]
for name, row_vals in rf_data:
    print(f"Processing {name} with RF values: {row_vals}")
    raw_row = [name] + [str(rf) for rf in row_vals]
    if max_rf == 0:
        norm_row = [name] + ["0.00000" for _ in row_vals]
    else:
        norm_row = [name] + [f"{rf / max_rf:.5f}" for rf in row_vals]
    raw_matrix.append(raw_row)
    norm_matrix.append(norm_row)

# Step 4: Write to CSV
with open("rf_matrix.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(raw_matrix)

with open("pyamilyseq_80_rf_matrix_normalised_0_1.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(norm_matrix)

print("✓ Raw RF matrix written to rf_matrix.csv")
print("✓ Min-max normalized RF matrix written to xxx_rf_matrix_normalised_0_1.csv")

