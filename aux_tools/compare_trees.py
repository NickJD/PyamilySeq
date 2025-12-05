from ete3 import Tree

# Load trees from Newick format
tree1 = Tree("/mnt/Internal/Nextcloud/Current_Work/PyamilySeq/E-coli/74_genome/PyamilySeq_Full_cds_AA_90_0/core_gene_alignment.nwk")  # Replace with your file path
tree2 = Tree("/mnt/Internal/Nextcloud/Current_Work/PyamilySeq/E-coli/74_genome/PyamilySeq_Full_DNA_90_0/core_gene_alignment.nwk")  # Replace with your file path


rf_result = tree1.robinson_foulds(tree2, unrooted_trees=True)

# Unpack based on 7 returned values
rf_distance, max_rf, common_clades, parts_t1, parts_t2, extra1, extra2 = rf_result

# Compute normalised R-F distance
normalised_rf = rf_distance / max_rf if max_rf > 0 else 0

# Print results
print(f"Robinson-Foulds Distance: {rf_distance}")
print(f"Maximum Possible RF Distance: {max_rf}")
print(f"Normalised RF Distance: {normalised_rf:.6f}")
# print(f"Maximum Possible RF Distance: {max_rf}")
# print(f"Common Clades: {common_clades}")
# print(f"Unique Clades in Tree 1: {parts_t1}")
# print(f"Unique Clades in Tree 2: {parts_t2}")
# print(f"Extra Value 1: {extra1}")
# print(f"Extra Value 2: {extra2}")


