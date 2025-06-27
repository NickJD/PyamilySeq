#!/usr/bin/env python3

import sys
from ete3 import Tree

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} tree1.nwk tree2.nwk")
    sys.exit(1)

# Load trees
t1 = Tree(open(sys.argv[1]).read(), format=1, quoted_node_names=True)



t2 = Tree(open(sys.argv[2]).read(), format=1, quoted_node_names=True)



# Compute RF distance
rf_result = t1.robinson_foulds(t2, unrooted_trees=True)

# Handle both dict and tuple return types
if isinstance(rf_result, dict):
    rf = rf_result.get("rf", 0)
    max_rf = rf_result.get("max_rf", 1)
else:
    rf, max_rf = rf_result[0], rf_result[1]

# Print results
print(f"Raw RF distance: {rf}")
print(f"Max RF distance: {max_rf}")
if max_rf > 0:
    print(f"Normalised RF: {rf / max_rf:.5f}")
else:
    print("Normalised RF: N/A (max_rf is 0)")
