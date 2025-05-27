#!/usr/bin/env python3
import re
from collections import defaultdict
import pandas as pd
from upsetplot import UpSet, from_memberships
import matplotlib.pyplot as plt

CLSTR_FILE = "combined_aa_100_100.clstr"
TOOLS      = ["roary", "roary_dsp", "panaroo", "pyamilyseq_0", "pyamilyseq_80"]

def parse_cdhit_clstr(path):
    cluster_map = defaultdict(set)
    current = None
    with open(path) as f:
        for line in f:
            if line.startswith(">Cluster"):
                current = int(line.split()[1])
            else:
                m = re.search(r">([^|>]+)\|", line)
                if m and current is not None:
                    cluster_map[current].add(m.group(1))
    return cluster_map

def main():
    cluster_map = parse_cdhit_clstr(CLSTR_FILE)

    # build list of memberships
    memberships = []
    for cid, present in cluster_map.items():
        # only include tools that are in our defined list
        tools_here = tuple(sorted(p for p in present if p in TOOLS))
        memberships.append(tools_here)

    # build the UpSet data
    data = from_memberships(memberships)
    upset = UpSet(data, subset_size='count', show_counts=True)
    upset.plot()
    plt.title("Cluster presence/absence across tools")
    plt.show()

if __name__ == "__main__":
    main()
