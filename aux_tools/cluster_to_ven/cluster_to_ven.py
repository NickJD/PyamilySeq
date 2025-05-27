#!/usr/bin/env python3
import re
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib_venn import venn3

# === CONFIGURATION ===
CLSTR_FILE = "rpp_combined_aa_100_100.clstr"    # your CD-HIT .clstr file
TOOLS = ["roary", "panaroo", "pyamilyseq_80"]


def parse_cdhit_clstr(path):
    """
    Returns a dict:
       cluster_id (int) -> set of tool names present in that cluster
    """
    cluster_map = defaultdict(set)
    current = None
    with open(path) as f:
        for line in f:
            if line.startswith(">Cluster"):
                # e.g. ">Cluster 1"
                current = int(line.split()[1])
            else:
                # e.g. "1    234aa, >panaroo|group_10157... at 100.00%"
                m = re.search(r">([^|>]+)\|", line)
                if m and current is not None:
                    tool = m.group(1)
                    cluster_map[current].add(tool)
    return cluster_map


def build_tool_sets(cluster_map, tools):
    """
    Given cluster_map and list of tools,
    returns a dict tool->set(cluster_ids).
    """
    sets = {tool: set() for tool in tools}
    for cid, present in cluster_map.items():
        for tool in tools:
            if tool in present:
                sets[tool].add(cid)
    return sets


def plot_venn3(sets_dict, tools):
    """
    Given dict {tool: set(cids)} for exactly 3 tools,
    draws and shows a venn3.
    """
    a, b, c = [sets_dict[t] for t in tools]
    plt.figure(figsize=(6,6))
    v = venn3([a, b, c], set_labels=tools)
    plt.title("Clusters present per tool")
    plt.tight_layout()
    plt.show()


def main():
    cluster_map = parse_cdhit_clstr(CLSTR_FILE)
    tool_sets   = build_tool_sets(cluster_map, TOOLS)
    plot_venn3(tool_sets, TOOLS)


if __name__ == "__main__":
    main()
