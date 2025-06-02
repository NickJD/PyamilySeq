#!/bin/bash

#set -euo pipefail

ALIGN_DIR="alignments"
TREE_DIR="trees"
RF_MATRIX="rf_matrix.csv"
NSEQ_EXPECTED=74

#mkdir -p "$TREE_DIR"

echo "Step 1: Building trees with FastTree..."

for aln in "$ALIGN_DIR"/*.fa*; do
    unique_genomes=$(grep '^>' "$aln" | sed -E 's/^>([^_]+_combined).*$/\1/' | sort -u | wc -l | xargs)

    if [ "$unique_genomes" -eq "$NSEQ_EXPECTED" ]; then
        filename=$(basename "$aln")
        basename="${filename%%.*}"
        tree_file="$TREE_DIR/${basename}.tree"

        echo "Processing $basename: $unique_genomes unique genomes → building tree"

        # Clean headers: keep only everything before "_combined"
        tmpfile=$(mktemp)
        sed -E 's/_combined.*//' "$aln" > "$tmpfile"

        # Build tree
        FastTree -quiet -nt "$tmpfile" > "$tree_file"
        rm "$tmpfile"
    else
        echo "Skipping $aln: found $unique_genomes unique genomes (expected $NSEQ_EXPECTED)"
    fi
done
