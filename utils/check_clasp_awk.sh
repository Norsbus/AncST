#!/bin/bash
# Script to check clasp_out AWK processing without permission prompts
# Replaces: diff <(grep -v "^#" clasp_file) <(cat awk_file) | head -40

if [ $# -lt 1 ]; then
    echo "Usage: $0 <clasp_file_base> [lines_to_show]"
    echo "Example: $0 clasp_out_forward/GCF_000001215.4/GCF_000001215.4.part-081.fasta"
    exit 1
fi

clasp_file="$1"
awk_file="${clasp_file}_awk"
lines="${2:-40}"

echo "Comparing files:"
echo "  Original: $clasp_file (without comment lines)"
echo "  AWK processed: $awk_file"
echo "---"

# Use the safe_diff script
./utils/safe_diff.sh "$clasp_file" "$awk_file" "$lines"

# Also show counts for verification
echo "---"
echo "Line counts:"
echo -n "  Original (no comments): "
grep -v "^#" "$clasp_file" 2>/dev/null | wc -l || echo "File not found"
echo -n "  AWK processed: "
wc -l < "$awk_file" 2>/dev/null || echo "File not found"