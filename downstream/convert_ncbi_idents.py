#!/usr/bin/env python3
import re
import json
import subprocess
import sys
import time
import os

def get_species_name(accession):
    """Return organism name for a given accession using NCBI datasets CLI."""
    try:
        result = subprocess.run(
            [os.path.expanduser("~/datasets"), "summary", "genome", "accession", accession],
            capture_output=True, text=True, check=True
        )
        data = json.loads(result.stdout)
        return data["reports"][0]["organism"]["organism_name"]
    except Exception as e:
        print(f"⚠️  Failed to fetch {accession}: {e}", file=sys.stderr)
        return None

def main():
    input_file = "tree"
    output_file = "tree_with_species.txt"

    with open(input_file) as f:
        tree_text = f.read()

    accessions = sorted(set(re.findall(r"(GC[FA]_\d+\.\d+)", tree_text)))
    print(f"🔍 Found {len(accessions)} accessions.")

    mapping = {}
    for acc in accessions:
        print(f"→ {acc}", end=": ")
        species = get_species_name(acc)
        if species:
            formatted = f"{species.replace(' ', '_')}[{acc}]"
            mapping[acc] = formatted
            print(formatted)
        else:
            mapping[acc] = acc
            print("⚠️  (kept accession)")
        time.sleep(0.2)

    def replace_func(m):
        acc = m.group(1)
        return mapping.get(acc, acc)

    new_tree = re.sub(r"(GC[FA]_\d+\.\d+)", replace_func, tree_text)

    with open(output_file, "w") as f:
        f.write(new_tree)

    print(f"\n✅ Done! Wrote output to '{output_file}'.")
    for k, v in list(mapping.items())[:5]:
        print(f"  {k} → {v}")

if __name__ == "__main__":
    main()

