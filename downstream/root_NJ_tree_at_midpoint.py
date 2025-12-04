#! /usr/bin/env python3

from ete3 import Tree
from sys import argv

tree = Tree(argv[1], format=1)

# Set midpoint root
tree.set_outgroup(tree.get_midpoint_outgroup())

# Write rooted tree
tree.write(outfile=argv[2], format=1)

