- PREREQUISITE: the blossom binary in blossom5-v2.05.src has to be executable. otherwise compile from https://pub.ista.ac.at/~vnk/software.html
- also some python packages may be needed but not too many and they are common and stable (numpy,networkx,...)

- code to scaffold a target genome based on one or more reference genomes
- if using a single reference just run "python3 single_scaff.py {ref} {target}"
- otherwise run: "python3 AncST_scaff.py {target}" but you need:
																	- a "refs.txt" file with the name of each ref for which anchors exist with the same name in ../anchors/candidates and ../anchors/aligned (line separated)
																	- a "ref_weights.txt" file with the ref and its weight separated by one or more whitespaces per line
																	(- each ref needs to have a weight)
- this will produce the output files for single ref runs in singles_out and the multi-ref results in multi_out
- there will be the newly assembled scaffolds with orientations of the contigs as well as files documenting contigs which show evidence of larger rearrangement events, particularly they show alignments to different ref chromosomes and stretches of differentially oriented alignment, this indicating translocations breakpoints and inversions, respectively
