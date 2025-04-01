# AncST - Anchor Synteny Tool   

This project provides a pipeline and auxiliary code to find potential synteny anchor candidates from genomes and making pairwise alignments. Those can then be used by downstream programmes to study macro- and microsynteny.
It is only a prototype version at the moment: works in principle; almost certainly needs programming skills and patience; bioinformatics experience helps.

## Webserver/Precomputed Anchors
There is an automated webserver for anchor computation as well as precomputed datasets (anchors + utility scripts) available at [http://anchored.bioinf.uni-leipzig.de:8080/](http://anchored.bioinf.uni-leipzig.de:8080/).

## Description

This project can roughly be subdivided into 2 main features:
1. anchor candidate computation and pairwise match making
* the core pipeline takes a set of genomes and computes high confidence pairwise alignments
2. downstream programmes for macro- and microsyntenic analyses, e.g. build colinear chains and find syntelogs
* the pairwise alignments can then be fed into MCScanX or i-ADHoRe to build (transitive) colinear chains
* further custom analysis and plotting scripts are provided to draw the alignments for visual inspection and include genetic elements of interest whose potential syntenic orthologs ("syntelogs") can be found using the pairwise and transitive multiple alignments as anchors to reduce the search space


## Getting Started

### Dependencies

* the pipeline needs to be run on linux and is using snakemake
* the version of the anchor making pipeline finalized and tested the most thoroughly is a version running on slurm clusters
* there is a requirements.txt file with necessary python libraries to be used by conda - this pipeline uses was developed in an environment called snakemake. please make one called like that and if you really cannot debug the errors of where it fails because of that specific name which should not take too long

### Installing

* make sure you are running on a slurm cluster or prepare to be brave and debug the local verion (more below) for your system
* install the python packages in requirements.txt

### Executing program

1. you need to prepare some things:
* some file paths need to be added (your home, not really important, just add the parent directory of the repo; and your conda home where the envrionment you made with requirements is located (see above)). replace accordingly and run (order of commands is important):
```
grep -rlZ "CONDAHOMEDIR" . | xargs -0 sed -i 's#CONDAHOMEDIR#your\_conda\_home#g'
grep -rlZ "HOMEDIR" . | xargs -0 sed -i 's#HOMEDIR#your\_home#g'
```
* put the genomes in a file called "orgs" with each genome's name on one line (without suffix fa/fasta/whatever). this file should exist in the directories utils,utils/util\_code and template
* go to utils and execute:
```
python3 make_directories.py
```
* go to utils/util\_code and execute:
```
python3 make_anchor_directories.py
```
* either put your genomes in format "name\_as\_in\_orgs\_file.fasta" or if your are using genomes pubclicly available in NCBI you can go to utils/util\_code and execute (alternatively execute the ...\_mp.py versions if you have multiple cores available and configure in the code how many):
```
python3 get_genomes.py
python3 extract_fastas.py
```
2. now go to the template directory (make sure you have the same orgs file there):
* make sure you have the same orgs file there or one with a subset of the species (which you then only want to consider)
* make a "compute\_anchors\_for" file which is the same format as orgs and defines the subset of species for which you want to compute new anchor candidates for. thus, if you already have anchors for some of the species and you do not want to compute new ones for them, leave them out of this file but include them in the orgs file, so that their anchors will be matched against the other species.
* configure you parameters, you need: "filter\_params" and "orig\_params". "filter\_params" defines which k and e values will be used for GenMap, how big the window size and pitch are to create potential anchor candidates and which best percentile of the candidates will further be processed. thus, put these numbers separated by whitespaces per line in "filter\_params" (see the downloaded one for an example and modify). "orig\_params" includes the same 5 numbers plus two for blast word sizes and a number which is irrelevant now but better include it (again see the example and modify). the first blast word size is for searches against the own genome and the other for pairwise matching.
* then execute:
```
python3 run_complete_snakemake.py
```
3. you can know use the anchors/alignments stored in anchors/aligned to conduct downstream analyses:
* (if you have a big dataset and resources available contact me for multiprocessing or slurm versions of this but usually its fine to run small datasets on one core (lets say in total genomic size <= 10GB))
* go to downstream analyses, put the same orgs file (or a subset of the species if you want to consider only those) and a "ref\_org" file which is just one of the species in orgs that will serve as the reference orientation of sequences; then execute:
```
python3 run_local.py {anchor_directory_path==normally now just ../anchors/}
```
* this will create an "out" directory containing various possibilities and READMEs
* e.g. you can run colinearity tools pretty much out of the box and draw some synteny pictures
* you can also use the anchors and colinearity chains to cluster potential homologs into likely syntelog groups

## Help

karl.kaether@posteo.de

* this project is still in the prototype stage, thus there will likely be errorrs and bugs. most of thenm are due to specificities of different datasets and my inability to make this pipeline completely automatic/generic. thus i have likely seen most errors and can help you.
* i have additional scripts which can be helpful (converting formats,extracting syntenic regions, improving assemblies with the alignments, etc.), so dont hesitate to ask

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This project is licensed under the GNU GPL v3 License - see the LICENSE.md file for details
