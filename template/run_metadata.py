#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
from multiprocessing import Pool,cpu_count
import os

def thread_process(org):
    if not (pathlib.Path(f'{work_dir}/../utils/small_meta/{org}').exists() and pathlib.Path(f'{work_dir}/../utils/metadata_genomes/{org}').exists() and pathlib.Path(f'{work_dir}/../utils/blastdbs/{org}.nsq').exists()):
        cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH --job-name=metadata_{org}\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --mem=64000\n\
#SBATCH --time=3-00:00:00\n\
#SBATCH --error log/metadata/stderr_metadata_{org}\n\
#SBATCH --output log/metadata/stdout_metadata_{org}\n\
    eval "$(conda shell.bash hook)"\n\
    conda activate CONDAHOMEDIR/miniconda3/envs/snakemake\n\
    ./make_blastdb_and_metadata.py {org} && touch touch/{org}_metadata_done\n\
    """
        with open('run.slurm','w') as f:
            f.write(cmd)
        run(f'sbatch run.slurm',shell=True)
    else:
        run(f'touch touch/{org}_metadata_done',shell=True)
        print(f'blastdbs and metadata exists for {org}')
    return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = 'HOMEDIR/tmp'
    orgs = set()
    with open('compute_anchors_for') as f:
        for line in f:
            orgs.add(line.strip())
    for org in orgs:
        thread_process(org)
