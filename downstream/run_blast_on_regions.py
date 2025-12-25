#! /usr/bin/env python3

from subprocess import run
import os
import multiprocessing as mp
from os.path import isfile
import argparse

def tblastn(subject,query,outfile,evalue,word_size):
    cmd = 'tblastn -subject {} -query {} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq sframe" -evalue {} -word_size {} -out {}'.format(subject,query,evalue,word_size,outfile)
    run(cmd,shell=True)
    return(0)

def blastn(subject,query,outfile,evalue,word_size):
    cmd = 'blastn -subject {} -query {} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -evalue {} -word_size {} -out {}'.format(subject,query,evalue,word_size,outfile)
    run(cmd,shell=True)
    return(0)

def clasp(infile,outfile,clasp_path,l,e):
    cmd = '{} -m -f -O -i {} -c 7 8 9 10 12 -C 1 2 -l {} -e {} -o {}'.format(clasp_path,infile,l,e,outfile).split()
    run(cmd)
    return(0)

def exec_blast_clasp(org, query, blast_type, evalue, word_size, fastas_dir, output_dir, use_clasp, clasp_path, clasp_l, clasp_e):

    if blast_type == 'tblastn':
        tblastn(f'{fastas_dir}/{org}/all.fasta',query,f'{output_dir}/blast_out_both/{org}_regions',evalue,word_size)
    else:
        blastn(f'{fastas_dir}/{org}/all.fasta',query,f'{output_dir}/blast_out_both/{org}_regions',evalue,word_size)

    # separate forward/reverse
    frame_info = {}
    expected_blast_cols = 14 if blast_type == 'tblastn' else 13
    with open(f'{output_dir}/blast_out_both/{org}_regions','r') as f:
        newlines_forward = []
        newlines_reverse = []
        for line in f:
            if not line.strip():
                continue
            line_split = line.split()
            if line_split[0] == '#':
                # BLAST comment lines (usually empty for tabular format)
                continue
            cols = line.strip().split('\t')

            # Strict column count validation
            if len(cols) != expected_blast_cols:
                raise ValueError(f"BLAST {blast_type} output has {len(cols)} columns, expected {expected_blast_cols}. Line: {line.strip()}")

            qseqid,sseqid = cols[0],cols[1]
            sstart,send = int(cols[8]),int(cols[9])

            # get frame for tblastn (column 13, the last column)
            frame = None
            if blast_type == 'tblastn':
                frame = int(cols[13])

            if sstart > send:
                cols[8],cols[9] = cols[9],cols[8]
                newlines_reverse.append('\t'.join(cols))
                if frame is not None:
                    frame_info[(qseqid,sseqid,send,sstart)] = frame
            else:
                newlines_forward.append('\t'.join(cols))
                if frame is not None:
                    frame_info[(qseqid,sseqid,sstart,send)] = frame

    with open(f'{output_dir}/blast_out_forward/{org}_regions','w') as f:
        f.write('\n'.join(newlines_forward))
    with open(f'{output_dir}/blast_out_reverse/{org}_regions','w') as f:
        f.write('\n'.join(newlines_reverse))

    if use_clasp:
        # run clasp on both orientations
        for orientation in ['forward','reverse']:
            if isfile(f'{output_dir}/blast_out_{orientation}/{org}_regions'):
                if os.path.getsize(f'{output_dir}/blast_out_{orientation}/{org}_regions') > 0:
                    clasp(f'{output_dir}/blast_out_{orientation}/{org}_regions',f'{output_dir}/clasp_out_{orientation}/{org}_regions',clasp_path,clasp_l,clasp_e)

        # write temp_coords from clasp output
        hit_no = 1
        with open(f'{output_dir}/temp_coords_{org}','w') as out:
            for orientation in ['forward','reverse']:
                clasp_file = f'{output_dir}/clasp_out_{orientation}/{org}_regions'
                if not isfile(clasp_file):
                    continue

                with open(clasp_file) as f:
                    lines = f.readlines()

                # parse C lines and collect their F lines
                i = 0
                while i < len(lines):
                    line = lines[i]
                    if not line.strip():
                        i += 1
                        continue
                    if line[0] == '#':
                        # Expected CLASP header/comment lines
                        i += 1
                        continue

                    cols = line.strip().split('\t')

                    # Only C and F lines are expected in CLASP output body
                    if cols[0] not in ['C', 'F']:
                        raise ValueError(f"Unexpected CLASP line type '{cols[0]}', expected 'C' or 'F'. Line: {line.strip()}")

                    # look for C lines (chains)
                    if cols[0] == 'C':
                        # C-lines must have exactly 8 columns
                        if len(cols) != 8:
                            raise ValueError(f"CLASP C-line has {len(cols)} columns, expected 8. Line: {line.strip()}")
                        orig_prot = cols[1]
                        id_region = cols[2]

                        # Parse BLAST subject ID format: {org}_{chr}_{start}_{end}_{orientation}
                        # Robust parsing: strip known org prefix, parse backwards from end

                        # Verify subject ID starts with expected org name
                        if not id_region.startswith(org + '_'):
                            raise ValueError(f"BLAST subject ID doesn't start with org '{org}': {id_region}")

                        # Remove org prefix
                        after_org = id_region[len(org)+1:]  # Strip "org_"

                        # Now after_org is: {chr}_{start}_{end}_{orientation}
                        # Parse from right: last 3 parts are start, end, orientation
                        parts = after_org.split('_')

                        try:
                            shift = int(parts[-3])  # start position
                            end_pos = int(parts[-2])  # end position (for validation)
                            orientation_suffix = parts[-1]  # should be 'forward' or 'reverse'
                        except (ValueError, IndexError):
                            raise ValueError(f"Cannot parse positions from BLAST subject ID: {id_region}")

                        # Chromosome is everything before the last 3 parts
                        chromo = '_'.join(parts[:-3])

                        sstart = int(cols[5])
                        send = int(cols[6])
                        start = sstart + shift
                        end = send + shift

                        gene_id = f'hit-{hit_no}_{org}_{orig_prot}'

                        if blast_type == 'tblastn':
                            # For tblastn: concatenate protein fragment sequences
                            fragment_seqs = []
                            j = i + 1
                            while j < len(lines):
                                fline = lines[j]
                                if not fline.strip():
                                    j += 1
                                    continue
                                if fline[0] == '#':
                                    # Expected CLASP header/comment lines
                                    j += 1
                                    continue
                                fcols = fline.strip().split('\t')

                                # Only C and F lines expected
                                if fcols[0] not in ['C', 'F']:
                                    raise ValueError(f"Unexpected CLASP line type '{fcols[0]}', expected 'C' or 'F'. Line: {fline.strip()}")

                                if fcols[0] == 'C':
                                    # next chain, stop
                                    break
                                if fcols[0] == 'F':
                                    # F-lines must have exactly 15 columns
                                    if len(fcols) != 15:
                                        raise ValueError(f"CLASP F-line has {len(fcols)} columns, expected 15. Line: {fline.strip()}")
                                    # fragment line - extract protein sequence from column 13 (not 12 which is bitscore)
                                    fragment_seqs.append(fcols[13])
                                j += 1
                            # concatenate protein fragments
                            full_seq = ''.join(fragment_seqs)
                            # write with protein sequence
                            out.write(f'{org}\t{chromo}\t{start}\t{end}\t{orientation}\t{gene_id}\t{full_seq}\n')
                        else:
                            # For blastn: write full genomic range without sequence
                            # blast_results_to_synthology.py will extract full range from genome
                            out.write(f'{org}\t{chromo}\t{start}\t{end}\t{orientation}\t{gene_id}\t\n')

                        hit_no += 1

                    i += 1

    return(0)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run BLAST on extracted syntenic regions',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('query', type=str, help="query protein FASTA file")
    parser.add_argument('blast_type', type=str, choices=['tblastn','blastn'], help="BLAST type (tblastn or blastn)")
    parser.add_argument('--use_clasp', type=str, nargs='?', default='y', choices=['y','n'], help="use CLASP for chaining (y or n)")
    parser.add_argument('--evalue', type=float, nargs='?', default=1e-3, help="e-value threshold")
    parser.add_argument('--word_size', type=int, nargs='?', help="word size (default: 3 for tblastn, 7 for blastn)")
    parser.add_argument('--fastas_dir', type=str, nargs='?', default='fastas', help="directory with extracted FASTAs")
    parser.add_argument('--output_dir', type=str, nargs='?', default='blast_results', help="output directory")
    parser.add_argument('--clasp_path', type=str, nargs='?', default='../clasp.x', help="path to clasp.x")
    parser.add_argument('--clasp_l', type=float, nargs='?', default=0.5, help="CLASP -l parameter")
    parser.add_argument('--clasp_e', type=float, nargs='?', default=0, help="CLASP -e parameter")
    parser.add_argument('--cores', type=int, nargs='?', default=1, help="number of cores")

    args = parser.parse_args()

    query = args.query
    blast_type = args.blast_type
    use_clasp = args.use_clasp == 'y'
    evalue = args.evalue
    word_size = args.word_size if args.word_size else (3 if blast_type == 'tblastn' else 7)
    fastas_dir = args.fastas_dir
    output_dir = args.output_dir
    clasp_path = args.clasp_path
    clasp_l = args.clasp_l
    clasp_e = args.clasp_e
    cores = args.cores

    path = os.getcwd()

    run(f'mkdir -p {output_dir}/blast_out_both',shell=True)
    if use_clasp:
        run(f'mkdir -p {output_dir}/blast_out_forward',shell=True)
        run(f'mkdir -p {output_dir}/blast_out_reverse',shell=True)
        run(f'mkdir -p {output_dir}/clasp_out_forward',shell=True)
        run(f'mkdir -p {output_dir}/clasp_out_reverse',shell=True)

    orgs = []
    if os.path.isfile('orgs'):
        with open('orgs') as f:
            for line in f:
                org = line.strip()
                if isfile(path + f'/{fastas_dir}/{org}/all.fasta'):
                    orgs.append(org)
    else:
        # fallback: list all subdirectories in fastas_dir that contain all.fasta
        if os.path.isdir(fastas_dir):
            for org in os.listdir(fastas_dir):
                if os.path.isdir(f'{fastas_dir}/{org}') and isfile(f'{fastas_dir}/{org}/all.fasta'):
                    orgs.append(org)

    print(f'{len(orgs)} organisms with FASTA files')
    if len(orgs) > 0:
        print(f'Organisms: {orgs}')

    if cores > 1:
        with mp.Pool(processes=cores) as pool:
            res = pool.starmap(exec_blast_clasp, [(org, query, blast_type, evalue, word_size, fastas_dir, output_dir, use_clasp, clasp_path, clasp_l, clasp_e) for org in orgs])
    else:
        for org in orgs:
            exec_blast_clasp(org, query, blast_type, evalue, word_size, fastas_dir, output_dir, use_clasp, clasp_path, clasp_l, clasp_e)

    print('Done')
