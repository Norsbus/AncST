"""
Stage 9: Downstream processing
Converts anchor alignments to various output formats (GFF, MCScanX, phylogeny, etc.)

Entry point: run_local.py ../anchors/
Depends on: Stage 7.5 (finalize_aligned) completing for ALL organisms

Parallelization:
- aligned_succinct files: created by finalize_aligned rule in Stage 7.5
- make_pw_table, make_gff3: parallel after mapping
- get_gff_pairwise, phylogeny: parallel after all aligned_succinct files ready
"""

import os

# Downstream directory
DOWNSTREAM_DIR = f"{ROOT_DIR}/downstream"

rule downstream_setup:
    """
    Initialize downstream directories and copy orgs file.
    Equivalent to: make_directories.py + orgs/ref_org setup

    REQUIRES: Both verbose aligned checks AND succinct checks to pass.
    If any check fails, downstream processing is blocked.
    """
    input:
        # Wait for all verbose aligned validation to complete
        [ancient(f) for f in expand(f"{WORK_DIR}/touch/checks_done_{{organism}}", organism=ALL_ORGS)],
        [ancient(f) for f in expand(f"{WORK_DIR}/touch/statistics_done_{{organism}}", organism=ALL_ORGS)],
        # Wait for all succinct validation to complete
        [ancient(f) for f in expand(f"{WORK_DIR}/touch/checks_succinct_done_{{organism}}", organism=ALL_ORGS)],
        [ancient(f) for f in expand(f"{WORK_DIR}/touch/statistics_succinct_done_{{organism}}", organism=ALL_ORGS)],
        # Need orgs file from work dir
        orgs=f"{WORK_DIR}/orgs"
    output:
        touch_file=f"{DOWNSTREAM_DIR}/touch/setup_done",
        orgs_copy=f"{DOWNSTREAM_DIR}/orgs",
        ref_org=f"{DOWNSTREAM_DIR}/ref_org"
    params:
        downstream_dir=DOWNSTREAM_DIR
    shell:
        """
        cd {params.downstream_dir}

        # Copy orgs file from template
        cp {input.orgs} orgs

        # Create ref_org (last org in list)
        tail -1 orgs > ref_org

        # Run make_directories.py to create output structure
        {PYTHON} make_directories.py

        touch {output.touch_file}
        """

rule get_all_orgs:
    """
    Extract all organisms that have matches in aligned data.
    Creates orgs_all_matches file.

    Reads: orgs (downstream copy), aligned/{org}
    """
    input:
        orgs=f"{DOWNSTREAM_DIR}/orgs",
        # ancient() prevents triggering update_candidates when anchors already exist
        aligned=ancient(expand(f"{ROOT_DIR}/anchors/aligned/{{organism}}", organism=ALL_ORGS))
    output:
        orgs_all=f"{DOWNSTREAM_DIR}/orgs_all_matches"
    params:
        downstream_dir=DOWNSTREAM_DIR,
        anchor_dir=f"{ROOT_DIR}/anchors"
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} get_all_orgs.py {params.anchor_dir}
        """

rule make_mapping:
    """
    Create chromosome and organism mapping for GFF format.
    Outputs mapping and mapping_pickle files.

    Reads: orgs (downstream copy), aligned/{org}
    Note: Does NOT actually need orgs_all_matches, but we wait for it as a sequencing point.
    """
    input:
        orgs=f"{DOWNSTREAM_DIR}/orgs",
        orgs_all=f"{DOWNSTREAM_DIR}/orgs_all_matches",  # Sequencing dependency
        # ancient() prevents triggering update_candidates when anchors already exist
        aligned=ancient(expand(f"{ROOT_DIR}/anchors/aligned/{{organism}}", organism=ALL_ORGS))
    output:
        mapping=f"{DOWNSTREAM_DIR}/mapping",
        mapping_pickle=f"{DOWNSTREAM_DIR}/mapping_pickle"
    params:
        downstream_dir=DOWNSTREAM_DIR,
        anchor_dir=f"{ROOT_DIR}/anchors"
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} make_chr_and_org_mapping_for_gff.py {params.anchor_dir}
        """

# Parallel outputs (can run after mapping). compress_maps has been replaced
# by the finalize_aligned rule in 07_update_aligned.smk which creates
# aligned_succinct/{org} directly in the anchors directory.

rule make_pw_table:
    """
    Create pairwise alignments table from aligned data.
    Can run in parallel - only needs mapping and aligned/candidates.

    Reads: orgs, mapping_pickle (via get_mapping), candidates/{org}, aligned/{org}
    """
    input:
        orgs=f"{DOWNSTREAM_DIR}/orgs",
        mapping=f"{DOWNSTREAM_DIR}/mapping_pickle",
        # ancient() prevents Snakemake from triggering update_candidates when anchors already exist
        aligned=ancient(expand(f"{ROOT_DIR}/anchors/aligned/{{organism}}", organism=ALL_ORGS)),
        candidates=ancient(expand(f"{ROOT_DIR}/anchors/candidates/{{organism}}", organism=ALL_ORGS))
    output:
        table=f"{DOWNSTREAM_DIR}/pairwise_alignments_table"
    params:
        downstream_dir=DOWNSTREAM_DIR,
        anchor_dir=f"{ROOT_DIR}/anchors"
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} make_pw_table_from_aligned.py {params.anchor_dir}
        """

rule make_gff3:
    """
    Create GFF3 format output files.
    Can run in parallel - only needs orgs file and aligned data.

    Reads: orgs, aligned/{org}
    Writes: gff/all_anchors.gff3, gff/all_alignments.gff3, gff/anchors_{org}.gff3, gff/alignments_{org}.gff3
    """
    input:
        orgs=f"{DOWNSTREAM_DIR}/orgs",
        # ancient() prevents triggering update_candidates when anchors already exist
        aligned=ancient(expand(f"{ROOT_DIR}/anchors/aligned/{{organism}}", organism=ALL_ORGS))
    output:
        all_anchors=f"{DOWNSTREAM_DIR}/gff/all_anchors.gff3",
        all_alignments=f"{DOWNSTREAM_DIR}/gff/all_alignments.gff3"
    params:
        downstream_dir=DOWNSTREAM_DIR,
        anchor_dir=f"{ROOT_DIR}/anchors"
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} make_gff3.py {params.anchor_dir}
        """

# Outputs needing all aligned_succinct files.

rule get_gff_pairwise:
    """
    Create MCScanX, DIALIGN, Smore, and i-ADHoRe format files.
    Needs ALL aligned_succinct files to be complete.

    Reads: ref_org, mapping_pickle, orgs, aligned_succinct/{org}
    Writes: MCScanX.gff, MCScanX.homology, dialign.homology, smore_anchors/, adhore_*
    """
    input:
        aligned_succinct=expand(f"{ROOT_DIR}/anchors/aligned_succinct/{{organism}}", organism=ALL_ORGS),
        orgs=f"{DOWNSTREAM_DIR}/orgs",
        ref_org=f"{DOWNSTREAM_DIR}/ref_org",
        mapping=f"{DOWNSTREAM_DIR}/mapping_pickle"
    output:
        mcscanx_gff=f"{DOWNSTREAM_DIR}/MCScanX.gff",
        mcscanx_homology=f"{DOWNSTREAM_DIR}/MCScanX.homology",
        dialign_homology=f"{DOWNSTREAM_DIR}/dialign.homology",
        adhore_config=f"{DOWNSTREAM_DIR}/adhore_config",
        adhore_table=f"{DOWNSTREAM_DIR}/adhore_pairwise.table"
    params:
        downstream_dir=DOWNSTREAM_DIR
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} get_gff_pairwise_succinct.py
        """

rule phylogeny_from_anchors:
    """
    Build UPGMA and NJ phylogenetic trees from anchor data.
    Needs ALL aligned_succinct files to be complete.

    Reads: orgs, aligned_succinct/{org}, small_meta/{org}
    """
    input:
        aligned_succinct=expand(f"{ROOT_DIR}/anchors/aligned_succinct/{{organism}}", organism=ALL_ORGS),
        orgs=f"{DOWNSTREAM_DIR}/orgs",
        small_meta=expand(f"{ROOT_DIR}/utils/small_meta/{{organism}}", organism=ALL_ORGS)
    output:
        upgma=f"{DOWNSTREAM_DIR}/UPGMATree.nwk",
        nj=f"{DOWNSTREAM_DIR}/NJTree.nwk"
    params:
        downstream_dir=DOWNSTREAM_DIR
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} phylogeny_from_anchors.py
        """

rule replace_tree_names:
    """
    Replace organism IDs with actual names in tree files.
    """
    input:
        upgma=f"{DOWNSTREAM_DIR}/UPGMATree.nwk",
        nj=f"{DOWNSTREAM_DIR}/NJTree.nwk",
        mapping=f"{DOWNSTREAM_DIR}/mapping_pickle"
    output:
        touch_file=f"{DOWNSTREAM_DIR}/touch/trees_renamed"
    params:
        downstream_dir=DOWNSTREAM_DIR
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} replace_names.py UPGMATree.nwk UPGMATree.nwk
        {PYTHON} replace_names.py NJTree.nwk NJTree.nwk
        touch {output.touch_file}
        """

rule package_output:
    """
    Package all outputs into final distribution structure.
    This is the final downstream rule.
    """
    input:
        # All required outputs
        trees_done=f"{DOWNSTREAM_DIR}/touch/trees_renamed",
        gff_pairwise=f"{DOWNSTREAM_DIR}/MCScanX.gff",
        pw_table=f"{DOWNSTREAM_DIR}/pairwise_alignments_table",
        gff3=f"{DOWNSTREAM_DIR}/gff/all_anchors.gff3"
    output:
        touch_file=f"{DOWNSTREAM_DIR}/touch/package_done",
        out_dir=directory(f"{DOWNSTREAM_DIR}/out")
    params:
        downstream_dir=DOWNSTREAM_DIR,
        anchor_dir=f"{ROOT_DIR}/anchors"
    shell:
        """
        cd {params.downstream_dir}
        {PYTHON} cp_multis_to_one.py {params.anchor_dir}
        touch {output.touch_file}
        """

rule create_full_archive:
    """
    Create full.tar.gz archive with everything in out/.
    """
    input:
        package_done=f"{DOWNSTREAM_DIR}/touch/package_done"
    output:
        archive=f"{DOWNSTREAM_DIR}/full.tar.gz"
    params:
        downstream_dir=DOWNSTREAM_DIR
    shell:
        """
        cd {params.downstream_dir}
        tar -czvf full.tar.gz out/
        """

rule create_main_archive:
    """
    Create main.tar.gz archive without gff/ and experimental/ folders.
    Contains: stable/, utils/, scripts_which_produced_these_results/, scaffolder/, README, requirements.txt
    """
    input:
        package_done=f"{DOWNSTREAM_DIR}/touch/package_done"
    output:
        archive=f"{DOWNSTREAM_DIR}/main.tar.gz"
    params:
        downstream_dir=DOWNSTREAM_DIR
    shell:
        """
        cd {params.downstream_dir}
        tar -czvf main.tar.gz --exclude='out/gff' --exclude='out/experimental' out/
        """

rule create_archives:
    """
    Checkpoint: both archives created.
    """
    input:
        full=f"{DOWNSTREAM_DIR}/full.tar.gz",
        main=f"{DOWNSTREAM_DIR}/main.tar.gz"
    output:
        touch(f"{DOWNSTREAM_DIR}/touch/archives_done")

rule downstream_all:
    """
    Complete downstream processing target (includes archives).
    """
    input:
        f"{DOWNSTREAM_DIR}/touch/archives_done"
