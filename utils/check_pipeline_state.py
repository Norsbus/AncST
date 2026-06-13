#!/usr/bin/env python3
"""
check_pipeline_state.py -- read-only structural integrity diagnostic for a
AncST pipeline working tree.

WHEN TO RUN:
    * snakemake aborts with IncompleteFilesException -- before deciding what to
      clean up, get a per-file health report
    * after an interrupted run, before --continue-run
    * any time you suspect data corruption

WHAT IT DOES:
    Walks the project's known output trees and validates each file structurally:
        anchors/aligned/{org}            - pickle stream (pickletools.genops)
        anchors/aligned_succinct/{org}   - pickle stream
        anchors/candidates/{org}         - pickle stream
        template/touch/<rule>_done_<org> - existence (a touch file has no content)
        template/parse_bcamm/{o1}/{...}  - FASTA header check
        template/sequences_to_compare/{o}/{forward,reverse}{,_split/*.fasta}
                                         - FASTA header check
        utils/genmap_out/{org}/*.freq16  - uint16 stride parity (even byte count)
        template/macle_out/{org}*        - first-line sanity for macle output

    Cross-references all per-org files against template/orgs.

WHAT IT DOES NOT DO:
    * No deletes, no moves, no writes anywhere.
    * Does not touch .snakemake/ at all (we don't depend on snakemake internals).
    * Does not reach out to a network.

OUTPUT:
    Per-file lines (tab-separated, suitable for grep / awk):
        <STATUS>\\t<size>\\t<mtime_iso>\\t<path>\\t<detail>
    where STATUS in {PASS, FAIL, MISSING, UNKNOWN}.

    Then a summary table by category, plus per-org coverage.

EXIT CODES:
    0 - no FAIL anywhere (MISSING / UNKNOWN are informational only)
    1 - at least one FAIL detected (data corruption, manual decision needed)
    2 - usage / invocation error

PORTABILITY:
    Python 3.8+. macOS + Linux, x86_64 + arm64.

    Uses concurrent.futures.ThreadPoolExecutor for parallelism, not
    multiprocessing -- the work is I/O-bound (reading bytes from disk and
    streaming through pickletools) so threads avoid the spawn-vs-fork gotcha
    on macOS (default start method changed to 'spawn' in 3.8).

NEVER DELETES ANYTHING. Per the project's never-remove policy, this
script reports problems and points the user at backups; it does not act on
its own findings.
"""

from __future__ import annotations

import argparse
import os
import pickletools
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Iterator, NamedTuple, Tuple


# ---------------------------------------------------------------------------
# Per-file structural checks. Each returns (status, detail). Pure functions.
# ---------------------------------------------------------------------------

class Report(NamedTuple):
    status: str    # PASS / FAIL / MISSING / UNKNOWN
    size: int      # bytes (0 if MISSING)
    mtime: str     # iso-formatted (- if MISSING)
    path: str
    category: str  # e.g. "aligned", "candidates", "touch", "freq16", "fasta_split"
    detail: str    # human-readable note


def _stat_or_none(path: Path):
    try:
        return path.stat()
    except OSError:
        return None


def _check_pickle_stream(path: Path) -> Tuple[str, str]:
    """Stream the pickle ops without materializing the object graph.

    pickletools.genops yields (opcode, arg, pos) for each op. The pickle
    stream ends with a STOP opcode (b'.') -- if we don't see one, the file
    is truncated.  Materializing the object (pickle.load) would cost RAM
    proportional to the data; this costs only a per-op tuple, irrespective
    of payload size, so it's safe to run on multi-GB aligned/{org} files.
    """
    try:
        with path.open('rb') as fh:
            saw_stop = False
            for op, arg, pos in pickletools.genops(fh):
                if op.name == 'STOP':
                    saw_stop = True
            if not saw_stop:
                return 'FAIL', 'pickle stream has no STOP opcode (truncated)'
            return 'PASS', 'pickle stream walked OK (genops)'
    except ValueError as e:
        # pickletools raises ValueError on malformed opcodes
        return 'FAIL', f'pickletools.ValueError: {e}'
    except EOFError as e:
        return 'FAIL', f'EOF mid-pickle: {e}'
    except OSError as e:
        return 'FAIL', f'OSError reading pickle: {e}'
    except Exception as e:
        return 'FAIL', f'{type(e).__name__}: {e}'


def _check_fasta_header(path: Path) -> Tuple[str, str]:
    """A FASTA file must start with a '>' character on its first line."""
    try:
        with path.open('rb') as fh:
            first = fh.readline()
        if not first:
            return 'FAIL', 'empty file'
        if not first.startswith(b'>'):
            return 'FAIL', "first line does not start with '>'"
        return 'PASS', 'FASTA header present'
    except OSError as e:
        return 'FAIL', f'OSError: {e}'


def _check_freq16(path: Path) -> Tuple[str, str]:
    """genmap *.freq16 files are uint16 streams -- size must be even.

    Cheaper than parsing the full array; catches truncation at odd byte counts.
    """
    try:
        sz = path.stat().st_size
    except OSError as e:
        return 'FAIL', f'OSError: {e}'
    if sz == 0:
        return 'FAIL', 'zero-length freq16'
    if sz % 2 != 0:
        return 'FAIL', f'odd byte count {sz} (not uint16 stride)'
    return 'PASS', f'even byte count ({sz} bytes)'


def _check_existence(path: Path) -> Tuple[str, str]:
    """Touch files: presence is all that matters (content empty by design)."""
    return 'PASS', 'present'


def _check_macle_txt(path: Path) -> Tuple[str, str]:
    """macle output: tab-separated, exactly 3 fields per non-blank line
    (chrom, position-as-int, complexity-as-float). Only sample the first
    N lines and the last line -- we want truncation detection, not a full
    re-parse of multi-million-line files.
    """
    try:
        with path.open('rb') as fh:
            sample = fh.readlines(8192)
        if not sample:
            return 'FAIL', 'empty macle output'
        for i, raw in enumerate(sample):
            if not raw.strip():
                continue
            parts = raw.rstrip(b'\n').split(b'\t')
            if len(parts) != 3:
                return 'FAIL', f'line {i+1}: expected 3 tab-separated fields, got {len(parts)}'
            try:
                int(parts[1])
                float(parts[2])
            except ValueError as e:
                return 'FAIL', f'line {i+1}: numeric parse failed: {e}'
        # last-line truncation check: if file is bigger than 8KB, look at tail
        sz = path.stat().st_size
        if sz > 8192:
            with path.open('rb') as fh:
                fh.seek(max(0, sz - 4096))
                tail_bytes = fh.read()
            # Final newline should exist; last non-empty line should still parse.
            tail_lines = [ln for ln in tail_bytes.split(b'\n') if ln.strip()]
            if tail_lines:
                parts = tail_lines[-1].split(b'\t')
                if len(parts) != 3:
                    return 'FAIL', f'tail: expected 3 fields, got {len(parts)} (truncated mid-line?)'
                try:
                    int(parts[1])
                    float(parts[2])
                except ValueError as e:
                    return 'FAIL', f'tail line parse failed: {e}'
        return 'PASS', f'macle .txt looks well-formed ({sz} bytes)'
    except OSError as e:
        return 'FAIL', f'OSError: {e}'


# ---------------------------------------------------------------------------
# Output-tree discovery. Each generator yields (category, path) tuples.
# Pure file-system walks; no validation. The validator dispatches by category.
# ---------------------------------------------------------------------------

def iter_anchor_pickles(root: Path) -> Iterator[Tuple[str, Path]]:
    for sub in ('aligned', 'aligned_succinct', 'candidates'):
        d = root / 'anchors' / sub
        if d.is_dir():
            for f in sorted(d.iterdir()):
                if f.is_file():
                    yield sub, f


def iter_touch_files(work_dir: Path) -> Iterator[Tuple[str, Path]]:
    d = work_dir / 'touch'
    if d.is_dir():
        for f in sorted(d.iterdir()):
            if f.is_file():
                yield 'touch', f


def iter_freq16(root: Path) -> Iterator[Tuple[str, Path]]:
    # Tested layouts: utils/genmap_out/{org}/{k}_{e}.freq16
    base = root / 'utils' / 'genmap_out'
    if base.is_dir():
        for org_dir in sorted(base.iterdir()):
            if org_dir.is_dir():
                for f in sorted(org_dir.glob('*.freq16')):
                    yield 'freq16', f


def iter_parse_bcamm(work_dir: Path) -> Iterator[Tuple[str, Path]]:
    # NOTE: parse_bcamm output files are named *.fasta by historical accident
    # but their content is PICKLE (first byte 0x80 = pickle protocol marker),
    # not FASTA. The category name reflects content, not file extension.
    base = work_dir / 'parse_bcamm'
    if base.is_dir():
        for org_dir in sorted(base.iterdir()):
            if org_dir.is_dir():
                for f in sorted(org_dir.glob('*.fasta')):
                    yield 'parse_bcamm_pickle', f


def iter_sequences_to_compare(work_dir: Path) -> Iterator[Tuple[str, Path]]:
    base = work_dir / 'sequences_to_compare'
    if base.is_dir():
        # Top-level forward.fasta / reverse.fasta per org dir, plus split subdirs.
        for org_dir in sorted(base.iterdir()):
            if not org_dir.is_dir():
                continue
            for f in sorted(org_dir.glob('*.fasta')):
                if f.is_file():
                    yield 'seqcmp_fasta', f
            for split_sub in ('forward_split', 'reverse_split'):
                sd = org_dir / split_sub
                if sd.is_dir():
                    for f in sorted(sd.glob('*.fasta')):
                        yield 'seqcmp_split', f


def iter_macle_out(root: Path) -> Iterator[Tuple[str, Path]]:
    # macle output lives at utils/macle_out/{org}/{window}_{percentile}.txt
    # Tab-separated: chrom \t position \t complexity.
    base = root / 'utils' / 'macle_out'
    if base.is_dir():
        for org_dir in sorted(base.iterdir()):
            if org_dir.is_dir():
                for f in sorted(org_dir.glob('*.txt')):
                    if f.is_file():
                        yield 'macle_out', f


# Validator dispatch by category -- thread-safe (no shared state).
# Caveat: parse_bcamm files are named *.fasta but contain pickle data.
_VALIDATORS = {
    'aligned':            _check_pickle_stream,
    'aligned_succinct':   _check_pickle_stream,
    'candidates':         _check_pickle_stream,
    'touch':              _check_existence,
    'freq16':             _check_freq16,
    'parse_bcamm_pickle': _check_pickle_stream,  # .fasta extension, pickle content
    'seqcmp_fasta':       _check_fasta_header,
    'seqcmp_split':       _check_fasta_header,
    'macle_out':          _check_macle_txt,
}


def validate_one(category: str, path: Path) -> Report:
    st = _stat_or_none(path)
    if st is None:
        return Report('MISSING', 0, '-', str(path), category, 'file does not exist')
    size = st.st_size
    mtime = time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(st.st_mtime))
    fn = _VALIDATORS.get(category)
    if fn is None:
        return Report('UNKNOWN', size, mtime, str(path), category,
                      'no structural check for this category')
    status, detail = fn(path)
    return Report(status, size, mtime, str(path), category, detail)


# ---------------------------------------------------------------------------
# Cross-org coverage check: every org in template/orgs should have a triplet
# of candidates / aligned / aligned_succinct.
# ---------------------------------------------------------------------------

def read_orgs_file(orgs_path: Path) -> list[str]:
    orgs = []
    if not orgs_path.is_file():
        return orgs
    for line in orgs_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        orgs.append(line)
    return orgs


def per_org_coverage(root: Path, orgs: Iterable[str]) -> list[Report]:
    """For each org, MISSING-report any expected anchor file that isn't on disk."""
    reports = []
    for org in orgs:
        for sub in ('candidates', 'aligned', 'aligned_succinct'):
            expected = root / 'anchors' / sub / org
            if not expected.exists():
                reports.append(Report(
                    'MISSING', 0, '-', str(expected), sub,
                    f'expected per-org anchor file for {org} not on disk'
                ))
    return reports


# ---------------------------------------------------------------------------
# CLI + driver.
# ---------------------------------------------------------------------------

def collect_files(root: Path, work_dir: Path) -> list[Tuple[str, Path]]:
    items: list[Tuple[str, Path]] = []
    items.extend(iter_anchor_pickles(root))
    items.extend(iter_touch_files(work_dir))
    items.extend(iter_freq16(root))
    items.extend(iter_parse_bcamm(work_dir))
    items.extend(iter_sequences_to_compare(work_dir))
    items.extend(iter_macle_out(root))
    return items


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog='check_pipeline_state.py',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('root',
                   help='Project root (directory containing anchors/, template/, utils/, ...)')
    p.add_argument('--cores', type=int, default=max(1, (os.cpu_count() or 1) // 2),
                   help='Validation parallelism (ThreadPoolExecutor workers). '
                        'Default: half the visible CPUs (min 1).')
    p.add_argument('--work-dir', default=None,
                   help='Work dir holding transient files (touch/, parse_bcamm/, '
                        'sequences_to_compare/, orgs). Default: <root>/template.')
    p.add_argument('--orgs', default=None,
                   help='Override path to orgs file. Default: <work-dir>/orgs.')
    p.add_argument('--quiet', '-q', action='store_true',
                   help='Suppress per-file lines; print only the summary.')
    p.add_argument('--only-fail', action='store_true',
                   help='In per-file output, print only FAIL lines.')
    return p.parse_args(argv)


def emit_per_file(rep: Report, out=sys.stdout) -> None:
    out.write('\t'.join((rep.status, str(rep.size), rep.mtime,
                         rep.path, rep.detail)) + '\n')


def summarize(reports: list[Report]) -> dict:
    by_status: dict[str, int] = {'PASS': 0, 'FAIL': 0, 'MISSING': 0, 'UNKNOWN': 0}
    by_category: dict[str, dict[str, int]] = {}
    for r in reports:
        by_status[r.status] = by_status.get(r.status, 0) + 1
        cat = by_category.setdefault(r.category,
                                     {'PASS': 0, 'FAIL': 0, 'MISSING': 0, 'UNKNOWN': 0})
        cat[r.status] = cat.get(r.status, 0) + 1
    return {'total': by_status, 'per_category': by_category}


def print_summary(summary: dict, out=sys.stderr) -> None:
    out.write('\n')
    out.write('=' * 64 + '\n')
    out.write('check_pipeline_state.py -- summary\n')
    out.write('=' * 64 + '\n')
    t = summary['total']
    out.write(f"  total checked: {sum(t.values())}\n")
    out.write(f"    PASS:    {t['PASS']}\n")
    out.write(f"    FAIL:    {t['FAIL']}\n")
    out.write(f"    MISSING: {t['MISSING']}\n")
    out.write(f"    UNKNOWN: {t['UNKNOWN']}\n")
    out.write('\n  by category:\n')
    cats = sorted(summary['per_category'].keys())
    for c in cats:
        row = summary['per_category'][c]
        out.write(f"    {c:24s}  PASS={row['PASS']:5d} "
                  f"FAIL={row['FAIL']:5d} "
                  f"MISSING={row['MISSING']:5d} "
                  f"UNKNOWN={row['UNKNOWN']:5d}\n")
    out.write('=' * 64 + '\n')


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"error: root is not a directory: {root}", file=sys.stderr)
        return 2

    work_dir = Path(args.work_dir).resolve() if args.work_dir else (root / 'template')
    orgs_path = Path(args.orgs) if args.orgs else (work_dir / 'orgs')
    orgs = read_orgs_file(orgs_path)
    if not orgs:
        print(f"warning: no organisms read from {orgs_path}", file=sys.stderr)

    items = collect_files(root, work_dir)
    if not items and not orgs:
        print(f"warning: no output files discovered under {root}", file=sys.stderr)

    # ThreadPoolExecutor: I/O-bound work, portable across mac/linux/intel/arm,
    # no spawn-vs-fork start-method issue. Cap workers to len(items) to avoid
    # creating threads for nothing.
    cores = max(1, min(args.cores, max(1, len(items))))
    reports: list[Report] = []
    with ThreadPoolExecutor(max_workers=cores) as ex:
        futures = [ex.submit(validate_one, cat, p) for (cat, p) in items]
        for fut in as_completed(futures):
            try:
                rep = fut.result()
            except Exception as e:
                # validator itself crashed -- record so we don't lose the signal
                rep = Report('FAIL', 0, '-', '<validator-internal-error>',
                             '?', f'{type(e).__name__}: {e}')
            reports.append(rep)

    # Add cross-org coverage gaps (these are MISSING, not FAIL -- informational)
    if orgs:
        reports.extend(per_org_coverage(root, orgs))

    # Per-file output (stable sort by status then path so FAILs cluster at top
    # when piping through grep -E '^FAIL|^MISSING')
    status_order = {'FAIL': 0, 'MISSING': 1, 'UNKNOWN': 2, 'PASS': 3}
    reports.sort(key=lambda r: (status_order.get(r.status, 99), r.path))

    if not args.quiet:
        for r in reports:
            if args.only_fail and r.status != 'FAIL':
                continue
            emit_per_file(r)

    summary = summarize(reports)
    print_summary(summary)

    # Exit code: 1 only on real corruption (FAIL). MISSING is informational
    # (often expected when a stage hasn't run yet); UNKNOWN means "we don't
    # have a structural check for this file type" and is not failure.
    return 1 if summary['total']['FAIL'] > 0 else 0


if __name__ == '__main__':
    sys.exit(main())
