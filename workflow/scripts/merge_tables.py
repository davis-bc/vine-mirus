"""
merge_tables.py — concatenate per-chunk tool tables into one table, keeping a
single header. Used by the *_merge rules so scatter width never changes any
downstream path.

Chunk files may legitimately be header-only (a tool called nothing in that
chunk) or completely empty (more chunks than contigs).
"""

import sys

inputs = list(snakemake.input.chunks)
out = snakemake.output[0]

with open(snakemake.log[0], "w") as log, open(out, "w") as fh:
    header = None
    n_rows = 0
    for path in inputs:
        with open(path) as src:
            lines = src.read().splitlines()
        if not lines:
            print(f"{path}: empty", file=log)
            continue
        if header is None:
            header = lines[0]
            fh.write(header + "\n")
        elif lines[0] != header:
            sys.exit(f"Header mismatch in {path}:\n  {lines[0]!r}\n  expected {header!r}")
        for line in lines[1:]:
            fh.write(line + "\n")
        n_rows += len(lines) - 1
        print(f"{path}: {len(lines) - 1} rows", file=log)
    if header is None:
        print("WARNING: all chunk files were empty", file=log)
    print(f"Merged {len(inputs)} chunk(s), {n_rows} rows -> {out}", file=log)
