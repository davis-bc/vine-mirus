#!/usr/bin/env bash
# Post-deploy hook for envs/cenote_taker3.yaml.
#
# Snakemake runs this automatically after it creates the environment, with
# $CONDA_PREFIX pointing at it. That matters: a manual binary swap gets silently
# reverted the next time Snakemake rebuilds the env (a changed yaml, a new
# --conda-prefix, a fresh clone), and the symptom is Cenote-Taker3 dying with
# "Illegal instruction" halfway through a long run.
#
# What it works around: current MMseqs2 conda builds use AVX2. On a CPU that
# only has SSE4.2 — goldberg-r710, and anything else of that generation — the
# AVX2 build dies with SIGILL / "Illegal instruction". MMseqs2 publishes a
# static SSE4.1 build, which runs there; it is slower, and pointless on a
# machine with AVX2.
#
# Opt in per machine, since the rest of the pipeline is meant to run on Slurm
# nodes that do have AVX2:
#
#     export VINE_LEGACY_CPU=1
#
# Same family of problem as the megahit=1.2.8 pin in envs/qc_assembly.yaml and
# the diamond=2.1.8 pin in envs/checkv.yaml, both of which are harmless on
# modern hardware and so are pinned unconditionally. This one is a binary
# download rather than a pin, so it stays opt-in.

set -euo pipefail

MMSEQS_VERSION="18-8cc5c"
MMSEQS_URL="https://mmseqs.com/archive/${MMSEQS_VERSION}/mmseqs-linux-sse41.tar.gz"

if [ "${VINE_LEGACY_CPU:-0}" != "1" ]; then
    echo "[post-deploy] VINE_LEGACY_CPU is not 1; keeping the conda MMseqs2 build."
    echo "[post-deploy] Set VINE_LEGACY_CPU=1 before the run on a CPU without AVX2."
    exit 0
fi

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "[post-deploy] ERROR: CONDA_PREFIX is unset; cannot locate the env." >&2
    exit 1
fi

target="${CONDA_PREFIX}/bin/mmseqs"

if grep -qm1 avx2 /proc/cpuinfo 2>/dev/null; then
    echo "[post-deploy] NOTE: this CPU reports avx2, so the swap is not needed"
    echo "[post-deploy] here. Proceeding because VINE_LEGACY_CPU=1 was set."
fi

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

echo "[post-deploy] fetching MMseqs2 ${MMSEQS_VERSION} (static SSE4.1)"
curl -fsSL "$MMSEQS_URL" -o "$workdir/mmseqs.tar.gz"
tar -xzf "$workdir/mmseqs.tar.gz" -C "$workdir"

binary=$(find "$workdir" -type f -name mmseqs -perm -u+x -print -quit)
if [ -z "$binary" ]; then
    echo "[post-deploy] ERROR: no mmseqs binary in the archive from $MMSEQS_URL" >&2
    exit 1
fi

# Verify it actually runs on THIS cpu before replacing a working binary.
if ! "$binary" version > "$workdir/version.txt" 2>&1; then
    echo "[post-deploy] ERROR: the SSE4.1 build does not run here:" >&2
    cat "$workdir/version.txt" >&2
    exit 1
fi

if [ -e "$target" ] && [ ! -e "${target}.conda" ]; then
    mv "$target" "${target}.conda"
fi
install -m 0755 "$binary" "$target"

echo "[post-deploy] mmseqs replaced with $(cat "$workdir/version.txt" | head -1)"
echo "[post-deploy] the conda build is kept at ${target}.conda"
