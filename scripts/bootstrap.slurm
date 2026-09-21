#!/usr/bin/env bash
#SBATCH --job-name=protist-bootstrap
#SBATCH --account=grp-org-sc-mgs
#SBATCH --partition=dori
#SBATCH --qos=jgi_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tasks/logs/bootstrap-%j.out
#SBATCH --error=tasks/logs/bootstrap-%j.err
set -euo pipefail

project_root=$PWD
root_prefix="$project_root/.pixi/envs/default"
root_bash="$root_prefix/bin/bash"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=2
export PIXI_CONCURRENT_DOWNLOADS=2 PIXI_CONCURRENT_SOLVES=1
export PIXI_CACHE_DIR="${SLURM_TMPDIR:-/tmp}/protist-meta-nf-pixi-${USER}"
if [[ ${PROTIST_BOOTSTRAP_REEXEC:-0} != 1 ]]; then
    [[ -f "$project_root/pixi.toml" ]]
    pixi install --manifest-path "$project_root/pixi.toml"
    if [[ ! -x "$root_bash" ]]; then
        printf 'Missing installed root Pixi Bash: %s\n' "$root_bash" >&2
        exit 2
    fi
    export PROTIST_BOOTSTRAP_REEXEC=1 PIXI_NO_INSTALL=true PIXI_FROZEN=true
    pixi run --as-is --manifest-path "$project_root/pixi.toml" \
        "$root_bash" "$0" "$@"
    exit "$?"
fi
if [[ $BASH != "$root_bash" || ${CONDA_PREFIX:-} != "$root_prefix" \
    || ${PIXI_PROJECT_ROOT:-} != "$project_root" ]]; then
    printf 'Bootstrap did not enter the root Pixi environment: %s\n' "$root_prefix" >&2
    exit 2
fi
unset PROTIST_BOOTSTRAP_REEXEC PIXI_NO_INSTALL PIXI_FROZEN

pixi install --manifest-path "$project_root/pixi.toml" --environment core --locked
export PIXI_NO_INSTALL=true PIXI_FROZEN=true
pixi run --as-is --manifest-path "$project_root/pixi.toml" \
    bash "$project_root/scripts/require-pixi-tools.sh" \
    python nextflow java awk bash basename cat cmp comm cp date diff dirname env \
    find grep head ln ls mkdir mktemp mv ps realpath rm sed sha256sum sleep sort \
    tail tee touch tr uname uniq wc
pixi run --as-is --environment core bash scripts/require-pixi-tools.sh \
    chopper seqkit pigz myloasm minimap2 samtools
pixi run --as-is python --version
pixi run --as-is nextflow -version
pixi run --as-is --environment core myloasm --help
