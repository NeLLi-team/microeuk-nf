#!/usr/bin/env bash
#SBATCH --job-name=protist-custom
#SBATCH --account=grp-org-sc-mgs
#SBATCH --partition=dori
#SBATCH --qos=jgi_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=tasks/logs/custom-%j.out
#SBATCH --error=tasks/logs/custom-%j.err
set -euo pipefail

project_root=$PWD
root_prefix="$project_root/.pixi/envs/default"
root_bash="$root_prefix/bin/bash"
if [[ $BASH != "$root_bash" || ${CONDA_PREFIX:-} != "$root_prefix" \
    || ${PIXI_PROJECT_ROOT:-} != "$project_root" ]]; then
    if [[ ! -x "$root_bash" ]]; then
        printf 'Missing installed root Pixi Bash: %s\n' "$root_bash" >&2
        exit 2
    fi
    export PIXI_NO_INSTALL=true PIXI_FROZEN=true
    pixi run --as-is --manifest-path "$project_root/pixi.toml" \
        "$root_bash" "$0" "$@"
    exit "$?"
fi
unset PIXI_NO_INSTALL PIXI_FROZEN

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=2
export PIXI_CONCURRENT_DOWNLOADS=2 PIXI_CONCURRENT_SOLVES=1
export PIXI_CACHE_DIR="${SLURM_TMPDIR:-/tmp}/protist-meta-nf-custom-${USER}"
proof="$PWD/tasks/custom-${SLURM_JOB_ID}"
mkdir -p "$proof/interpro/output" "$proof/interpro/tmp"

pixi install --manifest-path envs/gene/longrna/pixi.toml
pixi run --as-is --manifest-path envs/gene/longrna/pixi.toml \
    check-cli > "$proof/braker-cli.txt" 2>&1
pixi install --manifest-path envs/annotation/current/pixi.toml
export INTERPROSCAN_DATA_DIR=/clusterfs/jgi/scratch/science/mgs/nelli/databases/interproscan/data
pixi run --as-is --manifest-path envs/annotation/current/pixi.toml \
    check-cli > "$proof/interpro-identity.txt" 2>&1
pixi run --as-is --manifest-path envs/annotation/current/pixi.toml \
    interproscan.sh --version > "$proof/interpro-version.txt" 2>&1
cp /clusterfs/jgi/scratch/science/mgs/nelli/databases/interproscan/test_single_protein.fasta \
    "$proof/interpro/input.fasta"
export INTERPROSCAN_TEMP_DIR="$proof/interpro/tmp"
pixi run --as-is --manifest-path envs/annotation/current/pixi.toml \
    interproscan.sh -i "$proof/interpro/input.fasta" \
    -d "$proof/interpro/output" -f TSV,GFF3 -dp -cpu "$SLURM_CPUS_PER_TASK" \
    > "$proof/interpro-smoke.txt" 2>&1
test "$(find "$proof/interpro/output" -type f -name '*.tsv' -size +0c | wc -l)" -gt 0
mapfile -t positive_gffs < <(find "$proof/interpro/output" -type f -name '*.gff3')
[[ ${#positive_gffs[@]} -eq 1 ]]
grep -q '^##gff-version 3' "${positive_gffs[0]}"
awk -F '\t' '!/^#/ && NF==9 {features++} END {exit features==0}' "${positive_gffs[0]}"
mkdir -p "$proof/interpro/no-signature"
printf '>no_signature_probe\nMPEPTIDE\n' > "$proof/interpro/no-signature.fasta"
pixi run --as-is --manifest-path envs/annotation/current/pixi.toml \
    interproscan.sh -i "$proof/interpro/no-signature.fasta" \
    -d "$proof/interpro/no-signature" -f TSV,GFF3 -dp -cpu "$SLURM_CPUS_PER_TASK" \
    > "$proof/interpro-no-signature.txt" 2>&1
mapfile -t empty_reports < <(find "$proof/interpro/no-signature" -type f -name '*.tsv')
[[ ${#empty_reports[@]} -eq 1 ]]
[[ ! -s "${empty_reports[0]}" ]]
mapfile -t negative_gffs < <(find "$proof/interpro/no-signature" -type f -name '*.gff3')
[[ ${#negative_gffs[@]} -eq 1 ]]
grep -q '^##gff-version 3' "${negative_gffs[0]}"
pixi run --as-is --manifest-path envs/annotation/current/pixi.toml \
    check-cli > "$proof/interpro-identity-after.txt" 2>&1
cmp "$proof/interpro-identity.txt" "$proof/interpro-identity-after.txt"
bash scripts/check-genemark.sh "$proof/genemark" > "$proof/genemark-smoke.txt" 2>&1
printf 'passed\n' > "$proof/status.txt"
