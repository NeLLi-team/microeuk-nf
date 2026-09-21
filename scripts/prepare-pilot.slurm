#!/usr/bin/env bash
#SBATCH --job-name=protist-prepare-pilot
#SBATCH --account=grp-org-sc-mgs
#SBATCH --partition=dori
#SBATCH --qos=jgi_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=tasks/logs/subset-%j.out
#SBATCH --error=tasks/logs/subset-%j.err
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
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=1
export PIXI_NO_INSTALL=true PIXI_FROZEN=true
source_reads="$project_root/../09enrichments-ont/01_shared_preprocessing/results/X0341/reads/subsets/X0341.clean.min1kb.p10.fastq.gz"
output="$project_root/data/pilot"
[[ -s "$source_reads" ]]
[[ ! -e "$output/X0341.p01.fastq.gz" ]]
for executable in seqkit pigz; do
    [[ -x "$project_root/.pixi/envs/core/bin/$executable" ]] || {
        echo "Missing installed core Pixi executable: $executable" >&2
        exit 2
    }
done
pixi run --as-is --manifest-path "$project_root/pixi.toml" --environment core \
    bash "$project_root/scripts/require-pixi-tools.sh" seqkit pigz
mkdir -p "$output"
pixi run --as-is --manifest-path "$project_root/pixi.toml" --environment core \
    seqkit sample -j 1 -s 42 -p 0.1 "$source_reads" \
    | pixi run --as-is --manifest-path "$project_root/pixi.toml" --environment core \
        pigz -p 2 > "$output/X0341.p01.fastq.gz.partial"
pixi run --as-is --manifest-path "$project_root/pixi.toml" --environment core \
    seqkit stats -j 1 -T "$output/X0341.p01.fastq.gz.partial" > "$output/read-stats.tsv"
[[ $(awk 'NR==2 {print $4}' "$output/read-stats.tsv") -gt 0 ]]
mv "$output/X0341.p01.fastq.gz.partial" "$output/X0341.p01.fastq.gz"
printf 'sample_id\tplatform\treads\tbasecaller_model\nX0341_p01\tONT\tX0341.p01.fastq.gz\tdna_r10.4.1_e8.2_400bps_sup@v5.2.0\n' > "$output/samples.tsv"
sha256sum "$output/X0341.p01.fastq.gz" > "$output/sha256.txt"
printf 'source\t%s\nsource_fraction\t0.1\nretained_fraction\t0.1\nseed\t42\njob_id\t%s\n' \
    "$source_reads" "$SLURM_JOB_ID" > "$output/provenance.tsv"
