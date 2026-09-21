#!/usr/bin/env bash
#SBATCH --job-name=protist-tools
#SBATCH --account=grp-org-sc-mgs
#SBATCH --partition=dori
#SBATCH --qos=jgi_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tasks/logs/tools-%j.out
#SBATCH --error=tasks/logs/tools-%j.err
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
export RAYON_NUM_THREADS=2 PIXI_CONCURRENT_DOWNLOADS=2 PIXI_CONCURRENT_SOLVES=1
export PIXI_CACHE_DIR="${SLURM_TMPDIR:-/tmp}/protist-meta-nf-pixi-${USER}"
for tool_environment in bbtools:bbtools gene:repeatmasker gene:prodigal-gv annotation:eggnog; do
    manifest="envs/${tool_environment%%:*}/pixi.toml"
    environment="${tool_environment##*:}"
    pixi install --manifest-path "$manifest" --environment "$environment" --locked
done
export PIXI_NO_INSTALL=true PIXI_FROZEN=true
dfam_db=$(
    pixi run --as-is --quiet python - conf/databases.yaml <<'PYTHON'
import sys
from pathlib import Path

import yaml

registry_path = Path(sys.argv[1]).resolve()
registry = yaml.safe_load(registry_path.read_text(encoding="utf-8"))
value = registry["entries"]["dfam_db"]["path"]
path = Path(value).expanduser()
print((path if path.is_absolute() else registry_path.parent / path).resolve())
PYTHON
)
bash scripts/configure-famdb.sh write envs/gene/pixi.toml "$dfam_db"
bash scripts/configure-famdb.sh check envs/gene/pixi.toml "$dfam_db"
pixi run --as-is --quiet --manifest-path envs/bbtools/pixi.toml \
    --environment bbtools bbversion.sh
pixi run --as-is --quiet --manifest-path envs/gene/pixi.toml \
    --environment repeatmasker RepeatMasker -v
pixi run --as-is --quiet --manifest-path envs/gene/pixi.toml \
    --environment repeatmasker famdb.py info
pixi run --as-is --quiet --manifest-path envs/gene/pixi.toml \
    --environment prodigal-gv prodigal-gv -v
pixi run --as-is --quiet --manifest-path envs/annotation/pixi.toml --environment eggnog \
    emapper.py --version --data_dir /clusterfs/jgi/scratch/science/mgs/nelli/databases/eggnogmapper
