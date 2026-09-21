#!/usr/bin/env bash
#SBATCH --job-name=protist-dfam
#SBATCH --account=grp-org-sc-mgs
#SBATCH --partition=dori
#SBATCH --qos=jgi_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=tasks/logs/dfam-%j.out
#SBATCH --error=tasks/logs/dfam-%j.err
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

[[ -s "$project_root/conf/databases.yaml" ]]
[[ -s "$project_root/envs/gene/pixi.toml" ]]
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=1
export PIXI_CONCURRENT_DOWNLOADS=2 PIXI_CONCURRENT_SOLVES=1
export PIXI_CACHE_DIR="${SLURM_TMPDIR:-/tmp}/protist-meta-dfam-${USER}"

pixi install --manifest-path "$project_root/envs/gene/pixi.toml" \
    --environment repeatmasker --locked
export PIXI_NO_INSTALL=true PIXI_FROZEN=true

dfam_db=$(
    pixi run --as-is --quiet --manifest-path "$project_root/pixi.toml" \
        python - "$project_root/conf/databases.yaml" <<'PYTHON'
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

pixi run --as-is --quiet --manifest-path "$project_root/pixi.toml" \
    python "$project_root/scripts/install-dfam.py" "$dfam_db"
bash "$project_root/scripts/configure-famdb.sh" write \
    "$project_root/envs/gene/pixi.toml" "$dfam_db"
bash "$project_root/scripts/configure-famdb.sh" check \
    "$project_root/envs/gene/pixi.toml" "$dfam_db"
pixi run --as-is --quiet --manifest-path "$project_root/envs/gene/pixi.toml" \
    --environment repeatmasker famdb.py info > "$dfam_db/famdb-info.txt"
grep -Eq '^Database[[:space:]]*:[[:space:]]*Dfam$' "$dfam_db/famdb-info.txt"
grep -Eq '^Version[[:space:]]*:[[:space:]]*4\.0$' "$dfam_db/famdb-info.txt"
