#!/usr/bin/env bash
#SBATCH --job-name=protist-validate
#SBATCH --account=grp-org-sc-mgs
#SBATCH --partition=dori
#SBATCH --qos=jgi_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=tasks/logs/validate-%j.out
#SBATCH --error=tasks/logs/validate-%j.err
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
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=2
export NXF_OPTS='-XX:ActiveProcessorCount=2 -Xmx2g'
export MPLBACKEND=Agg
export PIXI_NO_INSTALL=true PIXI_FROZEN=true

for executable in python nextflow java ruff mkdocs; do
    [[ -x "$PWD/.pixi/envs/default/bin/$executable" ]] || {
        echo "Missing installed default Pixi executable: $executable" >&2
        exit 2
    }
done

proof="$project_root/tasks/validation-${SLURM_JOB_ID}"
mkdir -p "$proof"
pixi run --as-is bash scripts/require-pixi-tools.sh \
    python nextflow java ruff mkdocs awk bash basename cat chmod cmp comm cp date diff \
    dirname env find grep head ln ls mkdir mktemp mv ps realpath rm sed sha256sum \
    sleep sort tail tee touch tr uname uniq wc xargs \
    > "$proof/pixi-tools.tsv"
status=0
# Check the committed models before regenerating them: after a regeneration the check cannot fail.
if ! pixi run --as-is python scripts/generate_models.py --check > "$proof/schema-check.log" 2>&1; then
    status=1
fi
pixi run --as-is python scripts/generate_models.py > "$proof/generate.log" 2>&1
mapfile -t sources < <(find src scripts tests -type f -name '*.py' | LC_ALL=C sort)
if ! pixi run --as-is python /clusterfs/jgi/groups/science/homes/fschulz/.agents/skills/python-style/scripts/check.py "${sources[@]}" > "$proof/python-style.log" 2>&1; then
    status=1
fi
if ! pixi run --as-is python -m pytest -q > "$proof/pytest.log" 2>&1; then
    status=1
fi
if ! pixi run --as-is bash tests/check_prodigal_parallel.sh "$proof/prodigal-parallel" \
    > "$proof/prodigal-parallel.log" 2>&1; then
    status=1
fi
if ! pixi run --as-is bash tests/check_stage_merge.sh "$proof/stage-merge" \
    > "$proof/stage-merge.log" 2>&1; then
    status=1
fi
if ! pixi run --as-is mkdocs build --strict > "$proof/docs.log" 2>&1; then
    status=1
fi
if ! pixi run --as-is nextflow config -flat > "$proof/nextflow.config" 2> "$proof/nextflow-config.err"; then
    status=1
fi
if ! pixi run --as-is nextflow lint main.nf modules nextflow.config conf > "$proof/nextflow-lint.log" 2>&1; then
    status=1
fi
printf '@preview\nACGT\n+\nIIII\n' > "$proof/preview.fastq"
printf 'sample_id\tplatform\treads\nPREVIEW\tONT\t%s\n' \
    "$proof/preview.fastq" > "$proof/preview-samples.tsv"
# Gene calling is off by default. A boolean override needs a params file: the command line
# passes `--skip_gene_calling false` as a string, which main.nf rejects.
printf '{"skip_gene_calling": false}\n' > "$proof/preview-gene-calling.json"
printf '{"skip_gene_calling": false, "skip_repeat_masking": false}\n' \
    > "$proof/preview-repeat-masking.json"
run_preview() {
    local name=$1 mode=$2
    shift 2
    mkdir -p "$proof/preview-$name"
    if ! (
        cd "$proof/preview-$name"
        pixi run --as-is --manifest-path "$project_root/pixi.toml" \
            nextflow -log "$proof/preview-$name.log" run "$project_root/main.nf" \
            -preview -offline --input "$proof/preview-samples.tsv" --run_mode "$mode" "$@"
    ) > "$proof/preview-$name.stdout" 2>&1; then
        status=1
    fi
}
run_preview core core
run_preview full full
run_preview full-gene-calling full -params-file "$proof/preview-gene-calling.json"
run_preview full-repeat-masking full -params-file "$proof/preview-repeat-masking.json"
sha256sum "${sources[@]}" schema/catalog.yaml pixi.lock > "$proof/source.sha256"
printf '%s\n' "$status" > "$proof/exit-status.txt"
exit "$status"
