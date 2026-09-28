#!/usr/bin/env bash
#SBATCH --job-name=protist-allocation
#SBATCH --account=grp-org-sc-mgs
#SBATCH --partition=dori
#SBATCH --qos=jgi_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=384G
#SBATCH --time=72:00:00
#SBATCH --output=tasks/logs/allocation-%j.out
#SBATCH --error=tasks/logs/allocation-%j.err
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
skip_gene_calling=false
skip_annotation=false
prepare_args=()
positional_args=()
for argument in "$@"; do
    case "$argument" in
        --skip-gene-calling)
            if [[ $skip_gene_calling == true ]]; then
                echo 'Duplicate option: --skip-gene-calling' >&2
                exit 2
            fi
            skip_gene_calling=true
            prepare_args+=("$argument")
            ;;
        --skip-annotation)
            if [[ $skip_annotation == true ]]; then
                echo 'Duplicate option: --skip-annotation' >&2
                exit 2
            fi
            skip_annotation=true
            prepare_args+=("$argument")
            ;;
        --*)
            printf 'Unknown option: %s\n' "$argument" >&2
            exit 2
            ;;
        *) positional_args+=("$argument") ;;
    esac
done
set -- "${positional_args[@]}"
if [[ $# != 3 && $# != 4 && $# != 6 ]]; then
    echo 'Usage: run-allocation.slurm SAMPLE_TSV RUN_DIRECTORY core|full [resume | resume-reviewed OLD_DIGEST NEW_DIGEST] [--skip-gene-calling] [--skip-annotation]' >&2
    exit 2
fi
if [[ ! -f "$project_root/pixi.toml" ]]; then
    echo 'Run this script with the project root as the Slurm working directory.' >&2
    exit 2
fi
mode=$3
if [[ $mode != core && $mode != full ]]; then
    echo 'The third argument must be core or full.' >&2
    exit 2
fi
if [[ ($# == 4 && $4 != resume) || ($# == 6 && $4 != resume-reviewed) ]]; then
    echo 'Use resume alone, or resume-reviewed with OLD_DIGEST and NEW_DIGEST.' >&2
    exit 2
fi
if [[ ${4:-} == resume-reviewed ]]; then
    if [[ ! $5 =~ ^sha256:[0-9a-f]{64}$ || ! $6 =~ ^sha256:[0-9a-f]{64}$ ]]; then
        echo 'Reviewed source digests must be sha256: followed by 64 lowercase hex digits.' >&2
        exit 2
    fi
fi
samples=$(realpath "$1")
run_dir=$(realpath -m "$2")
registry="$project_root/conf/databases.yaml"
identity="$run_dir/prepared/resume-identity.tsv"
sample_count=$(awk 'NR > 1 && NF {count++} END {print count + 0}' "$samples")
if [[ $sample_count -ne 1 ]]; then
    echo 'The allocation launcher currently requires exactly one sample.' >&2
    exit 2
fi

: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is required}"
: "${SLURM_MEM_PER_NODE:?SLURM_MEM_PER_NODE is required}"
if [[ ! $SLURM_CPUS_PER_TASK =~ ^[0-9]+$ || ! $SLURM_MEM_PER_NODE =~ ^[0-9]+$ ]]; then
    echo 'Slurm CPU and memory values must be integer counts.' >&2
    exit 2
fi
if ((SLURM_CPUS_PER_TASK < 32)); then
    echo 'The allocation needs at least 32 CPUs.' >&2
    exit 2
fi
allocation_memory_gb=$((SLURM_MEM_PER_NODE / 1024))
executor_memory_gb=$((allocation_memory_gb - 24))
if [[ $mode == core ]]; then
    minimum_allocation_gb=216
    largest_process_memory_gb=192
else
    minimum_allocation_gb=384
    largest_process_memory_gb=256
fi
if ((allocation_memory_gb < minimum_allocation_gb)); then
    echo "$mode mode needs at least ${minimum_allocation_gb} GB of allocation memory." >&2
    exit 2
fi
if ((executor_memory_gb < largest_process_memory_gb)); then
    echo "$mode mode needs at least ${largest_process_memory_gb} GB after the 24 GB executor reserve." >&2
    exit 2
fi

resume_args=()
if [[ ${4:-} == resume || ${4:-} == resume-reviewed ]]; then
    [[ -d "$run_dir/launch/.nextflow" ]]
    [[ -s "$run_dir/prepared/samples.tsv" ]]
    [[ -s "$run_dir/prepared/dependencies.config" ]]
    [[ -s "$run_dir/prepared/inputs.json" ]]
    [[ -s "$identity" ]]
    resume_args=(-resume)
elif [[ -e "$run_dir" ]]; then
    echo 'Run directory already exists; use resume for its saved Nextflow session.' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=2
export NXF_OPTS='-XX:ActiveProcessorCount=2 -Xmx2g'
export SLURM_CLUSTERS=perceus-00
export MPLBACKEND=Agg
export PIXI_NO_INSTALL=true PIXI_FROZEN=true

for executable in python nextflow java; do
    [[ -x "$project_root/.pixi/envs/default/bin/$executable" ]] || {
        echo "Missing installed default Pixi executable: $executable" >&2
        exit 2
    }
done
for executable in chopper seqkit pigz myloasm minimap2 samtools; do
    [[ -x "$project_root/.pixi/envs/core/bin/$executable" ]] || {
        echo "Missing installed core Pixi executable: $executable" >&2
        exit 2
    }
done
for executable in quickbin.sh quickclade.sh bbversion.sh java; do
    [[ -x "$project_root/envs/bbtools/.pixi/envs/bbtools/bin/$executable" ]] || {
        echo "Missing installed BBTools Pixi executable: $executable" >&2
        exit 2
    }
done

pixi run --as-is --manifest-path "$project_root/pixi.toml" \
    bash "$project_root/scripts/require-pixi-tools.sh" \
    python nextflow java awk bash basename cat cmp comm cp date diff dirname env \
    find grep head ln ls mkdir mktemp mv ps realpath rm sed sha256sum sleep sort \
    tail tee touch tr uname uniq wc
pixi run --as-is --manifest-path "$project_root/pixi.toml" --environment core \
    bash "$project_root/scripts/require-pixi-tools.sh" \
    chopper seqkit pigz myloasm minimap2 samtools
pixi run --as-is --manifest-path "$project_root/envs/bbtools/pixi.toml" \
    --environment bbtools bash "$project_root/scripts/require-pixi-tools.sh" \
    quickbin.sh quickclade.sh bbversion.sh java
pixi run --as-is --manifest-path "$project_root/pixi.toml" \
    python - "$registry" "$mode" "$skip_gene_calling" "$skip_annotation" \
    "$run_dir/prepared/inputs.json" "${4:-}" <<'PYTHON'
import json
import sys
from pathlib import Path

from protist_meta.inputs import read_registry

skip_gene_calling = sys.argv[3] == "true"
skip_annotation = sys.argv[4] == "true"
if sys.argv[6] in ("resume", "resume-reviewed"):
    saved = json.loads(Path(sys.argv[5]).read_text(encoding="utf-8"))
    if (
        saved.get("skip_gene_calling", False) != skip_gene_calling
        or saved.get("skip_annotation", False) != skip_annotation
    ):
        raise SystemExit(
            "Resume stage choices differ from saved preparation; use a new run directory."
        )
read_registry(
    Path(sys.argv[1]),
    mode=sys.argv[2],
    skip_gene_calling=skip_gene_calling,
    skip_annotation=skip_annotation,
)
PYTHON

if [[ ${4:-} == resume || ${4:-} == resume-reviewed ]]; then
    current_source_revision=$(pixi run --as-is --manifest-path "$project_root/pixi.toml" \
        python -c 'import sys; from pathlib import Path; from protist_meta.inputs import _source_digest; print(_source_digest(Path(sys.argv[1])))' \
        "$project_root")
    expected_source_revision=$current_source_revision
    if [[ ${4:-} == resume-reviewed ]]; then
        if [[ $6 != "$current_source_revision" ]]; then
            echo 'Reviewed execution source digest does not match the current source.' >&2
            exit 2
        fi
        expected_source_revision=$5
    fi
    while IFS=$'\t' read -r key current; do
        saved=$(awk -F '\t' -v key="$key" '$1 == key {print $2}' "$identity")
        if [[ -z "$saved" || "$saved" != "$current" ]]; then
            echo "Resume identity mismatch: $key" >&2
            exit 2
        fi
    done <<IDENTITY
run_mode	$mode
input_samples_path	$samples
input_samples_sha256	$(sha256sum "$samples" | awk '{print $1}')
registry_sha256	$(sha256sum "$registry" | awk '{print $1}')
prepared_samples_sha256	$(sha256sum "$run_dir/prepared/samples.tsv" | awk '{print $1}')
dependencies_config_sha256	$(sha256sum "$run_dir/prepared/dependencies.config" | awk '{print $1}')
inputs_json_sha256	$(sha256sum "$run_dir/prepared/inputs.json" | awk '{print $1}')
source_revision	$expected_source_revision
IDENTITY
else
    mkdir -p "$run_dir/launch"
    pixi run --as-is --manifest-path "$project_root/pixi.toml" \
        python -m protist_meta.cli prepare \
        --input "$samples" --registry "$registry" \
        --output-dir "$run_dir/prepared" --mode "$mode" "${prepare_args[@]}"
    source_revision=$(awk -F "'" '/^params.source_revision = / {print $2}' \
        "$run_dir/prepared/dependencies.config")
    [[ -n "$source_revision" ]]
    {
        printf 'run_mode\t%s\n' "$mode"
        printf 'input_samples_path\t%s\n' "$samples"
        printf 'input_samples_sha256\t%s\n' "$(sha256sum "$samples" | awk '{print $1}')"
        printf 'registry_sha256\t%s\n' "$(sha256sum "$registry" | awk '{print $1}')"
        printf 'prepared_samples_sha256\t%s\n' \
            "$(sha256sum "$run_dir/prepared/samples.tsv" | awk '{print $1}')"
        printf 'dependencies_config_sha256\t%s\n' \
            "$(sha256sum "$run_dir/prepared/dependencies.config" | awk '{print $1}')"
        printf 'inputs_json_sha256\t%s\n' \
            "$(sha256sum "$run_dir/prepared/inputs.json" | awk '{print $1}')"
        printf 'source_revision\t%s\n' "$source_revision"
    } > "$identity"
fi
mkdir -p "$run_dir/executions"
execution_dir=$(mktemp -d "$run_dir/executions/${SLURM_JOB_ID}.XXXXXX")
cat > "$execution_dir/allocation.config" <<CONFIG
executor.cpus = ${SLURM_CPUS_PER_TASK}
executor.memory = '${executor_memory_gb} GB'
trace.file = '${execution_dir}/trace.tsv'
trace.fields = 'task_id,hash,native_id,name,status,exit,cpus,memory,realtime,peak_rss,workdir'
report.file = '${execution_dir}/nextflow-report.html'
timeline.file = '${execution_dir}/timeline.html'
dag.file = '${execution_dir}/dag.html'
CONFIG
if [[ ${4:-} == resume-reviewed ]]; then
    pixi run --as-is --manifest-path "$project_root/pixi.toml" \
        python - "$identity" "$execution_dir" "$5" "$6" <<'PYTHON'
import csv
import hashlib
import json
import re
import sys
from datetime import UTC, datetime
from pathlib import Path

identity = Path(sys.argv[1])
execution = Path(sys.argv[2])
with identity.open(encoding="utf-8", newline="") as handle:
    preparation_identity = dict(csv.reader(handle, delimiter="\t"))
inputs = json.loads((identity.parent / "inputs.json").read_text(encoding="utf-8"))
dependencies = (identity.parent / "dependencies.config").read_text(encoding="utf-8")
run_id = re.search(r"^params\.run_id = '([^']+)'$", dependencies, re.MULTILINE)
if run_id is None:
    raise ValueError("saved preparation has no run_id")
forced_processes = ["ROUTE_BINS", "COLLECT_RECORDS", "BUILD_CATALOG", "BUILD_REPORT"]
receipt_path = execution / "recovery.json"
config_path = execution / "allocation.config"
quoted_receipt = str(receipt_path).replace("\\", "\\\\").replace("'", "\\'")
with config_path.open("a", encoding="utf-8") as handle:
    for process in forced_processes:
        handle.write(f"process {{ withName: '{process}' {{ cache = false }} }}\n")
    handle.write(f"params.execution_provenance = '{quoted_receipt}'\n")
receipt = {
    "schema_version": 1,
    "recovery_mode": "reviewed-source",
    "sample_id": inputs["samples"][0]["sample_id"],
    "run_id": run_id.group(1),
    "original_source_revision": sys.argv[3],
    "execution_source_revision": sys.argv[4],
    "execution_directory": str(execution),
    "started_at": datetime.now(UTC).isoformat(),
    "forced_processes": forced_processes,
    "preparation_identity": preparation_identity,
    "preparation_identity_sha256": hashlib.sha256(identity.read_bytes()).hexdigest(),
    "allocation_config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
}
receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
PYTHON
fi
cd "$run_dir/launch"
pixi run --as-is --manifest-path "$project_root/pixi.toml" \
    nextflow -log "$execution_dir/nextflow.log" run "$project_root/main.nf" \
    -offline -profile integration \
    -c "$run_dir/prepared/dependencies.config" -c "$execution_dir/allocation.config" \
    -work-dir "$run_dir/work" \
    --input "$run_dir/prepared/samples.tsv" --outdir "$run_dir/results" \
    --run_mode "$mode" "${resume_args[@]}"
test -s "$run_dir/results/catalog/protist-meta.sqlite"
test -s "$run_dir/results/report/report.executed.ipynb"
test -s "$run_dir/results/report/index.html"
