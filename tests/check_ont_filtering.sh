#!/usr/bin/env bash
# Run inside a full-workflow allocation with an installed MICRO sample sheet.
set -euo pipefail
project_root=$PWD
samples=$(realpath "${1:?Usage: check_ont_filtering.sh MICRO_SAMPLE_TSV PROOF_DIRECTORY}")
proof=$(realpath -m "${2:?Supply a new proof directory}")
mkdir "$proof"
export PIXI_NO_INSTALL=true PIXI_FROZEN=true

awk '
    function read_record(name, length_bp, quality, sequence, qualities, i) {
        for (i = 0; i < length_bp; i++) {
            sequence = sequence "A"
            qualities = qualities quality
        }
        printf "@%s\n%s\n+\n%s\n", name, sequence, qualities
    }
    BEGIN {
        read_record("short_q40", 999, "I")
        read_record("boundary_q40", 1000, "I")
        read_record("long_q9", 2000, "*")
        read_record("long_q11", 2000, ",")
        read_record("long_q40", 2000, "I")
    }
' > "$proof/boundaries.fastq"
pixi run --as-is --manifest-path "$project_root/pixi.toml" --environment core \
    chopper -q 10 -l 1000 --threads 1 < "$proof/boundaries.fastq" \
    > "$proof/filtered.fastq" 2> "$proof/chopper.stderr"
awk 'NR % 4 == 1 {print substr($0, 2)}' "$proof/filtered.fastq" > "$proof/retained.txt"
printf 'boundary_q40\nlong_q11\nlong_q40\n' > "$proof/expected.txt"
cmp "$proof/expected.txt" "$proof/retained.txt"

run_dir="$proof/workflow"
bash scripts/run-allocation.slurm "$samples" "$run_dir" full --ont-min-length 1000
sample=$(awk -F '\t' 'NR == 2 {print $1}' "$samples")
qc="$run_dir/results/$sample/01_read_qc"
grep -F '"command": "chopper -q 10 -l 1000 --threads 1"' "$qc/stage.json"
awk -F '\t' 'NR == 2 {ok=($4 > 0 && $6 >= 1000)} END {exit !ok}' "$qc/filtered_stats.tsv"
grep -Fx $'ont_min_length\t1000' "$run_dir/prepared/resume-identity.tsv"
sha256sum "$run_dir/results/catalog/protist-meta.sqlite" \
    "$run_dir/results/report/report.executed.ipynb" \
    "$run_dir/results/report/index.html" > "$proof/before-resume.sha256"
bash scripts/run-allocation.slurm "$samples" "$run_dir" full resume --ont-min-length 1000
sha256sum --check "$proof/before-resume.sha256" > "$proof/resume-check.txt"

for requested in omitted 500; do
    options=()
    if [[ "$requested" != omitted ]]; then
        options=(--ont-min-length "$requested")
    fi
    if bash scripts/run-allocation.slurm "$samples" "$run_dir" full resume \
        "${options[@]}" > "$proof/rejected-$requested.log" 2>&1; then
        printf 'Unexpectedly resumed with %s override.\n' "$requested" >&2
        exit 1
    fi
    grep -F 'Resume ONT minimum length differs' "$proof/rejected-$requested.log"
done
for invalid in 0 -1 1.5 ''; do
    if bash scripts/run-allocation.slurm "$samples" "$run_dir" full \
        --ont-min-length "$invalid" > "$proof/rejected-argument.log" 2>&1; then
        printf 'Unexpectedly accepted invalid minimum: %s\n' "$invalid" >&2
        exit 1
    fi
    grep -F 'requires a positive integer in bases' "$proof/rejected-argument.log"
done
printf 'PASS\n' > "$proof/status.txt"
