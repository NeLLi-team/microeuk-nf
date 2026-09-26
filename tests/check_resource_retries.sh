#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo 'Usage: bash tests/check_resource_retries.sh PROOF_DIR' >&2
    exit 2
fi
if [[ -z ${SLURM_JOB_ID:-} ]]; then
    echo 'Resource retry checks must run inside Slurm.' >&2
    exit 2
fi

root=$PWD
proof=$(realpath -m "$1")
ssuextract=/clusterfs/jgi/scratch/science/mgs/nelli/frederik/projects/apps/ssuextract
root_bash="$root/.pixi/envs/default/bin/bash"
module="$root/modules/characterization.nf"
nested_config="$proof/generated-nested.config"

[[ -f "$root/nextflow.config" && -f "$root/conf/resources.config" ]]
[[ -f "$module" && -f "$ssuextract/nextflow.config" ]]
[[ -x "$root_bash" ]]
mkdir -p "$proof"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=1
export NXF_ANSI_LOG=false NXF_OFFLINE=true
export NXF_OPTS='-XX:ActiveProcessorCount=1 -Xms128m -Xmx512m'
export PIXI_NO_INSTALL=true PIXI_FROZEN=true

run_root() {
    pixi run --as-is --quiet --manifest-path "$root/pixi.toml" "$@"
}

run_ssu() {
    pixi run --as-is --quiet --manifest-path "$ssuextract/pixi.toml" "$@"
}

run_root nextflow -version > "$proof/root-nextflow-version.txt"
run_ssu nextflow -version > "$proof/ssu-nextflow-version.txt"
printf 'Nextflow versions recorded in %s\n' "$proof"

# Extract the config that SSU_EXTRACT writes, replacing only its outer-task resources.
awk '
    /LAUNCH_ROOT.*nested\.config.*CONFIG/ {capture = 1; next}
    capture && /^CONFIG$/ {exit}
    capture {print}
' "$module" | sed \
        -e "s#\${projectDir}#$root#g" \
        -e 's/${nestedCpus}/1/g' \
        -e 's/${nestedMemoryGb}/4/g' \
        > "$nested_config"

nested_lines=$(wc -l < "$nested_config")
if ((nested_lines == 0)); then
    echo 'Failed to extract generated nested.config from SSU_EXTRACT.' >&2
    exit 1
fi
grep -F "withName: BLAST_ANNOTATE" "$nested_config" >/dev/null
grep -F 'time = { [8.h * task.attempt, params.max_time as nextflow.util.Duration].min() }' \
    "$nested_config" >/dev/null
grep -F 'maxRetries = 1' "$nested_config" >/dev/null
grep -F 'def error = task.previousException' "$nested_config" >/dev/null
grep -F 'def cause = error?.cause ?: error' "$nested_config" >/dev/null
grep -F 'cause instanceof nextflow.exception.ProcessException' "$nested_config" >/dev/null
grep -F "cause.message?.startsWith('Process exceeded running time limit')" \
    "$nested_config" >/dev/null
printf 'generated nested.config: %s lines with BLAST_ANNOTATE policy\n' "$nested_lines"

mkdir -p "$proof/config-root" "$proof/config-ssu"
cp "$nested_config" "$proof/config-ssu/nested.config"
(
    cd "$proof/config-root"
    run_root nextflow -log "$proof/root-config.log" \
        -C "$root/nextflow.config" config "$root" -profile integration -flat \
        > "$proof/root-resolved.config"
)
(
    cd "$proof/config-ssu"
    run_ssu nextflow -log "$proof/ssu-config.log" \
        -c "$proof/config-ssu/nested.config" \
        config "$ssuextract" -profile local -flat \
        > "$proof/ssu-resolved.config"
)

grep -F "process.'withName:READ_QC'.time = { 8.h * task.attempt }" \
    "$proof/root-resolved.config"
grep -F "process.'withName:READ_QC'.maxRetries = 1" \
    "$proof/root-resolved.config"
grep -F "process.'withName:BLAST_ANNOTATE'.time = { [8.h * task.attempt, params.max_time as nextflow.util.Duration].min() }" \
    "$proof/ssu-resolved.config"
grep -F "process.'withName:BLAST_ANNOTATE'.maxRetries = 1" \
    "$proof/ssu-resolved.config"
printf 'default budgets: READ_QC=8h/16h BLAST_ANNOTATE=8h/16h\n'

write_flow() {
    local process_name=$1
    local path=$2
    sed "s/PROCESS_NAME/$process_name/g" > "$path" <<'NEXTFLOW'
nextflow.enable.dsl = 2
params.mode = null

process PROCESS_NAME {
    output:
    path 'result.txt'

    script:
    """
    if [[ '${params.mode}' == 'exit2' ]]; then
        exit 2
    fi
    sleep 2
    printf 'attempt=%s\ntime_hours=%s\n' \
        '${task.attempt}' '${task.time.toHours()}' > result.txt
    """
}

workflow {
    PROCESS_NAME()
}
NEXTFLOW
}

write_case_config() {
    local process_name=$1
    local trace=$2
    local path=$3
    local time_override=$4
    cat > "$path" <<CONFIG
process.shell = ['$root_bash', '-euo', 'pipefail']
process.executor = 'local'
process.maxForks = 1
executor.cpus = 1
executor.memory = '4 GB'
executor.pollInterval = '100ms'
process {
    withName: $process_name {
        cpus = 1
        memory = '1 GB'
        $time_override
    }
}
trace.enabled = true
trace.file = '$trace'
trace.fields = 'name,status,exit,attempt,realtime'
timeline.enabled = false
report.enabled = false
dag.enabled = false
CONFIG
}

check_retry_trace() {
    awk -F '\t' '
        NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
        {n++; status[$c["attempt"]] = $c["status"]; code[$c["attempt"]] = $c["exit"]}
        END {
            ok = n == 2 && status[1] == "FAILED" && code[1] == "-" &&
                status[2] == "COMPLETED" && code[2] == "0"
            exit !ok
        }
    ' "$1"
}

check_exit_trace() {
    awk -F '\t' '
        NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
        {n++; attempt = $c["attempt"]; status = $c["status"]; code = $c["exit"]}
        END {exit !(n == 1 && attempt == 1 && status == "FAILED" && code == "2")}
    ' "$1"
}

check_default_result() {
    local case_dir=$1
    local trace=$2
    local expected_hours=$3
    local -a results
    awk -F '\t' '
        NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
        {n++; attempt = $c["attempt"]; status = $c["status"]; code = $c["exit"]}
        END {exit !(n == 1 && attempt == 1 && status == "COMPLETED" && code == "0")}
    ' "$trace"
    mapfile -t results < <(find "$case_dir/work" -type f -name result.txt)
    if [[ ${#results[@]} -ne 1 ]]; then
        echo "Expected one result.txt for default-budget case, found ${#results[@]}." >&2
        exit 1
    fi
    cp "${results[0]}" "$case_dir/result.txt"
    grep -Fx 'attempt=1' "$case_dir/result.txt" >/dev/null
    grep -Fx "time_hours=$expected_hours" "$case_dir/result.txt" >/dev/null
}

run_case() {
    local engine=$1
    local process_name=$2
    local mode=$3
    local expected_hours=${4:-8}
    local max_time=${5:-24.h}
    local case_dir="$proof/$engine-$mode"
    local flow="$case_dir/main.nf"
    local config="$case_dir/harness.config"
    local trace="$case_dir/trace.tsv"
    local time_override='time = { task.attempt == 1 ? 1.s : 4.s }'
    if [[ $mode == default || $mode == capped ]]; then
        time_override=
    fi
    mkdir -p "$case_dir/launch"
    write_flow "$process_name" "$flow"
    write_case_config "$process_name" "$trace" "$config" "$time_override"

    local -a command
    if [[ $engine == root ]]; then
        command=(run_root nextflow -log "$case_dir/nextflow.log"
            -C "$root/nextflow.config,$config" run "$flow" -offline
            -work-dir "$case_dir/work" --mode "$mode")
    else
        command=(run_ssu nextflow -log "$case_dir/nextflow.log"
            -C "$ssuextract/nextflow.config,$nested_config,$config"
            run "$flow" -offline -profile local -work-dir "$case_dir/work"
            --mode "$mode" --threads_per_job 1 --max_cpus 1
            --max_memory 4.GB --max_time "$max_time")
    fi

    if [[ $mode == default || $mode == capped ]]; then
        (cd "$case_dir/launch"; "${command[@]}") \
            > "$case_dir/stdout.log" 2> "$case_dir/stderr.log"
        check_default_result "$case_dir" "$trace" "$expected_hours"
        printf '%s %s: attempt1 COMPLETED with task.time=%sh\n' \
            "$process_name" "$mode" "$expected_hours"
    elif [[ $mode == timeout ]]; then
        (cd "$case_dir/launch"; "${command[@]}") \
            > "$case_dir/stdout.log" 2> "$case_dir/stderr.log"
        check_retry_trace "$trace"
        printf '%s timeout: attempt1 FAILED exit=-, attempt2 COMPLETED exit=0\n' "$process_name"
    else
        if (cd "$case_dir/launch"; "${command[@]}") \
            > "$case_dir/stdout.log" 2> "$case_dir/stderr.log"; then
            echo "$process_name exit2 case unexpectedly succeeded" >&2
            exit 1
        fi
        check_exit_trace "$trace"
        printf '%s exit2: one FAILED attempt with exit=2\n' "$process_name"
    fi
}

run_case root READ_QC default
run_case ssu BLAST_ANNOTATE default
run_case ssu BLAST_ANNOTATE capped 4 4.h
run_case root READ_QC timeout
run_case root READ_QC exit2
run_case ssu BLAST_ANNOTATE timeout
run_case ssu BLAST_ANNOTATE exit2

printf 'PASS\n' > "$proof/status.txt"
printf 'PASS resource retry policies\n'
