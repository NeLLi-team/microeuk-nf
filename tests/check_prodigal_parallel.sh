#!/usr/bin/env bash
# Native integration proof; run from the repository in an allocation of >=2 CPUs.
set -euo pipefail
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_CPUS_PER_TASK:-0} -ge 2 ]] || {
    echo 'Prodigal parallel checks require a Slurm allocation with at least 2 CPUs.' >&2
    exit 2
}
root=$PWD
proof=$(realpath -m "${1:?Usage: check_prodigal_parallel.sh NEW_PROOF_DIR [SERIAL_MODULE]}")
reference=$(realpath "${2:-$root/modules/gene_calling.nf}")
mkdir "$proof"
export PIXI_NO_INSTALL=true PIXI_FROZEN=true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=1
export NXF_ANSI_LOG=false NXF_OFFLINE=true
export NXF_OPTS='-XX:ActiveProcessorCount=2 -Xms128m -Xmx1g'
root_bash="$root/.pixi/envs/default/bin/bash"
real_pixi=$(command -v pixi)
mkdir -p "$proof/fixture/routing/"{prokaryotic,viral,eukaryotic,unresolved} \
    "$proof/fixture/viral/genomad/input_summary" "$proof/shim"
printf 'seq_name\tlength\ttopology\tcoordinates\n' \
    > "$proof/fixture/viral/genomad/input_summary/input_virus_summary.tsv"
: > "$proof/fixture/viral/genomad/input_summary/input_virus.fna"
# Several long coding regions train single mode; distinct contigs preserve routing exclusivity.
for source in 1 2 3 4; do
    route=viral
    [[ $source != 1 ]] || route=prokaryotic
    awk -v source="$source" 'BEGIN {
        srand(100 + source)
        split("GCT GCC GCA GCG GGT GGC GGA GGG CTT CTC CTG TTA TTT TTC ATT ATC GTG GTT CCG CCT ACC ACT AAT AAC GAT GAC CAG CAA GAA GAG", codons)
        print ">fixture_" source
        for (gene=0; gene<120; gene++) {
            printf "TAATAGTGAAGGAGGAAAAAATG"
            for (i=0; i<300; i++) printf "%s", codons[1 + int(rand()*30)]
            print "TAA"
        }
    }' > "$proof/fixture/routing/$route/$source.fna"
done
# Each 30-base block has TAA at offsets 0,5,10 and reverse-strand TAA
# (forward TTA) at 15,20,25, stopping all six frames. Poly-A permits partial CDS.
awk 'BEGIN {print ">zero_cds"; block="TAACCTAACCTAACCTTACCTTACCTTACC";
    for (i=0;i<2000;i++) printf "%s", substr(block,i%30+1,1); print ""}' \
    > "$proof/fixture/routing/viral/5.fna"

# Intercept only the existing Pixi boundary. Every successful scientific call stays native.
cat > "$proof/shim/pixi" <<'SHIM'
#!/usr/bin/env bash
set -euo pipefail
args=("$@")
input= output=
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) input=$2; shift ;;
        -o) output=$2; shift ;;
    esac
    shift
done
if [[ -z $input ]]; then
    exec "$REAL_PIXI" "${args[@]}"
fi
source=${input##*/}
source=${source%.fna}
printf 'start\t%s\n' "$source" >> "$EVENTS"
trap 'printf "end\t%s\n" "$source" >> "$EVENTS"' EXIT
if [[ $FAULT == nonzero && $source == prokaryotic-bin-1 ]]; then
    exit 23
fi
if [[ $FAULT == nonzero && $source == viral-bin-2 ]]; then
    sleep 4
fi
"$REAL_PIXI" "${args[@]}"
if [[ $FAULT == invalid && $source == viral-bin-5 ]]; then
    : > "$output"
fi
SHIM
chmod +x "$proof/shim/pixi"
cat > "$proof/cli" <<CLI
#!$root_bash
set -euo pipefail
exec "$root/.pixi/envs/default/bin/python" -m protist_meta.cli "\$@"
CLI
chmod +x "$proof/cli"

run_case() {
    local name=$1 cpus=$2 code=$3 fault=$4 module=$5 expected=$6
    local run="$proof/$name"
    local fixture=${7:-$proof/fixture} count=${8:-5}
    mkdir "$run"
    : > "$run/events.tsv"
    cat > "$run/main.nf" <<NF
nextflow.enable.dsl = 2
include { PRODIGAL_GV_GENES } from '$module'
workflow {
    PRODIGAL_GV_GENES(Channel.of(tuple([sample_id: 'CHECK', genetic_code: '$code'],
        file('$fixture/routing'), file('$fixture/viral'))))
}
NF
    cat > "$run/run.config" <<CONFIG
params.cli = '$proof/cli'
params.gene_manifest = '$root/envs/gene/pixi.toml'
params.outdir = '$run/results'
process {
    executor = 'local'
    shell = ['$root_bash', '-euo', 'pipefail']
    errorStrategy = 'terminate'
    withName: PRODIGAL_GV_GENES {
        cpus = $cpus
        memory = '4 GB'
        beforeScript = 'export PATH=$proof/shim:\$PATH REAL_PIXI=$real_pixi EVENTS=$run/events.tsv FAULT=$fault PYTHONPATH=$root/src'
    }
}
executor.cpus = $cpus
executor.memory = '8 GB'
executor.queueSize = 1
trace {
    enabled = true
    file = '$run/trace.tsv'
    fields = 'task_id,process,status,exit,realtime,cpus'
    raw = true
}
CONFIG
    local status=0
    (
        cd "$run"
        "$root/.pixi/envs/default/bin/nextflow" -log "$run/nextflow.log" \
            -C "$run/run.config" run "$run/main.nf" -offline -work-dir "$run/work"
    ) > "$run/stdout" 2>&1 || status=$?
    printf '%s\n' "$status" > "$run/exit-status.txt"
    if [[ $expected == pass ]]; then
        [[ $status -eq 0 ]]
        [[ -s "$run/results/CHECK/12_prodigal_gv/stage.json" ]]
        awk -F '\t' -v cpus="$cpus" 'NR==2 {ok=($3=="COMPLETED" && $4==0 && $6==cpus)} END {exit !ok}' "$run/trace.tsv"
    else
        [[ $status -ne 0 ]]
        [[ -z $(find "$run/work" -name stage.json -print -quit) ]]
        local failed_source=viral-bin-5
        [[ $fault != nonzero ]] || failed_source=prokaryotic-bin-1
        # Inspect the task stderr itself, not Nextflow's repeated error summary.
        local stderr
        stderr=$(find "$run/work" -name .command.err -print)
        [[ -f $stderr ]]
        awk -v expected="$failed_source" '
            /^Prodigal-gv source failed:/ {count++; if ($0 != "Prodigal-gv source failed: " expected) bad=1}
            END {exit (count!=1 || bad)}' "$stderr"
        awk -F '\t' 'NR==2 {ok=($3=="FAILED" && $4!=0)} END {exit !ok}' "$run/trace.tsv"
    fi
    # Observe command overlap and completion, not parent wait calls: Nextflow's
    # tee wrapper can also wait for a child that still holds its output pipes.
    local launched=$count peak=$cpus
    if [[ $fault == nonzero ]]; then launched=2; peak=0; fi
    awk -F '\t' -v expected="$launched" -v cap="$cpus" -v peak="$peak" '
        $1=="start" {if (started[$2]++) bad=1; active++; count++; if(active>max) max=active}
        $1=="end" {if (!started[$2] || ended[$2]++) bad=1; active--}
        active<0 || active>cap {bad=1}
        END {print "launched=" count, "peak=" max, "remaining=" active;
             exit (bad || active!=0 || count!=expected || (peak && max!=peak))}
    ' "$run/events.tsv" > "$run/concurrency.txt"
    if [[ $fault == nonzero ]]; then
        awk -F '\t' '$1=="start" {seen[$2]=1}
            END {exit !(seen["prokaryotic-bin-1"] && seen["viral-bin-2"])}' "$run/events.tsv"
    fi
}

run_case serial 1 auto none "$reference" pass
run_case odd 2 auto none "$root/modules/gene_calling.nf" pass
diff -r "$proof/serial/results" "$proof/odd/results" > "$proof/odd/equivalence.diff"
cp --no-preserve=mode -r "$proof/fixture" "$proof/four"
rm "$proof/four/routing/viral/5.fna"
run_case exact 2 auto none "$root/modules/gene_calling.nf" pass "$proof/four" 4
[[ ! -s "$proof/odd/results/CHECK/12_prodigal_gv/proteins/viral-bin-5.faa" ]]
awk '!/^#/ && NF {count++} END {exit (count!=0)}' \
    "$proof/odd/results/CHECK/12_prodigal_gv/gff/viral-bin-5.gff"
[[ -s "$proof/odd/results/CHECK/12_prodigal_gv/proteins/prokaryotic-bin-1.faa" ]]
run_case code11-serial 1 11 none "$reference" pass
run_case code11-odd 2 11 none "$root/modules/gene_calling.nf" pass
diff -r "$proof/code11-serial/results" "$proof/code11-odd/results" > "$proof/code11-odd/equivalence.diff"
grep -q '^# Model Data:.*run_type=Single;.*transl_table=11;' \
    "$proof/code11-odd/results/CHECK/12_prodigal_gv/gff/prokaryotic-bin-1.gff"
grep -q '^# Model Data:.*run_type=Metagenomic;' \
    "$proof/code11-odd/results/CHECK/12_prodigal_gv/gff/viral-bin-2.gff"
run_case fail-full 2 auto nonzero "$root/modules/gene_calling.nf" fail
run_case fail-tail 2 auto invalid "$root/modules/gene_calling.nf" fail
sha256sum "$reference" "$root/modules/gene_calling.nf" > "$proof/modules.sha256"
find "$proof/fixture" -type f -print0 | sort -z | xargs -0 sha256sum > "$proof/fixture.sha256"
printf 'PASS\n' > "$proof/status.txt"
