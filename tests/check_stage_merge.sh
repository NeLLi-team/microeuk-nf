#!/usr/bin/env bash
# Run the split screen, prokaryote QC, contig and eukaryote QC processes against a
# shimmed `pixi` and check the merged stage directories. Run from the repository root.
set -euo pipefail
root=$PWD
proof=$(realpath -m "${1:?Usage: check_stage_merge.sh NEW_PROOF_DIR}")
mkdir "$proof"
export PIXI_NO_INSTALL=true PIXI_FROZEN=true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NXF_ANSI_LOG=false NXF_OFFLINE=true
export NXF_OPTS='-XX:ActiveProcessorCount=1 -Xms128m -Xmx1g'
root_bash="$root/.pixi/envs/default/bin/bash"
real_python="$root/.pixi/envs/default/bin/python"
real_pixi=$(command -v pixi)
commit=0123456789abcdef0123456789abcdef01234567

fail() {
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

has() {
    grep -Fq -- "$2" "$1" || fail "$1 lacks: $2"
}

lacks() {
    ! grep -Fq -- "$2" "$1" || fail "$1 unexpectedly has: $2"
}

# Fixtures: bins, an assembly and routed bins, each with an empty twin.
fixture="$proof/fixture"
mkdir -p "$fixture/bins/03_binning/bins" "$fixture/nobins/03_binning/bins" \
    "$fixture/assembly/02_assembly" "$fixture/db/gvclass" "$fixture/db/busco" \
    "$fixture/routed/07_routing/"{prokaryotic,eukaryotic,viral,unresolved} \
    "$fixture/norouted/07_routing/"{prokaryotic,eukaryotic,viral,unresolved}
for bin in bin_1 bin_2; do
    printf '>%s_contig\nACGTACGTACGT\n' "$bin" > "$fixture/bins/03_binning/bins/$bin.fa"
done
printf '>contig_1\nACGTACGTACGT\n' > "$fixture/assembly/02_assembly/assembly_primary.fa"
for bin in p1 p2; do
    printf '>%s_contig\nACGTACGTACGT\n' "$bin" > "$fixture/routed/07_routing/prokaryotic/$bin.fna"
done
for bin in e1 e2; do
    printf '>%s_contig\nACGTACGTACGT\n' "$bin" > "$fixture/routed/07_routing/eukaryotic/$bin.fna"
done

# The only fake is the Pixi boundary; every other command in the process scripts runs natively.
mkdir -p "$proof/shim"
cat > "$proof/shim/pixi" <<'SHIM'
#!/usr/bin/env bash
set -euo pipefail
shift
environment=default
while [[ $# -gt 0 ]]; do
    case "$1" in
        --as-is | --quiet) shift ;;
        --manifest-path) shift 2 ;;
        --environment) environment=$2; shift 2 ;;
        *) break ;;
    esac
done
command_name=$1
shift
printf '%s\t%s %s\n' "$environment" "$command_name" "$*" >> "$SHIM_LOG"
option() {
    local key=$1
    shift
    while [[ $# -gt 0 ]]; do
        if [[ $1 == "$key" ]]; then
            printf '%s\n' "$2"
            return
        fi
        shift
    done
}
case "$command_name $*" in
    "checkeuk --version") echo 'checkeuk 9.9.9' ;;
    "checkeuk "*)
        out="$(option -o "$@")_checkeuk"
        mkdir -p "$out"
        printf 'genome\tstatus\n' > "$out/checkeuk_report.tsv"
        for fasta in "$(option -i "$@")"/*.fa; do
            printf '%s\tok\n' "$(basename "$fasta" .fa)" >> "$out/checkeuk_report.tsv"
        done
        printf 'profile\n' > "$(option --profile-out "$@")" ;;
    "checkm2 --version") echo '1.1.0' ;;
    "checkm2 predict "*)
        out=$(option --output-directory "$@")
        mkdir -p "$out/protein_files" "$out/diamond_output"
        echo 'checkm2 predict (shim)'
        for fasta in "$(option --input "$@")"/*.fa; do
            printf '>p\nMKV\n' > "$out/protein_files/$(basename "$fasta" .fa).faa"
        done
        if [[ ${SHIM_CHECKM2:-ok} == nodiamond ]]; then
            printf '[2026-01-01 00:00:00] ERROR: No DIAMOND annotation was generated. Exiting\n' \
                > "$out/checkm2.log"
            : > "$out/diamond_output/DIAMOND_RESULTS.tsv"
            exit 1
        fi
        printf 'Name\tCompleteness\tContamination\n' > "$out/quality_report.tsv"
        printf 'INFO\n' > "$out/checkm2.log" ;;
    "python -m protist_meta.checkm2 "*) ;;
    "python "*"check-gvclass-reference.py "*) exec "$REAL_PYTHON" "$@" ;;
    "python -c "*importlib.metadata*) echo '1.2.3' ;;
    "python -c "*checkm*) echo '1.2.10' ;;
    "gvclass --version")
        printf 'GVClass Pipeline\n  Software version: v9.9.9\n  Database version: not installed\n' ;;
    "gvclass "*)
        out=$(option -o "$@")
        mkdir -p "$out"
        printf 'query\tdomain\n' > "$out/gvclass_summary.tsv"
        for fasta in "$1"/*.fna; do
            printf '%s\td_EUK\n' "$(basename "$fasta" .fna)" >> "$out/gvclass_summary.tsv"
        done
        printf '{"database": {"path": "%s", "version": "%s"}}\n' \
            "$(option --database "$@")" "$SHIM_GVCLASS_DB_VERSION" > "$out/run_status.json" ;;
    "checkm lineage_wf "*)
        printf 'Bin Id\tMarker lineage\n' > "$(option --file "$@")"
        mkdir -p "${*: -1}" ;;
    "gtdbtk classify_wf "*)
        out=$(option --out_dir "$@")
        mkdir -p "$out/classify"
        printf 'user_genome\tclassification\n' > "$out/gtdbtk.bac120.summary.tsv" ;;
    "symclatron --version") echo 'symclatron version 9.9.9' ;;
    "symclatron classify "*)
        out=$(option --output-dir "$@")
        mkdir -p "$out"
        printf 'taxon_oid\tclassification\n' > "$out/symclatron_results.tsv"
        for fasta in "$(option --genome-dir "$@")"/*.fna; do
            printf '%s\tx\n' "$(basename "$fasta" .fna)" >> "$out/symclatron_results.tsv"
        done ;;
    "busco --version") echo 'BUSCO 6.1.0' ;;
    "busco "*)
        out="$(option --out_path "$@")/$(option -o "$@")"
        mkdir -p "$out"
        printf 'Input_file\tDataset\tComplete\n' > "$out/batch_summary.txt"
        for fasta in "$(option -i "$@")"/*.fna; do
            mkdir -p "$out/$(basename "$fasta")"
            printf '%s\teukaryota_odb12.2\t90.0\n' "$(basename "$fasta")" >> "$out/batch_summary.txt"
        done ;;
    "tidk --version")
        echo "Warning! No clades found in the database." >&2
        echo 'tidk 9.9.9' ;;
    "tidk explore "*)
        printf 'canonical_repeat_unit\tcopies\tcopies_as_unit\tcopies_as_revcomp\n'
        if [[ ${*: -1} == */e1.fna ]]; then
            printf 'AACCCT\t399\t199\t200\nAAACCC\t399\t199\t200\nCCCTAAT\t12\t6\t6\n'
        fi ;;
    "tiara "*)
        out=$(option -o "$@")
        printf 'sequence_id\tclass_fst_stage\tclass_snd_stage\n' > "$out"
        printf 'contig_1\teukarya\tn/a\n' >> "$out"
        printf 'log\n' > "$(dirname "$out")/log_$(basename "$out")" ;;
    "genomad --version") echo 'geNomad, version 9.9.9' ;;
    "genomad end-to-end "*)
        out="${*: -2:1}"
        mkdir -p "$out/input_summary"
        printf 'seq_name\tlength\n' > "$out/input_summary/input_virus_summary.tsv"
        printf '>v\nACGT\n' > "$out/input_summary/input_virus.fna" ;;
    "checkv end_to_end "*)
        mkdir -p "$3"
        printf 'contig_id\n' > "$3/quality_summary.tsv" ;;
    *) echo "unexpected shimmed command: $command_name $*" >&2; exit 97 ;;
esac
SHIM
chmod +x "$proof/shim/pixi"

# The harness repeats the wiring of main.nf from the binning, assembly and routing stages on.
mkdir -p "$proof/project"
ln -s "$root/scripts" "$proof/project/scripts"
cat > "$proof/project/main.nf" <<NF
include {
    CHECKEUK_SCREEN; CHECKM2_SCREEN; GVCLASS_SCREEN; ALL_BIN_SCREEN; CONTIG_CLASSIFY;
    VIRAL_SCREEN; CHECKM1_QC; GTDBTK_CLASSIFY; SYMCLATRON_CLASSIFY;
    PROKARYOTE_CHARACTERIZATION; EUKARYOTE_QC
} from '$root/modules/characterization'

workflow {
    def meta = [sample_id: 'T1']
    def binning = channel.of(tuple(meta, file(params.binning)))
    def assembly = channel.of(tuple(meta, file(params.assembly)))
    def routing = channel.of(tuple(meta, file(params.routing)))

    CHECKEUK_SCREEN(binning)
    CHECKM2_SCREEN(binning)
    GVCLASS_SCREEN(binning)
    def screenInput = binning
        .map { m, binningDir -> tuple(m.sample_id, m, binningDir) }
        .join(CHECKEUK_SCREEN.out.result.map { m, dir, report -> tuple(m.sample_id, dir, report) })
        .join(CHECKM2_SCREEN.out.result.map { m, dir, report -> tuple(m.sample_id, dir, report) })
        .join(GVCLASS_SCREEN.out.result.map { m, dir, report -> tuple(m.sample_id, dir, report) })
        .map { _id, m, binningDir, euDir, euReport, c2Dir, c2Report, gvDir, gvReport ->
            tuple(m, binningDir, euDir, euReport, c2Dir, c2Report, gvDir, gvReport)
        }
    ALL_BIN_SCREEN(screenInput)
    CONTIG_CLASSIFY(assembly)
    VIRAL_SCREEN(assembly)

    CHECKM1_QC(routing)
    GTDBTK_CLASSIFY(routing)
    SYMCLATRON_CLASSIFY(routing)
    def prokaryoteInput = routing
        .map { m, routingDir -> tuple(m.sample_id, m, routingDir) }
        .join(CHECKM1_QC.out.result.map { m, dir, report -> tuple(m.sample_id, dir, report) })
        .join(GTDBTK_CLASSIFY.out.result.map { m, dir, report -> tuple(m.sample_id, dir, report) })
        .join(SYMCLATRON_CLASSIFY.out.result.map { m, dir, report -> tuple(m.sample_id, dir, report) })
        .map { _id, m, routingDir, c1Dir, c1Report, gtDir, gtReport, syDir, syReport ->
            tuple(m, routingDir, c1Dir, c1Report, gtDir, gtReport, syDir, syReport)
        }
    PROKARYOTE_CHARACTERIZATION(prokaryoteInput)
    EUKARYOTE_QC(routing)
}
NF
cat > "$proof/project/test.config" <<CONFIG
process.shell = ['$root_bash', '-euo', 'pipefail']
process.executor = 'local'
process.maxForks = 1
executor.cpus = 32
executor.memory = '360 GB'
process.beforeScript = 'export PATH=$proof/shim:\$PATH'
params.core_manifest = '$root/pixi.toml'
params.checkeuk_app = '$fixture/db'
params.checkeuk_app_commit = '$commit'
params.checkeuk_db = '$fixture/db'
params.checkeuk_db_version = 'v-checkeuk'
params.gvclass_app = '$fixture/db'
params.gvclass_app_commit = '$commit'
params.gvclass_db = '$fixture/db/gvclass'
params.gvclass_db_version = 'v-gvclass'
params.checkm_manifest = '$fixture/db'
params.checkm1_db = '$fixture/db'
params.checkm1_db_version = 'v-checkm1'
params.checkm2_db = '$fixture/db'
params.checkm2_db_version = 'v-checkm2'
params.gtdbtk_db = '$fixture/db'
params.gtdbtk_db_version = 'v-gtdbtk'
params.viral_manifest = '$fixture/db'
params.checkv_manifest = '$fixture/db'
params.genomad_db = '$fixture/db'
params.genomad_db_version = 'v-genomad'
params.checkv_db = '$fixture/db'
params.checkv_db_version = 'v-checkv'
params.contig_manifest = '$fixture/db'
params.tiara_min_length = 1000
params.eukqc_manifest = '$fixture/db'
params.busco_db = '$fixture/db/busco'
params.busco_db_version = 'v-busco'
CONFIG

run_case() {
    local name=$1 binning=$2 routing=$3
    local run="$proof/$name"
    mkdir -p "$run"
    : > "$run/shim.log"
    (
        cd "$run"
        export SHIM_LOG="$run/shim.log" SHIM_CHECKM2=${SHIM_CHECKM2:-ok}
        export REAL_PYTHON="$real_python" SHIM_GVCLASS_DB_VERSION=v-gvclass
        "$real_pixi" run --as-is --manifest-path "$root/pixi.toml" \
            nextflow -log "$run/nextflow.log" \
            -C "$root/conf/resources.config,$proof/project/test.config" \
            run "$proof/project/main.nf" -offline -work-dir "$run/work" \
            --outdir "$run/results" --binning "$fixture/$binning/03_binning" \
            --assembly "$fixture/assembly/02_assembly" --routing "$fixture/$routing/07_routing"
    ) > "$run/stdout" 2>&1 || { cat "$run/stdout" >&2; fail "workflow case $name"; }
    # Rendered scripts are valid Bash, and merged directories hold real files only.
    for script in "$run"/work/*/*/.command.sh; do
        bash -n "$script" || fail "bash -n $script"
    done
    [[ -z $(find "$run/results" -type l -print -quit) ]] || fail "$name published a symlink"
    for stage in "$run"/results/T1/*/stage.json; do
        "$real_python" -m json.tool "$stage" > /dev/null || fail "invalid JSON: $stage"
    done
}

screen_entry() {
    # name version status reason output
    printf '"%s": {"version": "%s", "status": "%s", "reason": %s, "database_name": "%s_db", ' \
        "$1" "$2" "$3" "$4" "$1"
    printf '"database_version": "v-%s", "output": "%s", "command_file": "%s/command.sh"}' \
        "$1" "$5" "$1"
}

# Case 1: bins present, every tool completes.
run_case bins bins routed
screen=$proof/bins/results/T1/04_all_bin_screen
has "$screen/stage.json" '"status": "completed"'
has "$screen/stage.json" '"reason": null,'
has "$screen/stage.json" "\"tool_version\": \"CheckEUK checkeuk 9.9.9 (${commit:0:12}); CheckM2 1.1.0; GVClass 9.9.9 (${commit:0:12})\""
has "$screen/stage.json" "$(screen_entry checkeuk "checkeuk 9.9.9 (${commit:0:12})" completed null checkeuk/result_checkeuk/checkeuk_report.tsv)"
has "$screen/stage.json" "$(screen_entry checkm2 1.1.0 completed null checkm2/output/quality_report.tsv)"
has "$screen/stage.json" "$(screen_entry gvclass "9.9.9 (${commit:0:12})" completed null gvclass/output/gvclass_summary.tsv)"
has "$screen/stage.json" '"outputs": {'
for path in checkeuk/result_checkeuk/checkeuk_report.tsv checkeuk/profile.tsv \
    checkm2/output/quality_report.tsv checkm2/command.log checkm2/exit_code.txt \
    gvclass/output/gvclass_summary.tsv gvclass/output/run_status.json \
    gvclass/input/bin_1.fna gvclass/input/bin_2.fna \
    checkeuk/command.sh checkm2/command.sh gvclass/command.sh; do
    [[ -s $screen/$path ]] || fail "missing $path"
done
has "$screen/checkeuk/command.sh" '--checkpoint-every 10'
[[ $(stat -c %h "$screen/checkm2/output/quality_report.tsv") -ge 2 ]] || fail 'not hard linked'
has "$proof/bins/shim.log" '-t 20 -j 20 --threads-per-worker 1 --min-length 0 --plain-output'
has "$proof/bins/shim.log" 'genomad end-to-end --cleanup --splits 2 --threads 12'
has "$proof/bins/shim.log" 'checkv end_to_end'
has "$proof/bins/shim.log" 'tiara -i'
tiara=$proof/bins/results/T1/15_contig_classification
has "$tiara/tiara.tsv" 'sequence_id'
[[ -f $tiara/log_tiara.tsv ]] || fail 'missing Tiara log'
has "$tiara/stage.json" '"stage": "tiara"'
has "$tiara/stage.json" '"outputs": {'
has "$proof/bins/results/T1/06_viral/stage.json" '"checkv": {"status": "completed"'

# Case 2: no bins; all three tools skipped and the stage keeps today's semantics.
run_case nobins nobins norouted
screen=$proof/nobins/results/T1/04_all_bin_screen
has "$screen/stage.json" '"status": "skipped"'
has "$screen/stage.json" '"reason": "no_bins",'
has "$screen/stage.json" '"tool_version": "not_run"'
has "$screen/stage.json" '"command": "No bin-screen tools run because QuickBin produced no accepted bins."'
for tool in checkeuk checkm2 gvclass; do
    has "$screen/stage.json" "$(screen_entry "$tool" not_run skipped '"no_bins"' "$tool")"
    [[ -s $screen/$tool/command.sh ]] || fail "missing $tool command"
done
has "$screen/stage.json" '"checkeuk": "checkeuk",'
has "$screen/stage.json" '"gvclass": "gvclass"'

# Case 3: CheckM2 finds no DIAMOND annotations; only that tool is skipped.
SHIM_CHECKM2=nodiamond run_case nodiamond bins norouted
screen=$proof/nodiamond/results/T1/04_all_bin_screen
has "$screen/stage.json" '"status": "completed"'
has "$screen/stage.json" "$(screen_entry checkm2 1.1.0 skipped '"no_diamond_annotations"' checkm2)"
has "$screen/stage.json" '"checkm2": "checkm2",'
has "$screen/checkm2/exit_code.txt" '1'
[[ -s $screen/checkm2/output/checkm2.log ]] || fail 'missing CheckM2 log'

# Case 4: routed bins; prokaryote and eukaryote QC complete.
prok=$proof/bins/results/T1/08_prokaryote_qc
has "$prok/stage.json" '"status": "completed"'
has "$prok/stage.json" '"reason": null,'
has "$prok/stage.json" '"tool_version": "CheckM1 1.2.10; GTDB-Tk 1.2.3; Symclatron symclatron version 9.9.9"'
has "$prok/stage.json" "$(screen_entry checkm1 1.2.10 completed null checkm1/quality_report.tsv)"
has "$prok/stage.json" "$(screen_entry gtdbtk 1.2.3 completed null gtdbtk)"
has "$prok/stage.json" '"symclatron": {"version": "symclatron version 9.9.9", "status": "completed", "reason": null, "database_name": "symclatron_bundled", "database_version": "bundled", "output": "symclatron/symclatron_results.tsv", "command_file": "symclatron/command.sh"}'
for path in checkm1/quality_report.tsv gtdbtk/gtdbtk.bac120.summary.tsv \
    symclatron/symclatron_results.tsv checkm1/command.sh gtdbtk/command.sh symclatron/command.sh; do
    [[ -s $prok/$path ]] || fail "missing $path"
done
[[ $(awk 'END {print NR - 1}' "$prok/symclatron/symclatron_results.tsv") -eq 2 ]] || fail 'symclatron rows'
euk=$proof/bins/results/T1/16_eukaryote_qc
has "$euk/stage.json" '"status": "completed"'
has "$euk/stage.json" '"tool_version": "BUSCO 6.1.0; tidk 9.9.9"'
has "$euk/stage.json" '"busco": {"version": "BUSCO 6.1.0", "status": "completed", "reason": null, "database_name": "busco_db", "database_version": "v-busco", "output": "busco/batch_summary.txt"}'
has "$euk/stage.json" '"tidk": {"version": "tidk 9.9.9", "status": "completed", "reason": null, "database_name": null, "database_version": null, "output": "tidk/tidk_summary.tsv"}'
[[ $(awk 'END {print NR - 1}' "$euk/busco/batch_summary.txt") -eq 2 ]] || fail 'busco rows'
[[ -s $euk/tidk/e1.explore.tsv && -s $euk/tidk/e2.explore.tsv ]] || fail 'missing tidk tables'
printf 'bin_id\ttop_repeat\ttop_repeat_count\ne1\tAAACCC\t399\ne2\t\t\n' > "$proof/tidk-expected.tsv"
cmp "$proof/tidk-expected.tsv" "$euk/tidk/tidk_summary.tsv" || fail 'tidk summary'
has "$proof/bins/shim.log" 'busco -i'
has "$proof/bins/shim.log" '-m genome -l '
has "$proof/bins/shim.log" '--offline --opt-out-run-stats -c 16'
has "$proof/bins/shim.log" 'tidk explore --minimum 5 --maximum 12 --threshold 5'

# Case 5: no routed bins; both stages skip with their reasons.
prok=$proof/nobins/results/T1/08_prokaryote_qc
has "$prok/stage.json" '"status": "skipped"'
has "$prok/stage.json" '"reason": "no_routed_prokaryotic_bins",'
has "$prok/stage.json" "$(screen_entry checkm1 1.2.10 skipped '"no_routed_prokaryotic_bins"' checkm1/quality_report.tsv)"
euk=$proof/nobins/results/T1/16_eukaryote_qc
has "$euk/stage.json" '"status": "skipped"'
has "$euk/stage.json" '"reason": "no_routed_eukaryotic_bins",'
[[ -d $euk/busco && -d $euk/tidk ]] || fail 'empty eukaryote QC layout'
lacks "$proof/nobins/shim.log" 'busco -i'

# Case 6: tool-entry.sh refuses a value that would break the JSON it prints.
entry_dir=$proof/tool-entry
mkdir -p "$entry_dir"
printf '1.0\tcompleted\tnull\tout.tsv\n' > "$entry_dir/good.tsv"
printf '1.0\tcompleted\tbad"reason\tout.tsv\n' > "$entry_dir/quote.tsv"
printf '1.0\tcompleted\tbad\\reason\tout.tsv\n' > "$entry_dir/backslash.tsv"
bash "$root/scripts/tool-entry.sh" tool "$entry_dir/good.tsv" db v1 > /dev/null \
    || fail 'tool-entry.sh rejects a clean report'
for report in quote backslash; do
    ! bash "$root/scripts/tool-entry.sh" tool "$entry_dir/$report.tsv" db v1 \
        > /dev/null 2> "$entry_dir/$report.err" || fail "tool-entry.sh accepts $report"
    has "$entry_dir/$report.err" 'must not contain'
done
! bash "$root/scripts/tool-entry.sh" tool "$entry_dir/good.tsv" 'db"x' v1 \
    > /dev/null 2> /dev/null || fail 'tool-entry.sh accepts a quote in a database name'

printf 'PASS\n' > "$proof/status.txt"
printf 'PASS stage merge checks\n'
