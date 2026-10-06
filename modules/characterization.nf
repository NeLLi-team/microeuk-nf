import groovy.json.JsonOutput

process ALL_BIN_SCREEN {
    tag "${meta.sample_id}"
    label 'heavy_qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(binning_dir)

    output:
    tuple val(meta), path('04_all_bin_screen'), emit: result

    script:
    """
    mkdir -p 04_all_bin_screen/checkeuk 04_all_bin_screen/checkm2 \
        04_all_bin_screen/gvclass/input 04_all_bin_screen/gvclass/output
    BIN_COUNT=\$(find "${binning_dir}/bins" -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l)
    if [[ \${BIN_COUNT} -eq 0 ]]; then
        STATUS=skipped
        REASON='"no_bins"'
        TOOL_VERSION=not_run
        TOOL_STATUS=skipped
        TOOL_REASON='"no_bins"'
        CHECKM2_STATUS=skipped
        CHECKM2_REASON='"no_bins"'
        CHECKEUK_VERSION=not_run
        CHECKM2_VERSION=not_run
        GVCLASS_VERSION=not_run
        CHECKEUK_OUTPUT=checkeuk
        CHECKM2_OUTPUT=checkm2
        GVCLASS_OUTPUT=gvclass
        COMMAND='No bin-screen tools run because QuickBin produced no accepted bins.'
    else
        STATUS=completed
        REASON=null
        TOOL_STATUS=completed
        TOOL_REASON=null
        CHECKM2_STATUS=completed
        CHECKM2_REASON=null
        CHECKEUK_OUTPUT=checkeuk/result_checkeuk/checkeuk_report.tsv
        CHECKM2_OUTPUT=checkm2/output/quality_report.tsv
        GVCLASS_OUTPUT=gvclass/output/gvclass_summary.tsv
        COMMAND='all-bin CheckEUK, CheckM2, and GVClass; GVClass physical .fna copies'

        pixi run --as-is --manifest-path "${params.checkeuk_app}/pixi.toml" checkeuk \
            -i "\$PWD/${binning_dir}/bins" -o "\$PWD/04_all_bin_screen/checkeuk/result" \
            --resources "${params.checkeuk_db}" --threads ${task.cpus} --checkpoint-every 10 \
            --profile-out "\$PWD/04_all_bin_screen/checkeuk/profile.tsv"
        CHECKEUK_REPORT=04_all_bin_screen/checkeuk/result_checkeuk/checkeuk_report.tsv
        [[ -s \${CHECKEUK_REPORT} ]]
        [[ \$(awk 'END {print NR - 1}' \${CHECKEUK_REPORT}) -eq \${BIN_COUNT} ]]
        awk -F '\t' 'NR == 1 {for (i=1; i<=NF; i++) if (\$i=="status") status=i; next}
            status && (\$status == "" || \$status == "error") {failed++}
            END {exit status == 0 || failed > 0}' \${CHECKEUK_REPORT}

        CHECKM2_EXIT=0
        pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment checkm2 \
            checkm2 predict --input "${binning_dir}/bins" \
            --output-directory 04_all_bin_screen/checkm2/output --threads ${task.cpus} \
            --database_path "${params.checkm2_db}" --extension fa --tmpdir /tmp --force \
            > 04_all_bin_screen/checkm2/command.log 2>&1 || CHECKM2_EXIT=\$?
        printf '%s\n' "\${CHECKM2_EXIT}" > 04_all_bin_screen/checkm2/exit_code.txt
        cat 04_all_bin_screen/checkm2/command.log
        if [[ \${CHECKM2_EXIT} -eq 0 ]]; then
            [[ -s 04_all_bin_screen/checkm2/output/quality_report.tsv ]]
            [[ \$(awk 'END {print NR - 1}' 04_all_bin_screen/checkm2/output/quality_report.tsv) -eq \${BIN_COUNT} ]]
        elif bash "${projectDir}/scripts/checkm2-no-annotations.sh" "\${CHECKM2_EXIT}" \
            04_all_bin_screen/checkm2/output "${binning_dir}/bins"; then
            CHECKM2_STATUS=skipped
            CHECKM2_REASON='"no_diamond_annotations"'
            CHECKM2_OUTPUT=checkm2
        else
            exit "\${CHECKM2_EXIT}"
        fi

        for FASTA in "${binning_dir}"/bins/*.fa; do
            cp -L "\${FASTA}" "04_all_bin_screen/gvclass/input/\$(basename "\${FASTA}" .fa).fna"
        done
        [[ \$(find 04_all_bin_screen/gvclass/input -type f -name '*.fna' -size +0c | wc -l) -eq \${BIN_COUNT} ]]
        export GVCLASS_RESOURCE_CACHE="\$PWD/tmp/gvclass_cache"
        export GVCLASS_PLAIN_OUTPUT=1
        pixi run --as-is --manifest-path "${params.gvclass_app}/pixi.toml" gvclass \
            "\$PWD/04_all_bin_screen/gvclass/input" -o "\$PWD/04_all_bin_screen/gvclass/output" \
            --database "${params.gvclass_db}" -t ${task.cpus} -j 4 --min-length 0 --plain-output
        pixi run --as-is --quiet --manifest-path "${params.core_manifest}" \
            python "${projectDir}/scripts/check-gvclass-reference.py" \
            04_all_bin_screen/gvclass/output/run_status.json \
            "${params.gvclass_db}" "${params.gvclass_db_version}"
        [[ -s 04_all_bin_screen/gvclass/output/gvclass_summary.tsv ]]
        [[ \$(awk 'END {print NR - 1}' 04_all_bin_screen/gvclass/output/gvclass_summary.tsv) -eq \${BIN_COUNT} ]]
        if [[ -s 04_all_bin_screen/gvclass/output/gvclass_failed_queries.tsv ]]; then
            [[ \$(awk 'END {print NR - 1}' 04_all_bin_screen/gvclass/output/gvclass_failed_queries.tsv) -eq 0 ]]
        fi

        CHECKEUK_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkeuk_app}/pixi.toml" \
            checkeuk --version | tr '\n' ' ' | sed 's/"/\\\\"/g')
        CHECKM2_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
            --environment checkm2 checkm2 --version 2>&1 | tail -n 1 | sed 's/"/\\\\"/g')
        GVCLASS_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.gvclass_app}/pixi.toml" \
            python -c 'import sys; sys.path.insert(0, "${params.gvclass_app}"); from src.__version__ import __version__; print(__version__)')
        TOOL_VERSION="CheckEUK \${CHECKEUK_VERSION}; CheckM2 \${CHECKM2_VERSION}; GVClass \${GVCLASS_VERSION}"
    fi
    cat > 04_all_bin_screen/stage.json <<JSON
    {
      "stage": "all_bin_screen",
      "status": "\${STATUS}",
      "tool": "CheckEUK;CheckM2;GVClass",
      "tool_version": "\${TOOL_VERSION}",
      "database_name": "checkeuk_db;checkm2_db;gvclass_db",
      "database_version": "CheckEUK ${params.checkeuk_db_version}; CheckM2 ${params.checkm2_db_version}; GVClass ${params.gvclass_db_version}",
      "command": "\${COMMAND}",
      "reason": \${REASON},
      "tools": {
        "checkeuk": {"version": "\${CHECKEUK_VERSION}", "status": "\${TOOL_STATUS}", "reason": \${TOOL_REASON}, "database_name": "checkeuk_db", "database_version": "${params.checkeuk_db_version}", "output": "\${CHECKEUK_OUTPUT}"},
        "checkm2": {"version": "\${CHECKM2_VERSION}", "status": "\${CHECKM2_STATUS}", "reason": \${CHECKM2_REASON}, "database_name": "checkm2_db", "database_version": "${params.checkm2_db_version}", "output": "\${CHECKM2_OUTPUT}"},
        "gvclass": {"version": "\${GVCLASS_VERSION}", "status": "\${TOOL_STATUS}", "reason": \${TOOL_REASON}, "database_name": "gvclass_db", "database_version": "${params.gvclass_db_version}", "output": "\${GVCLASS_OUTPUT}"}
      },
      "outputs": {
        "checkeuk": "\${CHECKEUK_OUTPUT}",
        "checkm2": "\${CHECKM2_OUTPUT}",
        "gvclass": "\${GVCLASS_OUTPUT}"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 04_all_bin_screen/checkeuk/result_checkeuk \
        04_all_bin_screen/checkm2/output 04_all_bin_screen/gvclass/output
    printf 'genome\tstatus\n' > 04_all_bin_screen/checkeuk/result_checkeuk/checkeuk_report.tsv
    printf 'Name\tCompleteness\tContamination\n' > 04_all_bin_screen/checkm2/output/quality_report.tsv
    printf 'genome\tdomain\n' > 04_all_bin_screen/gvclass/output/gvclass_summary.tsv
    printf '%s\n' '{"stage":"all_bin_screen","status":"completed","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"checkeuk":"checkeuk/result_checkeuk/checkeuk_report.tsv","checkm2":"checkm2/output/quality_report.tsv","gvclass":"gvclass/output/gvclass_summary.tsv"}}' > 04_all_bin_screen/stage.json
    """
}

process SSU_EXTRACT {
    tag "${meta.sample_id}"
    label 'heavy_qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly_dir)

    output:
    tuple val(meta), path('05_ssu'), emit: result

    script:
    if (task.cpus < 2 || task.memory.toGiga() < 5) {
        throw new IllegalArgumentException('SSU_EXTRACT requires at least 2 CPUs and 5 GB')
    }
    def nestedCpus = task.cpus - 1
    def nestedMemoryGb = task.memory.toGiga() - 4
    """
    TASK_ROOT=\$(pwd -P)
    LAUNCH_ROOT="\${TASK_ROOT}/ssu-launch"
    NESTED_WORK_ROOT="\${TASK_ROOT}/ssu-work"
    OUTPUT_ROOT="\${TASK_ROOT}/05_ssu/output"
    mkdir -p "\${LAUNCH_ROOT}" "\${NESTED_WORK_ROOT}" "\${OUTPUT_ROOT}"
    cat > "\${LAUNCH_ROOT}/nested.config" <<'CONFIG'
    process.shell = ['${projectDir}/.pixi/envs/default/bin/bash', '-euo', 'pipefail']
    executor {
        name = 'local'
        cpus = ${nestedCpus}
        memory = '${nestedMemoryGb} GB'
    }
    process {
        withName: BLAST_ANNOTATE {
            cpus = ${Math.min(8, nestedCpus)}
            time = { [8.h * task.attempt, params.max_time as nextflow.util.Duration].min() }
            maxRetries = 1
            errorStrategy = {
                def error = task.previousException
                def cause = error?.cause ?: error
                (cause instanceof nextflow.exception.ProcessException &&
                    cause.message?.startsWith('Process exceeded running time limit')) ||
                    task.exitStatus in ((130..145) + 104)
                    ? 'retry'
                    : 'finish'
            }
        }
    }
CONFIG
    SSU_DB_VERSION=\$(pixi run --as-is --manifest-path "${params.ssuextract_app}/pixi.toml" \
        python "${params.ssuextract_app}/scripts/database_manager.py" version \
        --root "${params.ssuextract_db}" --profile curated)
    [[ "\${SSU_DB_VERSION}" == "${params.ssuextract_db_version}" ]]
    export NXF_OFFLINE=true
    export NXF_OPTS='-XX:ActiveProcessorCount=1 -Xms256m -Xmx2g'
    cd "\${LAUNCH_ROOT}"
    pixi run --as-is --manifest-path "${params.ssuextract_app}/pixi.toml" \
        nextflow -c "\${LAUNCH_ROOT}/nested.config" \
        run "${params.ssuextract_app}/main.nf" -offline \
        -profile local -work-dir "\${NESTED_WORK_ROOT}" \
        --query "\${TASK_ROOT}/${assembly_dir}/assembly_primary.fa" \
        --modeldir "${params.ssuextract_models}" --outdir "\${OUTPUT_ROOT}" \
        --database_path "${params.ssuextract_db}" --database_profile curated \
        --threads_per_job 4 --max_cpus ${nestedCpus} \
        --max_memory '${nestedMemoryGb}.GB' --max_time ${task.time.toHours()}.h
    cd "\${TASK_ROOT}"
    [[ \$(pwd -P) == "\${TASK_ROOT}" ]]
    [[ -s 05_ssu/output/cmsearch_summary.tsv ]]
    [[ -s 05_ssu/output/cmsearch_summary.tab ]]
    [[ -s 05_ssu/output/blast_top_hits.tsv ]]
    cat > 05_ssu/stage.json <<JSON
    {
      "stage": "ssu",
      "status": "completed",
      "tool": "SSUextract",
      "tool_version": "${params.ssuextract_app_version}",
      "database_name": "ssuextract_db",
      "database_version": "${params.ssuextract_db_version}",
      "command": "SSUextract local nested workflow with bounded resources",
      "outputs": {
        "cmsearch_summary": "output/cmsearch_summary.tsv",
        "cmsearch_table": "output/cmsearch_summary.tab",
        "blast_hits": "output/blast_top_hits.tsv",
        "extracted": "output/extracted",
        "stats": "output/stats"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 05_ssu/output/extracted 05_ssu/output/stats
    printf 'sequence\tmodel\n' > 05_ssu/output/cmsearch_summary.tsv
    printf '# stub\n' > 05_ssu/output/cmsearch_summary.tab
    printf 'query\tsubject\n' > 05_ssu/output/blast_top_hits.tsv
    printf '%s\n' '{"stage":"ssu","status":"completed","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"cmsearch_summary":"output/cmsearch_summary.tsv","cmsearch_table":"output/cmsearch_summary.tab","blast_hits":"output/blast_top_hits.tsv","extracted":"output/extracted","stats":"output/stats"}}' > 05_ssu/stage.json
    """
}

process VIRAL_SCREEN {
    tag "${meta.sample_id}"
    label 'heavy_qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly_dir)

    output:
    tuple val(meta), path('06_viral'), emit: result

    script:
    def checkvManifest = params.checkv_manifest ?: params.viral_manifest
    """
    mkdir -p 06_viral/genomad 06_viral/checkv
    cp "${assembly_dir}/assembly_primary.fa" input.fna
    pixi run --as-is --manifest-path "${params.viral_manifest}" \
        genomad end-to-end --cleanup --splits 2 --threads ${task.cpus} \
        input.fna 06_viral/genomad "${params.genomad_db}"
    [[ -s 06_viral/genomad/input_summary/input_virus_summary.tsv ]]

    VIRUS_FASTA=06_viral/genomad/input_summary/input_virus.fna
    CHECKV_STATUS=skipped
    if [[ -s \${VIRUS_FASTA} ]]; then
        pixi run --as-is --manifest-path "${checkvManifest}" \
            checkv end_to_end "\${VIRUS_FASTA}" 06_viral/checkv \
            -t ${task.cpus} -d "${params.checkv_db}"
        [[ -s 06_viral/checkv/quality_summary.tsv ]]
        CHECKV_STATUS=completed
    fi
    GENOMAD_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.viral_manifest}" \
        genomad --version | tr '\n' ' ' | sed 's/"/\\\\"/g')
    CHECKV_VERSION=\$(pixi run --as-is --quiet --manifest-path "${checkvManifest}" \
        python -c 'from importlib.metadata import version; print(version("checkv"))')
    cat > 06_viral/stage.json <<JSON
    {
      "stage": "viral",
      "status": "completed",
      "tool": "geNomad;CheckV",
      "tool_version": "geNomad \${GENOMAD_VERSION}; CheckV \${CHECKV_VERSION}",
      "database_name": "genomad_db;checkv_db",
      "database_version": "geNomad ${params.genomad_db_version}; CheckV ${params.checkv_db_version}",
      "command": "genomad end-to-end; CheckV on nonempty geNomad viral contigs",
      "tools": {
        "genomad": {"status": "completed", "version": "\${GENOMAD_VERSION}", "database_name": "genomad_db", "database_version": "${params.genomad_db_version}", "output": "genomad/input_summary"},
        "checkv": {"status": "\${CHECKV_STATUS}", "version": "\${CHECKV_VERSION}", "database_name": "checkv_db", "database_version": "${params.checkv_db_version}", "output": "checkv"}
      },
      "outputs": {
        "genomad": "genomad/input_summary",
        "viral_fasta": "genomad/input_summary/input_virus.fna",
        "checkv": "checkv"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 06_viral/genomad/input_summary 06_viral/checkv
    printf 'seq_name\tlength\n' > 06_viral/genomad/input_summary/input_virus_summary.tsv
    : > 06_viral/genomad/input_summary/input_virus.fna
    printf '%s\n' '{"stage":"viral","status":"completed","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"genomad":"genomad/input_summary","viral_fasta":"genomad/input_summary/input_virus.fna","checkv":"checkv"}}' > 06_viral/stage.json
    """
}

process ROUTE_BINS {
    tag "${meta.sample_id}"
    label 'qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(binning_dir), path(screen_dir), path(ssu_dir), path(viral_dir)

    output:
    tuple val(meta), path('07_routing'), emit: result

    script:
    def sampleJson = JsonOutput.toJson(meta)
    """
    cat > sample.json <<'JSON'
    ${sampleJson}
JSON
    "${params.cli}" route \
        --sample-json sample.json \
        --bins-dir "${binning_dir}/bins" \
        --quickclade-dir "${binning_dir}" \
        --checkeuk-dir "${screen_dir}/checkeuk" \
        --gvclass-dir "${screen_dir}/gvclass" \
        --ssuextract-dir "${ssu_dir}" \
        --viral-dir "${viral_dir}" \
        --output-dir 07_routing
    [[ -s 07_routing/stage.json ]]
    [[ -s 07_routing/evidence.tsv ]]
    for ROUTE in prokaryotic eukaryotic viral unresolved; do
        [[ -d "07_routing/\${ROUTE}" ]]
    done
    """

    stub:
    """
    mkdir -p 07_routing/{prokaryotic,eukaryotic,viral,unresolved}
    printf 'sample_id\tbin_id\troute\n${meta.sample_id}\tbin_1\tunresolved\n' > 07_routing/evidence.tsv
    printf '>stub_contig\nACGTACGTACGT\n' > 07_routing/unresolved/bin_1.fna
    printf '%s\n' '{"stage":"routing","status":"completed","tool":"protist-meta","tool_version":"stub","database_version":null,"command":"stub","outputs":{"evidence":"evidence.tsv","prokaryotic":"prokaryotic","eukaryotic":"eukaryotic","viral":"viral","unresolved":"unresolved"}}' > 07_routing/stage.json
    """
}
