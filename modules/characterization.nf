process CHECKEUK_SCREEN {
    tag { meta.sample_id }
    label 'checkeuk'

    input:
    tuple val(meta), path(binning_dir)

    output:
    tuple val(meta), path('checkeuk'), path('checkeuk.tool.tsv'), emit: result

    script:
    def commit = params.checkeuk_app_commit.take(12)
    """
    mkdir -p checkeuk
    BIN_COUNT=\$(find "${binning_dir}/bins" -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l)
    if [[ \${BIN_COUNT} -eq 0 ]]; then
        STATUS=skipped
        REASON=no_bins
        VERSION=not_run
        OUTPUT=checkeuk
    else
        STATUS=completed
        REASON=null
        OUTPUT=checkeuk/result_checkeuk/checkeuk_report.tsv
        pixi run --as-is --manifest-path "${params.checkeuk_app}/pixi.toml" checkeuk \
            -i "\$PWD/${binning_dir}/bins" -o "\$PWD/checkeuk/result" \
            --resources "${params.checkeuk_db}" --threads ${task.cpus} --checkpoint-every 10 \
            --profile-out "\$PWD/checkeuk/profile.tsv"
        [[ -s "\${OUTPUT}" ]]
        [[ \$(awk 'END {print NR - 1}' "\${OUTPUT}") -eq \${BIN_COUNT} ]]
        awk -F '\t' 'NR == 1 {for (i=1; i<=NF; i++) if (\$i=="status") status=i; next}
            status && (\$status == "" || \$status == "error") {failed++}
            END {exit status == 0 || failed > 0}' "\${OUTPUT}"

        VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkeuk_app}/pixi.toml" \
            checkeuk --version | tr '\n\t' '  ' | tr -d '"' | sed 's/[[:space:]]*\$//')
        [[ -n "\${VERSION}" ]]
        VERSION="\${VERSION} (${commit})"
    fi
    cp .command.sh checkeuk/command.sh
    printf '%s\t%s\t%s\t%s\n' "\${VERSION}" "\${STATUS}" "\${REASON}" "\${OUTPUT}" \
        > checkeuk.tool.tsv
    """

    stub:
    """
    mkdir -p checkeuk/result_checkeuk
    printf 'stub\n' > checkeuk/command.sh
    printf 'genome\tstatus\n' > checkeuk/result_checkeuk/checkeuk_report.tsv
    printf 'stub\tcompleted\tnull\tcheckeuk/result_checkeuk/checkeuk_report.tsv\n' > checkeuk.tool.tsv
    """
}

process CHECKM2_SCREEN {
    tag { meta.sample_id }
    label 'checkm2'

    input:
    tuple val(meta), path(binning_dir)

    output:
    tuple val(meta), path('checkm2'), path('checkm2.tool.tsv'), emit: result

    script:
    """
    mkdir -p checkm2
    BIN_COUNT=\$(find "${binning_dir}/bins" -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l)
    if [[ \${BIN_COUNT} -eq 0 ]]; then
        STATUS=skipped
        REASON=no_bins
        VERSION=not_run
        OUTPUT=checkm2
    else
        STATUS=completed
        REASON=null
        OUTPUT=checkm2/output/quality_report.tsv
        CHECKM2_EXIT=0
        pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment checkm2 \
            checkm2 predict --input "${binning_dir}/bins" \
            --output-directory checkm2/output --threads ${task.cpus} \
            --database_path "${params.checkm2_db}" --extension fa --tmpdir /tmp --force \
            > checkm2/command.log 2>&1 || CHECKM2_EXIT=\$?
        printf '%s\n' "\${CHECKM2_EXIT}" > checkm2/exit_code.txt
        cat checkm2/command.log
        if [[ \${CHECKM2_EXIT} -eq 0 ]]; then
            pixi run --as-is --manifest-path "${params.core_manifest}" \
                python -m protist_meta.checkm2 checkm2/output "${binning_dir}/bins"
        elif bash "${projectDir}/scripts/checkm2-no-annotations.sh" "\${CHECKM2_EXIT}" \
            checkm2/output "${binning_dir}/bins"; then
            STATUS=skipped
            REASON=no_diamond_annotations
            OUTPUT=checkm2
        else
            exit "\${CHECKM2_EXIT}"
        fi

        VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
            --environment checkm2 checkm2 --version 2>&1 | tail -n 1 | tr -d '"')
        [[ -n "\${VERSION}" ]]
    fi
    cp .command.sh checkm2/command.sh
    printf '%s\t%s\t%s\t%s\n' "\${VERSION}" "\${STATUS}" "\${REASON}" "\${OUTPUT}" \
        > checkm2.tool.tsv
    """

    stub:
    """
    mkdir -p checkm2/output
    printf 'stub\n' > checkm2/command.sh
    printf 'Name\tCompleteness\tContamination\n' > checkm2/output/quality_report.tsv
    printf 'stub\tcompleted\tnull\tcheckm2/output/quality_report.tsv\n' > checkm2.tool.tsv
    """
}

process GVCLASS_SCREEN {
    tag { meta.sample_id }
    label 'gvclass'

    input:
    tuple val(meta), path(binning_dir)

    output:
    tuple val(meta), path('gvclass'), path('gvclass.tool.tsv'), emit: result

    script:
    def workers = task.cpus as int
    def commit = params.gvclass_app_commit.take(12)
    """
    mkdir -p gvclass/input gvclass/output
    BIN_COUNT=\$(find "${binning_dir}/bins" -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l)
    if [[ \${BIN_COUNT} -eq 0 ]]; then
        STATUS=skipped
        REASON=no_bins
        VERSION=not_run
        OUTPUT=gvclass
    else
        STATUS=completed
        REASON=null
        OUTPUT=gvclass/output/gvclass_summary.tsv
        for FASTA in "${binning_dir}"/bins/*.fa; do
            cp -L "\${FASTA}" "gvclass/input/\$(basename "\${FASTA}" .fa).fna"
        done
        [[ \$(find gvclass/input -type f -name '*.fna' -size +0c | wc -l) -eq \${BIN_COUNT} ]]
        export GVCLASS_RESOURCE_CACHE="\$PWD/tmp/gvclass_cache"
        export GVCLASS_PLAIN_OUTPUT=1
        pixi run --as-is --manifest-path "${params.gvclass_app}/pixi.toml" gvclass \
            "\$PWD/gvclass/input" -o "\$PWD/gvclass/output" \
            --database "${params.gvclass_db}" -t ${task.cpus} -j ${workers} \
            --threads-per-worker 1 --min-length 0 --plain-output
        pixi run --as-is --quiet --manifest-path "${params.core_manifest}" \
            python "${projectDir}/scripts/check-gvclass-reference.py" \
            gvclass/output/run_status.json \
            "${params.gvclass_db}" "${params.gvclass_db_version}"
        [[ -s "\${OUTPUT}" ]]
        [[ \$(awk 'END {print NR - 1}' "\${OUTPUT}") -eq \${BIN_COUNT} ]]
        if [[ -s gvclass/output/gvclass_failed_queries.tsv ]]; then
            [[ \$(awk 'END {print NR - 1}' gvclass/output/gvclass_failed_queries.tsv) -eq 0 ]]
        fi

        VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.gvclass_app}/pixi.toml" \
            gvclass --version | awk '/Software version:/ {sub(/^v/, "", \$NF); print \$NF}')
        [[ -n "\${VERSION}" ]]
        VERSION="\${VERSION} (${commit})"
    fi
    cp .command.sh gvclass/command.sh
    printf '%s\t%s\t%s\t%s\n' "\${VERSION}" "\${STATUS}" "\${REASON}" "\${OUTPUT}" \
        > gvclass.tool.tsv
    """

    stub:
    """
    mkdir -p gvclass/input gvclass/output
    printf 'stub\n' > gvclass/command.sh
    printf 'genome\tdomain\n' > gvclass/output/gvclass_summary.tsv
    printf 'stub\tcompleted\tnull\tgvclass/output/gvclass_summary.tsv\n' > gvclass.tool.tsv
    """
}

process ALL_BIN_SCREEN {
    tag { meta.sample_id }
    label 'report'
    publishDir { "${params.outdir}/${meta.sample_id}" }, mode: 'link'

    input:
    tuple val(meta), path(binning_dir), path(checkeuk_dir), path(checkeuk_report), path(checkm2_dir), path(checkm2_report), path(gvclass_dir), path(gvclass_report)

    output:
    tuple val(meta), path('04_all_bin_screen'), emit: result

    script:
    """
    mkdir -p 04_all_bin_screen/checkeuk 04_all_bin_screen/checkm2 04_all_bin_screen/gvclass
    cp -rL --link "${checkeuk_dir}/." 04_all_bin_screen/checkeuk/
    cp -rL --link "${checkm2_dir}/." 04_all_bin_screen/checkm2/
    cp -rL --link "${gvclass_dir}/." 04_all_bin_screen/gvclass/
    for TOOL in checkeuk checkm2 gvclass; do
        [[ -s "04_all_bin_screen/\${TOOL}/command.sh" ]]
    done

    BIN_COUNT=\$(find "${binning_dir}/bins" -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l)
    IFS=\$'\t' read -r CHECKEUK_VERSION _ _ CHECKEUK_OUTPUT < "${checkeuk_report}"
    IFS=\$'\t' read -r CHECKM2_VERSION _ _ CHECKM2_OUTPUT < "${checkm2_report}"
    IFS=\$'\t' read -r GVCLASS_VERSION _ _ GVCLASS_OUTPUT < "${gvclass_report}"
    if [[ \${BIN_COUNT} -eq 0 ]]; then
        STATUS=skipped
        REASON='"no_bins"'
        TOOL_VERSION=not_run
        COMMAND='No bin-screen tools run because QuickBin produced no accepted bins.'
    else
        STATUS=completed
        REASON=null
        TOOL_VERSION="CheckEUK \${CHECKEUK_VERSION}; CheckM2 \${CHECKM2_VERSION}; GVClass \${GVCLASS_VERSION}"
        COMMAND='all-bin CheckEUK, CheckM2, and GVClass; GVClass physical .fna copies'
    fi
    CHECKEUK_ENTRY=\$(bash "${projectDir}/scripts/tool-entry.sh" checkeuk "${checkeuk_report}" checkeuk_db "${params.checkeuk_db_version}")
    CHECKM2_ENTRY=\$(bash "${projectDir}/scripts/tool-entry.sh" checkm2 "${checkm2_report}" checkm2_db "${params.checkm2_db_version}")
    GVCLASS_ENTRY=\$(bash "${projectDir}/scripts/tool-entry.sh" gvclass "${gvclass_report}" gvclass_db "${params.gvclass_db_version}")
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
        \${CHECKEUK_ENTRY},
        \${CHECKM2_ENTRY},
        \${GVCLASS_ENTRY}
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
    for TOOL in checkeuk checkm2 gvclass; do
        printf 'stub\n' > "04_all_bin_screen/\${TOOL}/command.sh"
    done
    printf '%s\n' '{"stage":"all_bin_screen","status":"completed","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","reason":null,"tools":{"checkeuk":{"version":"stub","status":"completed","reason":null,"database_name":"checkeuk_db","database_version":"stub","output":"checkeuk/result_checkeuk/checkeuk_report.tsv","command_file":"checkeuk/command.sh"},"checkm2":{"version":"stub","status":"completed","reason":null,"database_name":"checkm2_db","database_version":"stub","output":"checkm2/output/quality_report.tsv","command_file":"checkm2/command.sh"},"gvclass":{"version":"stub","status":"completed","reason":null,"database_name":"gvclass_db","database_version":"stub","output":"gvclass/output/gvclass_summary.tsv","command_file":"gvclass/command.sh"}},"outputs":{"checkeuk":"checkeuk/result_checkeuk/checkeuk_report.tsv","checkm2":"checkm2/output/quality_report.tsv","gvclass":"gvclass/output/gvclass_summary.tsv"}}' > 04_all_bin_screen/stage.json
    """
}

process SSU_EXTRACT {
    tag { meta.sample_id }
    label 'ssu'
    publishDir { "${params.outdir}/${meta.sample_id}" }, mode: 'link'

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
    def commit = params.ssuextract_app_commit.take(12)
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
      "tool_version": "${params.ssuextract_app_version} (${commit})",
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
    tag { meta.sample_id }
    label 'viral'
    publishDir { "${params.outdir}/${meta.sample_id}" }, mode: 'link'

    input:
    tuple val(meta), path(assembly_dir)

    output:
    tuple val(meta), path('06_viral'), emit: result

    script:
    def checkvManifest = params.checkv_manifest ?: params.viral_manifest
    """
    mkdir -p 06_viral/genomad 06_viral/checkv
    cp "${assembly_dir}/assembly_primary.fa" input.fna
    pixi run --as-is --manifest-path "${params.viral_manifest}" --environment genomad \
        genomad end-to-end --cleanup --splits 2 --threads ${task.cpus} \
        input.fna 06_viral/genomad "${params.genomad_db}"
    [[ -s 06_viral/genomad/input_summary/input_virus_summary.tsv ]]

    VIRUS_FASTA=06_viral/genomad/input_summary/input_virus.fna
    CHECKV_STATUS=skipped
    if [[ -s \${VIRUS_FASTA} ]]; then
        pixi run --as-is --manifest-path "${checkvManifest}" --environment checkv \
            checkv end_to_end "\${VIRUS_FASTA}" 06_viral/checkv \
            -t ${task.cpus} -d "${params.checkv_db}"
        [[ -s 06_viral/checkv/quality_summary.tsv ]]
        CHECKV_STATUS=completed
    fi
    GENOMAD_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.viral_manifest}" \
        --environment genomad genomad --version | tr '\n' ' ' | sed 's/"/\\\\"/g')
    CHECKV_VERSION=\$(pixi run --as-is --quiet --manifest-path "${checkvManifest}" \
        --environment checkv \
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

process CONTIG_CLASSIFY {
    tag { meta.sample_id }
    label 'contig'
    publishDir { "${params.outdir}/${meta.sample_id}" }, mode: 'link'

    input:
    tuple val(meta), path(assembly_dir)

    output:
    tuple val(meta), path('15_contig_classification'), emit: result

    script:
    """
    mkdir -p 15_contig_classification
    pixi run --as-is --manifest-path "${params.contig_manifest}" --environment tiara \
        tiara -i "${assembly_dir}/assembly_primary.fa" \
        -o 15_contig_classification/tiara.tsv -m ${params.tiara_min_length} -t ${task.cpus}
    [[ \$(head -n 1 15_contig_classification/tiara.tsv) == \$'sequence_id\tclass_fst_stage\tclass_snd_stage' ]]
    [[ -f 15_contig_classification/log_tiara.tsv ]]
    TIARA_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.contig_manifest}" \
        --environment tiara \
        python -c 'from importlib.metadata import version; print(version("tiara"))')
    [[ -n "\${TIARA_VERSION}" ]]
    cat > 15_contig_classification/stage.json <<JSON
    {
      "stage": "tiara",
      "status": "completed",
      "tool": "Tiara",
      "tool_version": "\${TIARA_VERSION}",
      "database_name": null,
      "database_version": null,
      "command": "tiara -i assembly_primary.fa -o tiara.tsv -m ${params.tiara_min_length} -t ${task.cpus}",
      "outputs": {
        "classification": "tiara.tsv",
        "log": "log_tiara.tsv"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 15_contig_classification
    printf 'sequence_id\tclass_fst_stage\tclass_snd_stage\n' > 15_contig_classification/tiara.tsv
    : > 15_contig_classification/log_tiara.tsv
    printf '%s\n' '{"stage":"tiara","status":"completed","tool":"Tiara","tool_version":"stub","database_name":null,"database_version":null,"command":"stub","outputs":{"classification":"tiara.tsv","log":"log_tiara.tsv"}}' > 15_contig_classification/stage.json
    """
}

process ROUTE_BINS {
    tag { meta.sample_id }
    label 'qc'
    publishDir { "${params.outdir}/${meta.sample_id}" }, mode: 'link'

    input:
    tuple val(meta), path(binning_dir), path(screen_dir), path(ssu_dir), path(viral_dir), path(tiara_dir)

    output:
    tuple val(meta), path('07_routing'), emit: result

    script:
    def sampleJson = groovy.json.JsonOutput.toJson(meta)
    """
    cat > sample.json <<'JSON'
    ${sampleJson}
JSON
    "${params.cli}" route \
        --sample-json sample.json \
        --bins-dir "${binning_dir}/bins" \
        --quickclade-dir "${binning_dir}" \
        --checkeuk-dir "${screen_dir}/checkeuk" \
        --checkm2-dir "${screen_dir}/checkm2" \
        --gvclass-dir "${screen_dir}/gvclass" \
        --ssuextract-dir "${ssu_dir}" \
        --viral-dir "${viral_dir}" \
        --tiara-dir "${tiara_dir}" \
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

process CHECKM1_QC {
    tag { meta.sample_id }
    label 'checkm1'

    input:
    tuple val(meta), path(routing_dir)

    output:
    tuple val(meta), path('checkm1'), path('checkm1.tool.tsv'), emit: result

    script:
    def pplacerThreads = Math.min(8, task.cpus as int)
    """
    mkdir -p checkm1
    BIN_COUNT=\$(find "${routing_dir}/prokaryotic" -maxdepth 1 -type f -name '*.fna' -size +0c | wc -l)
    STATUS=skipped
    REASON=no_routed_prokaryotic_bins
    if [[ \${BIN_COUNT} -gt 0 ]]; then
        STATUS=completed
        REASON=null
        export CHECKM_DATA_PATH="${params.checkm1_db}"
        TMPDIR=/tmp pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment checkm1 \
            checkm lineage_wf --reduced_tree --extension fna --threads ${task.cpus} \
            --pplacer_threads ${pplacerThreads} --tab_table --file "\$PWD/checkm1/quality_report.tsv" \
            --tmpdir /tmp "${routing_dir}/prokaryotic" \
            "\$PWD/checkm1/output"
        [[ -s checkm1/quality_report.tsv ]]
    fi

    VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
        --environment checkm1 python -c \
        'from pathlib import Path; import checkm; print((Path(checkm.__file__).parent / "VERSION").read_text().splitlines()[0])')
    [[ -n "\${VERSION}" ]]
    cp .command.sh checkm1/command.sh
    printf '%s\t%s\t%s\t%s\n' "\${VERSION}" "\${STATUS}" "\${REASON}" checkm1/quality_report.tsv \
        > checkm1.tool.tsv
    """

    stub:
    """
    mkdir -p checkm1
    printf 'stub\n' > checkm1/command.sh
    printf 'Bin Id\tMarker lineage\n' > checkm1/quality_report.tsv
    printf 'stub\tskipped\tno_routed_prokaryotic_bins\tcheckm1/quality_report.tsv\n' > checkm1.tool.tsv
    """
}

process GTDBTK_CLASSIFY {
    tag { meta.sample_id }
    label 'gtdbtk'

    input:
    tuple val(meta), path(routing_dir)

    output:
    tuple val(meta), path('gtdbtk'), path('gtdbtk.tool.tsv'), emit: result

    script:
    def pplacerCpus = Math.min(8, task.cpus as int)
    """
    mkdir -p gtdbtk
    BIN_COUNT=\$(find "${routing_dir}/prokaryotic" -maxdepth 1 -type f -name '*.fna' -size +0c | wc -l)
    STATUS=skipped
    REASON=no_routed_prokaryotic_bins
    if [[ \${BIN_COUNT} -gt 0 ]]; then
        STATUS=completed
        REASON=null
        export GTDBTK_DATA_PATH="${params.gtdbtk_db}"
        TMPDIR=/tmp pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment gtdbtk \
            gtdbtk classify_wf --genome_dir "${routing_dir}/prokaryotic" \
            --out_dir gtdbtk --extension fna \
            --cpus ${task.cpus} --pplacer_cpus ${pplacerCpus} --tmpdir /tmp
        [[ \$(find gtdbtk -type f -name 'gtdbtk.*.summary.tsv' -size +0c | wc -l) -gt 0 ]]
    fi

    VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
        --environment gtdbtk python -c \
        'from importlib.metadata import version; print(version("gtdbtk"))')
    [[ -n "\${VERSION}" ]]
    cp .command.sh gtdbtk/command.sh
    printf '%s\t%s\t%s\t%s\n' "\${VERSION}" "\${STATUS}" "\${REASON}" gtdbtk \
        > gtdbtk.tool.tsv
    """

    stub:
    """
    mkdir -p gtdbtk
    printf 'stub\n' > gtdbtk/command.sh
    printf 'stub\tskipped\tno_routed_prokaryotic_bins\tgtdbtk\n' > gtdbtk.tool.tsv
    """
}

process SYMCLATRON_CLASSIFY {
    tag { meta.sample_id }
    label 'symclatron'

    input:
    tuple val(meta), path(routing_dir)

    output:
    tuple val(meta), path('symclatron'), path('symclatron.tool.tsv'), emit: result

    script:
    """
    mkdir -p symclatron
    BIN_COUNT=\$(find "${routing_dir}/prokaryotic" -maxdepth 1 -type f -name '*.fna' -size +0c | wc -l)
    STATUS=skipped
    REASON=no_routed_prokaryotic_bins
    if [[ \${BIN_COUNT} -gt 0 ]]; then
        STATUS=completed
        REASON=null
        pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment symclatron \
            symclatron classify --genome-dir "${routing_dir}/prokaryotic" \
            --input-kind contigs --input-ext .fna --threads ${task.cpus} \
            --confidence-threshold 0.725 --output-dir symclatron
        [[ -s symclatron/symclatron_results.tsv ]]
        [[ \$(awk 'END {print NR - 1}' symclatron/symclatron_results.tsv) -eq \${BIN_COUNT} ]]
    fi

    VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
        --environment symclatron symclatron --version | tr '\n\t' '  ' | tr -d '"' \
        | sed 's/[[:space:]]*\$//')
    [[ -n "\${VERSION}" ]]
    cp .command.sh symclatron/command.sh
    printf '%s\t%s\t%s\t%s\n' "\${VERSION}" "\${STATUS}" "\${REASON}" \
        symclatron/symclatron_results.tsv > symclatron.tool.tsv
    """

    stub:
    """
    mkdir -p symclatron
    printf 'stub\n' > symclatron/command.sh
    printf 'taxon_oid\tcompleteness_UNI56\tclassification\tconfidence\tpasses_confidence_threshold\tclassification_thresholded\n' > symclatron/symclatron_results.tsv
    printf 'stub\tskipped\tno_routed_prokaryotic_bins\tsymclatron/symclatron_results.tsv\n' > symclatron.tool.tsv
    """
}

process PROKARYOTE_CHARACTERIZATION {
    tag { meta.sample_id }
    label 'report'
    publishDir { "${params.outdir}/${meta.sample_id}" }, mode: 'link'

    input:
    tuple val(meta), path(routing_dir), path(checkm1_dir), path(checkm1_report), path(gtdbtk_dir), path(gtdbtk_report), path(symclatron_dir), path(symclatron_report)

    output:
    tuple val(meta), path('08_prokaryote_qc'), emit: result

    script:
    """
    mkdir -p 08_prokaryote_qc/checkm1 08_prokaryote_qc/gtdbtk 08_prokaryote_qc/symclatron
    cp -rL --link "${checkm1_dir}/." 08_prokaryote_qc/checkm1/
    cp -rL --link "${gtdbtk_dir}/." 08_prokaryote_qc/gtdbtk/
    cp -rL --link "${symclatron_dir}/." 08_prokaryote_qc/symclatron/
    for TOOL in checkm1 gtdbtk symclatron; do
        [[ -s "08_prokaryote_qc/\${TOOL}/command.sh" ]]
    done

    BIN_COUNT=\$(find "${routing_dir}/prokaryotic" -maxdepth 1 -type f -name '*.fna' -size +0c | wc -l)
    IFS=\$'\t' read -r CHECKM1_VERSION _ _ CHECKM1_OUTPUT < "${checkm1_report}"
    IFS=\$'\t' read -r GTDBTK_VERSION _ _ GTDBTK_OUTPUT < "${gtdbtk_report}"
    IFS=\$'\t' read -r SYMCLATRON_VERSION _ _ SYMCLATRON_OUTPUT < "${symclatron_report}"
    STATUS=completed
    REASON=null
    if [[ \${BIN_COUNT} -eq 0 ]]; then
        STATUS=skipped
        REASON='"no_routed_prokaryotic_bins"'
    fi
    CHECKM1_ENTRY=\$(bash "${projectDir}/scripts/tool-entry.sh" checkm1 "${checkm1_report}" checkm1_db "${params.checkm1_db_version}")
    GTDBTK_ENTRY=\$(bash "${projectDir}/scripts/tool-entry.sh" gtdbtk "${gtdbtk_report}" gtdbtk_db "${params.gtdbtk_db_version}")
    SYMCLATRON_ENTRY=\$(bash "${projectDir}/scripts/tool-entry.sh" symclatron "${symclatron_report}" symclatron_bundled bundled)
    cat > 08_prokaryote_qc/stage.json <<JSON
    {
      "stage": "prokaryote_qc",
      "status": "\${STATUS}",
      "tool": "CheckM1;GTDB-Tk;Symclatron",
      "tool_version": "CheckM1 \${CHECKM1_VERSION}; GTDB-Tk \${GTDBTK_VERSION}; Symclatron \${SYMCLATRON_VERSION}",
      "database_name": "checkm1_db;gtdbtk_db;symclatron_bundled",
      "database_version": "CheckM1 ${params.checkm1_db_version}; GTDB-Tk ${params.gtdbtk_db_version}; Symclatron bundled",
      "command": "bounded per-sample batch over routed prokaryotic bins",
      "reason": \${REASON},
      "tools": {
        \${CHECKM1_ENTRY},
        \${GTDBTK_ENTRY},
        \${SYMCLATRON_ENTRY}
      },
      "outputs": {
        "checkm1": "\${CHECKM1_OUTPUT}",
        "gtdbtk": "\${GTDBTK_OUTPUT}",
        "symclatron": "\${SYMCLATRON_OUTPUT}"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 08_prokaryote_qc/{checkm1,gtdbtk,symclatron}
    for TOOL in checkm1 gtdbtk symclatron; do
        printf 'stub\n' > "08_prokaryote_qc/\${TOOL}/command.sh"
    done
    printf 'Bin Id\tMarker lineage\n' > 08_prokaryote_qc/checkm1/quality_report.tsv
    printf 'taxon_oid\tcompleteness_UNI56\tclassification\tconfidence\tpasses_confidence_threshold\tclassification_thresholded\n' > 08_prokaryote_qc/symclatron/symclatron_results.tsv
    printf '%s\n' '{"stage":"prokaryote_qc","status":"skipped","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","reason":"no_routed_prokaryotic_bins","tools":{"checkm1":{"version":"stub","status":"skipped","reason":"no_routed_prokaryotic_bins","database_name":"checkm1_db","database_version":"stub","output":"checkm1/quality_report.tsv","command_file":"checkm1/command.sh"},"gtdbtk":{"version":"stub","status":"skipped","reason":"no_routed_prokaryotic_bins","database_name":"gtdbtk_db","database_version":"stub","output":"gtdbtk","command_file":"gtdbtk/command.sh"},"symclatron":{"version":"stub","status":"skipped","reason":"no_routed_prokaryotic_bins","database_name":"symclatron_bundled","database_version":"bundled","output":"symclatron/symclatron_results.tsv","command_file":"symclatron/command.sh"}},"outputs":{"checkm1":"checkm1/quality_report.tsv","gtdbtk":"gtdbtk","symclatron":"symclatron/symclatron_results.tsv"}}' > 08_prokaryote_qc/stage.json
    """
}

process EUKARYOTE_QC {
    tag { meta.sample_id }
    label 'eukqc'
    publishDir { "${params.outdir}/${meta.sample_id}" }, mode: 'link'

    input:
    tuple val(meta), path(routing_dir)

    output:
    tuple val(meta), path('16_eukaryote_qc'), emit: result

    script:
    """
    mkdir -p 16_eukaryote_qc/tidk
    shopt -s nullglob
    BINS=("${routing_dir}"/eukaryotic/*.fna)
    if [[ \${#BINS[@]} -eq 0 ]]; then
        STATUS=skipped
        REASON='"no_routed_eukaryotic_bins"'
        mkdir -p 16_eukaryote_qc/busco
    else
        STATUS=completed
        REASON=null
        pixi run --as-is --manifest-path "${params.eukqc_manifest}" --environment busco \
            busco -i "\$PWD/${routing_dir}/eukaryotic" -o busco --out_path "\$PWD/16_eukaryote_qc" \
            -m genome -l "${params.busco_db}" --offline --opt-out-run-stats -c ${task.cpus}
        [[ -s 16_eukaryote_qc/busco/batch_summary.txt ]]
        [[ \$(awk 'END {print NR - 1}' 16_eukaryote_qc/busco/batch_summary.txt) -eq \${#BINS[@]} ]]

        printf 'bin_id\ttop_repeat\ttop_repeat_count\n' > 16_eukaryote_qc/tidk/tidk_summary.tsv
        for FASTA in "\${BINS[@]}"; do
            BIN=\$(basename "\${FASTA}" .fna)
            pixi run --as-is --manifest-path "${params.eukqc_manifest}" --environment tidk \
                tidk explore --minimum 5 --maximum 12 --threshold 5 "\${FASTA}" \
                > "16_eukaryote_qc/tidk/\${BIN}.explore.tsv"
            awk -F '\t' -v bin="\${BIN}" '
                NR == 1 {next}
                \$2 + 0 > best || (best > 0 && \$2 + 0 == best && \$1 < repeat) {best = \$2 + 0; repeat = \$1}
                END {if (best > 0) printf "%s\\t%s\\t%s\\n", bin, repeat, best; else printf "%s\\t\\t\\n", bin}
            ' "16_eukaryote_qc/tidk/\${BIN}.explore.tsv" >> 16_eukaryote_qc/tidk/tidk_summary.tsv
        done
    fi

    BUSCO_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.eukqc_manifest}" \
        --environment busco busco --version | tr '\n\t' '  ' | tr -d '"' | sed 's/[[:space:]]*\$//')
    TIDK_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.eukqc_manifest}" \
        --environment tidk tidk --version 2>/dev/null | tail -n 1 | tr -d '"')
    [[ -n "\${BUSCO_VERSION}" && -n "\${TIDK_VERSION}" ]]
    cat > 16_eukaryote_qc/stage.json <<JSON
    {
      "stage": "eukaryote_qc",
      "status": "\${STATUS}",
      "tool": "BUSCO;tidk",
      "tool_version": "\${BUSCO_VERSION}; \${TIDK_VERSION}",
      "database_name": "busco_db",
      "database_version": "BUSCO ${params.busco_db_version}",
      "command": "BUSCO batch mode over routed eukaryotic bins (genome mode, offline); tidk explore per bin",
      "reason": \${REASON},
      "tools": {
        "busco": {"version": "\${BUSCO_VERSION}", "status": "\${STATUS}", "reason": \${REASON}, "database_name": "busco_db", "database_version": "${params.busco_db_version}", "output": "busco/batch_summary.txt"},
        "tidk": {"version": "\${TIDK_VERSION}", "status": "\${STATUS}", "reason": \${REASON}, "database_name": null, "database_version": null, "output": "tidk/tidk_summary.tsv"}
      },
      "outputs": {
        "busco": "busco/batch_summary.txt",
        "tidk": "tidk/tidk_summary.tsv"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 16_eukaryote_qc/busco 16_eukaryote_qc/tidk
    printf '%s\n' '{"stage":"eukaryote_qc","status":"skipped","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","reason":"no_routed_eukaryotic_bins","tools":{"busco":{"version":"stub","status":"skipped","reason":"no_routed_eukaryotic_bins","database_name":"busco_db","database_version":"stub","output":"busco/batch_summary.txt"},"tidk":{"version":"stub","status":"skipped","reason":"no_routed_eukaryotic_bins","database_name":null,"database_version":null,"output":"tidk/tidk_summary.tsv"}},"outputs":{"busco":"busco/batch_summary.txt","tidk":"tidk/tidk_summary.tsv"}}' > 16_eukaryote_qc/stage.json
    """
}
