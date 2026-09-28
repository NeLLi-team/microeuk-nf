process CORE_PREFLIGHT {
    tag "${meta.sample_id}"
    label 'qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    val meta

    output:
    tuple val(meta), path('00_preflight'), emit: result

    script:
    """
    mkdir -p 00_preflight
    [[ -x "${params.cli}" ]]
    [[ -s "${params.quickclade_ref}" ]]
    [[ -s "${params.core_manifest}" ]]
    [[ -s "${params.bbtools_manifest}" ]]
    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        bash "${projectDir}/scripts/require-pixi-tools.sh" chopper seqkit myloasm minimap2 samtools pigz
    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        bash -c 'myloasm --help 2>&1 | awk '\''/--hifi/ {found=1} END {exit !found}'\'''
    pixi run --as-is --manifest-path "${params.bbtools_manifest}" --environment bbtools \
        bash "${projectDir}/scripts/require-pixi-tools.sh" quickbin.sh quickclade.sh bbversion.sh java
    cat > 00_preflight/stage.json <<JSON
    {
      "stage": "preflight",
      "status": "completed",
      "tool": "protist-meta",
      "tool_version": "frozen",
      "database_name": null,
      "database_version": null,
      "command": "validate frozen core and BBTools environments plus local QuickClade reference",
      "outputs": {"report": "stage.json"}
    }
JSON
    """

    stub:
    """
    mkdir -p 00_preflight
    printf '%s\n' '{"stage":"preflight","status":"completed","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"report":"stage.json"}}' > 00_preflight/stage.json
    """
}

process FULL_PREFLIGHT {
    tag "${meta.sample_id}"
    label 'qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    val meta

    output:
    tuple val(meta), path('00_preflight'), emit: result

    script:
    def outputs = params.skip_gene_calling
        ? [report: 'stage.json']
        : [report: 'stage.json', famdb_info: 'famdb-info.txt', famdb_exports: 'famdb-exports.tsv']
    def databaseName = params.skip_gene_calling ? 'null' : '"Dfam"'
    def databaseVersion = params.skip_gene_calling ? 'null' : groovy.json.JsonOutput.toJson(params.dfam_db_version)
    def outputsJson = groovy.json.JsonOutput.toJson(outputs)
    """
    mkdir -p 00_preflight
    [[ -x "${params.cli}" ]]
    [[ -s "${params.quickclade_ref}" ]]
    [[ -s "${params.core_manifest}" ]]
    [[ -s "${params.bbtools_manifest}" ]]
    [[ -s "${params.checkm_manifest}" ]]
    [[ -s "${params.viral_manifest}" ]]
    [[ -s "${params.checkv_manifest}" ]]
    if [[ "${!params.skip_gene_calling}" == true ]]; then
        [[ -s "${params.gene_manifest}" ]]
        [[ -s "${params.braker_manifest}" ]]
    fi
    if [[ "${!params.skip_gene_calling && !params.skip_annotation}" == true ]]; then
        [[ -s "${params.annotation_manifest}" ]]
        [[ -s "${params.interpro_manifest}" ]]
    fi
    [[ -s "${params.checkeuk_app}/pixi.toml" ]]
    [[ -s "${params.gvclass_app}/pixi.toml" ]]
    [[ -s "${params.ssuextract_app}/pixi.toml" ]]
    [[ -s "${params.quickclade_ref}" ]]
    [[ -s "${params.checkm1_db}/hmms/checkm.hmm" ]]
    [[ -s "${params.checkm2_db}" ]]
    [[ -s "${params.gtdbtk_db}/metadata/metadata.txt" ]]
    [[ -d "${params.genomad_db}" ]]
    [[ -d "${params.checkv_db}" ]]
    if [[ "${!params.skip_gene_calling && !params.skip_annotation}" == true ]]; then
        [[ -d "${params.eggnog_db}" ]]
        [[ -d "${params.interproscan_db}" ]]
    fi
    if [[ "${!params.skip_gene_calling}" == true ]]; then
        [[ -d "${params.dfam_db}" ]]
        [[ -s "${params.dfam_db}/dfam40.0.h5" ]]
        [[ -s "${params.dfam_db}/dfam40.curated.consensus.0.h5" ]]
        [[ -s "${params.dfam_db}/PUBLISHER_MD5SUMS" ]]
        [[ -s "${params.dfam_db}/SHA256SUMS" ]]
        [[ -s "${params.dfam_db}/sources.tsv" ]]
    fi

    if [[ "${!params.skip_gene_calling && !params.skip_annotation}" == true ]]; then
        UNRESOLVED=0
        if [[ "${params.interproscan_db_version}" == *unverified* ]]; then
            printf '%s\n' 'InterProScan package/database compatibility is unresolved in conf/databases.yaml.' >&2
            UNRESOLVED=1
        fi
        [[ \${UNRESOLVED} -eq 0 ]] || exit 2
    fi

    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        bash "${projectDir}/scripts/require-pixi-tools.sh" chopper seqkit myloasm minimap2 samtools pigz
    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        bash -c 'myloasm --help 2>&1 | awk '\''/--hifi/ {found=1} END {exit !found}'\'''
    pixi run --as-is --manifest-path "${params.bbtools_manifest}" --environment bbtools \
        bash "${projectDir}/scripts/require-pixi-tools.sh" quickbin.sh quickclade.sh bbversion.sh java
    if [[ "${!params.skip_gene_calling}" == true ]]; then
        pixi run --as-is --manifest-path "${params.braker_manifest}" \
            braker3-longrna-preflight
        pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
            --environment repeatmasker bash "${projectDir}/scripts/require-pixi-tools.sh" \
            BuildDatabase RepeatModeler RepeatClassifier RepeatMasker famdb.py
        bash "${projectDir}/scripts/configure-famdb.sh" check \
            "${params.gene_manifest}" "${params.dfam_db}"
        pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
            --environment repeatmasker famdb.py info > 00_preflight/famdb-info.txt
        pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
            --environment repeatmasker famdb.py -i "${params.dfam_db}" info \
            > 00_preflight/famdb-info-registered.txt
        cmp 00_preflight/famdb-info.txt 00_preflight/famdb-info-registered.txt
        grep -Eq '^FamDB Creation Format Version[[:space:]]*:[[:space:]]*3[.]0[.]0\$' \
            00_preflight/famdb-info.txt
        grep -Eq '^Database[[:space:]]*:[[:space:]]*Dfam\$' 00_preflight/famdb-info.txt
        grep -Eq '^Version[[:space:]]*:[[:space:]]*${params.dfam_db_version}\$' \
            00_preflight/famdb-info.txt
        grep -Eq '^Date[[:space:]]*:[[:space:]]*[^[:space:]].*\$' 00_preflight/famdb-info.txt
        REPEAT_PEPTIDES=\$(
            pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
                --environment repeatmasker famdb.py repeat_peps \
                2> 00_preflight/famdb-repeat-peps.stderr \
            | awk '/^>/ {count++} END {if (count == 0) exit 1; print count}'
        )
        [[ ! -s 00_preflight/famdb-repeat-peps.stderr ]]
        CURATED_CONSENSUS=\$(
            pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
                --environment repeatmasker famdb.py fasta_all \
                2> 00_preflight/famdb-fasta-all.stderr \
            | awk '/^>/ {count++} END {if (count == 0) exit 1; print count}'
        )
        [[ ! -s 00_preflight/famdb-fasta-all.stderr ]]
        printf 'export\trecords\nrepeat_peps\t%s\nfasta_all\t%s\n' \
            "\${REPEAT_PEPTIDES}" "\${CURATED_CONSENSUS}" \
            > 00_preflight/famdb-exports.tsv
        pixi run --as-is --manifest-path "${params.gene_manifest}" --environment prodigal-gv \
            bash "${projectDir}/scripts/require-pixi-tools.sh" prodigal-gv
    fi
    if [[ "${!params.skip_gene_calling && !params.skip_annotation}" == true ]]; then
        pixi run --as-is --manifest-path "${params.annotation_manifest}" --environment eggnog \
            bash "${projectDir}/scripts/require-pixi-tools.sh" emapper.py
        export INTERPROSCAN_DATA_DIR="${params.interproscan_db}"
        pixi run --as-is --manifest-path "${params.interpro_manifest}" \
            interproscan-current-preflight
    fi
    "${params.cli}" route --help >/dev/null

    cat > 00_preflight/stage.json <<JSON
    {
      "stage": "preflight",
      "status": "completed",
      "tool": "protist-meta",
      "tool_version": "frozen full environments",
      "database_name": ${databaseName},
      "database_version": ${databaseVersion},
      "command": "validate enabled frozen environments, database sentinels, and tool dependencies",
      "outputs": ${outputsJson}
    }
JSON
    """

    stub:
    """
    mkdir -p 00_preflight
    printf '%s\n' '{"stage":"preflight","status":"completed","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"report":"stage.json"}}' > 00_preflight/stage.json
    """
}
