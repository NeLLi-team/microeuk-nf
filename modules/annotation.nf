process FUNCTIONAL_ANNOTATION {
    tag "${meta.sample_id}"
    label 'annotation'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(prodigal_dir), path(braker_dir)

    output:
    tuple val(meta), path('13_functional_annotation'), emit: result

    script:
    """
    export INTERPROSCAN_DATA_DIR="${params.interproscan_db}"
    PRODIGAL_COUNT=\$(awk -F '\t' \
        'NR > 1 && \$6 == "included" {count++} END {print count + 0}' \
        "${prodigal_dir}/inputs/sources.tsv")
    BRAKER_COUNT=\$(find "${braker_dir}/proteins" -maxdepth 1 \
        -type f -name '*.faa' -size +0c | wc -l)
    STATUS=skipped
    REASON=no_called_proteins
    if [[ \$((PRODIGAL_COUNT + BRAKER_COUNT)) -gt 0 ]]; then
        "${params.cli}" merge-proteins \
            --prodigal-dir "${prodigal_dir}" \
            --braker-dir "${braker_dir}" \
            --output-dir 13_functional_annotation
        [[ -e 13_functional_annotation/proteins.faa ]]
        [[ -s 13_functional_annotation/protein-map.tsv ]]
        mkdir -p 13_functional_annotation/{eggnog,interpro}

        if [[ -s 13_functional_annotation/proteins.faa ]]; then
            STATUS=completed
            REASON=
            pixi run --as-is --manifest-path "${params.annotation_manifest}" --environment eggnog \
                emapper.py -i 13_functional_annotation/proteins.faa \
                --output annotations --output_dir 13_functional_annotation/eggnog \
                --data_dir "${params.eggnog_db}" --cpu ${task.cpus} --override
            [[ -s 13_functional_annotation/eggnog/annotations.emapper.annotations ]]

            pixi run --as-is --manifest-path "${params.interpro_manifest}" \
                interproscan.sh -i 13_functional_annotation/proteins.faa \
                -d 13_functional_annotation/interpro -f TSV,GFF3 -dp -cpu ${task.cpus}
            [[ \$(find 13_functional_annotation/interpro -type f -name '*.tsv' | wc -l) -gt 0 ]]
        fi
    else
        mkdir -p 13_functional_annotation/{eggnog,interpro}
        : > 13_functional_annotation/proteins.faa
        printf 'normalized_id\tnative_id\tnative_gene_id\tcaller\tsource_id\tbin_id\n' \
            > 13_functional_annotation/protein-map.tsv
    fi
    EGGNOG_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.annotation_manifest}" \
        --environment eggnog emapper.py --version --data_dir "${params.eggnog_db}" \
        2>&1 | tr '\n' ' ' | sed 's/"/\\\\"/g')
    INTERPRO_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.interpro_manifest}" \
        interproscan.sh --version 2>&1 | tr '\n' ' ' | sed 's/"/\\\\"/g')

    cat > 13_functional_annotation/stage.json <<JSON
    {
      "stage": "functional_annotation",
      "status": "\${STATUS}",
      "tool": "eggNOG-mapper;InterProScan",
      "tool_version": "eggNOG-mapper \${EGGNOG_VERSION}; InterProScan \${INTERPRO_VERSION}",
      "database_name": "eggnog_db;interproscan_db",
      "database_version": "eggNOG ${params.eggnog_db_version}; InterPro ${params.interproscan_db_version}",
      "command": "namespace completed Prodigal-GV and BRAKER3 proteins; remove one terminal stop for annotation; reject remaining stop markers; batch eggNOG-mapper and InterProScan",
      "tools": {
        "eggnog_mapper": {"version": "\${EGGNOG_VERSION}", "database_name": "eggnog_db", "database_version": "${params.eggnog_db_version}", "output": "eggnog/annotations.emapper.annotations"},
        "interproscan": {"version": "\${INTERPRO_VERSION}", "database_name": "interproscan_db", "database_version": "${params.interproscan_db_version}", "output": "interpro"}
      },
      "outputs": {
        "proteins": "proteins.faa",
        "protein_map": "protein-map.tsv",
        "eggnog": "eggnog/annotations.emapper.annotations",
        "interpro": "interpro"
      },
      "reason": "\${REASON}"
    }
JSON
    """

    stub:
    """
    mkdir -p 13_functional_annotation/{eggnog,interpro}
    : > 13_functional_annotation/proteins.faa
    printf 'normalized_id\tnative_id\tnative_gene_id\tcaller\tsource_id\tbin_id\n' \
        > 13_functional_annotation/protein-map.tsv
    printf '%s\n' '{"stage":"functional_annotation","status":"skipped","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"proteins":"proteins.faa","protein_map":"protein-map.tsv","eggnog":"eggnog","interpro":"interpro"},"reason":"no_called_proteins"}' > 13_functional_annotation/stage.json
    """
}
