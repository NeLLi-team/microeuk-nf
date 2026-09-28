import groovy.json.JsonOutput

process COLLECT_RECORDS {
    tag "${meta.sample_id}"
    label 'report'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(stage_dirs)

    output:
    tuple val(meta), path("14_collection/${meta.sample_id}.records.json"), emit: result
    tuple val(meta), path('14_collection/sample.json'), emit: sample

    script:
    def sample = [
        sample_id: meta.sample_id,
        run_id: meta.run_id,
        platform: meta.platform,
        reads: meta.reads,
        run_mode: meta.run_mode,
        skip_gene_calling: params.skip_gene_calling,
        skip_annotation: params.skip_annotation,
        workflow_version: meta.workflow_version,
        source_revision: meta.source_revision,
        assembly: meta.assembly,
        rna_reads: meta.rna_reads,
        rna_platform: meta.rna_platform,
        protein_reference: meta.protein_reference,
        protein_lineage: meta.protein_lineage,
        genetic_code: meta.genetic_code,
        softmasked: meta.softmasked,
        basecaller: meta.basecaller,
        basecaller_model: meta.basecaller_model,
        library_prep: meta.library_prep
    ]
    def sampleJson = JsonOutput.toJson(sample)
    def stages = stage_dirs.collect { "\"${it}\"" }.join(' ')
    """
    mkdir -p 14_collection
    cat > 14_collection/sample.json <<'JSON'
    ${sampleJson}
JSON
    "${params.cli}" collect \
        --sample-json 14_collection/sample.json \
        --stages ${stages} \
        --output "14_collection/${meta.sample_id}.records.json"
    [[ -s "14_collection/${meta.sample_id}.records.json" ]]
    """

    stub:
    """
    mkdir -p 14_collection
    printf '%s\n' '{"sample_id":"${meta.sample_id}"}' > 14_collection/sample.json
    printf '%s\n' '{"sample":{"sample_id":"${meta.sample_id}"},"stages":[]}' > "14_collection/${meta.sample_id}.records.json"
    """
}

process BUILD_CATALOG {
    tag 'catalog'
    label 'report'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path records_jsons

    output:
    path 'catalog', emit: result

    script:
    def records = records_jsons.collect { "\"${it}\"" }.join(' ')
    """
    mkdir -p catalog
    "${params.cli}" catalog \
        --records ${records} \
        --output catalog/protist-meta.sqlite
    [[ -s catalog/protist-meta.sqlite ]]
    """

    stub:
    """
    mkdir -p catalog
    printf 'stub\n' > catalog/protist-meta.sqlite
    """
}

process BUILD_REPORT {
    tag 'report'
    label 'report'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path catalog_dir

    output:
    path 'report', emit: result

    script:
    def executionProvenance = params.execution_provenance
        ? "--execution-provenance '" + params.execution_provenance.toString().replace("'", "'\"'\"'") + "'"
        : ''
    """
    "${params.cli}" report \
        --catalog "${catalog_dir}/protist-meta.sqlite" \
        --output-dir report ${executionProvenance}
    [[ -s report/index.html ]]
    """

    stub:
    """
    mkdir -p report
    printf '<!doctype html><title>stub</title>\n' > report/index.html
    """
}
