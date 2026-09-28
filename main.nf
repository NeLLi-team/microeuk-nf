nextflow.enable.dsl = 2

include { CORE_PREFLIGHT; FULL_PREFLIGHT } from './modules/preflight'
include { READ_QC; ASSEMBLE; MAP_DEPTH_QUICKBIN } from './modules/core'
include {
    ALL_BIN_SCREEN;
    SSU_EXTRACT;
    VIRAL_SCREEN;
    ROUTE_BINS
} from './modules/characterization'
include {
    PROKARYOTE_CHARACTERIZATION;
    REPEAT_MASK_EUKARYOTES;
    MAP_RNA_EVIDENCE;
    BRAKER3_EUKARYOTES;
    PRODIGAL_GV_GENES
} from './modules/gene_calling'
include { FUNCTIONAL_ANNOTATION } from './modules/annotation'
include { COLLECT_RECORDS; BUILD_CATALOG; BUILD_REPORT } from './modules/catalog'

workflow {
    if (!params.input) {
        error 'Required parameter --input was not provided.'
    }
    if (!(params.run_mode in ['core', 'full'])) {
        error "--run_mode must be core or full, received: ${params.run_mode}"
    }
    for (name in ['skip_gene_calling', 'skip_annotation']) {
        if (!(params[name] instanceof Boolean)) {
            error "--${name} must be a boolean, received: ${params[name]}"
        }
        def prepared = params["prepared_${name}"]
        if (prepared != null && (!(prepared instanceof Boolean) || prepared != params[name])) {
            error "--${name} differs from the prepared run; prepare a new run with the requested options."
        }
    }

    runId = (params.run_id ?: workflow.runName).toString().replaceAll(/[^A-Za-z0-9_.-]/, '_')
    samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def platform = row.platform.toString().toLowerCase()
            if (!(platform in ['ont', 'pacbio_hifi'])) {
                error "Unsupported platform for ${row.sample_id}: ${row.platform}"
            }
            def meta = [
                sample_id: row.sample_id.toString(),
                run_id: runId,
                platform: platform,
                reads: row.reads.toString(),
                assembly: row.assembly?.toString() ?: '',
                rna_reads: row.rna_reads?.toString() ?: '',
                rna_platform: row.rna_platform?.toString()?.toLowerCase() ?: '',
                protein_reference: row.protein_reference?.toString() ?: '',
                protein_lineage: row.protein_lineage?.toString() ?: '',
                genetic_code: row.genetic_code?.toString() ?: 'auto',
                softmasked: row.softmasked?.toString() ?: '',
                basecaller: row.basecaller?.toString() ?: '',
                basecaller_model: row.basecaller_model?.toString() ?: '',
                library_prep: row.library_prep?.toString() ?: '',
                run_mode: params.run_mode,
                workflow_version: params.workflow_version,
                source_revision: params.source_revision
            ]
            tuple(
                meta,
                file(meta.reads, checkIfExists: true),
                meta.assembly ? file(meta.assembly, checkIfExists: true) : [],
                meta.rna_reads ? file(meta.rna_reads, checkIfExists: true) : [],
                meta.protein_reference ? file(meta.protein_reference, checkIfExists: true) : [],
                meta.softmasked ? file(meta.softmasked, checkIfExists: true) : []
            )
        }

    sampleMeta = samples.map { meta, reads, assembly, rna, proteins, softmasked -> meta }
    if (params.run_mode == 'full') {
        FULL_PREFLIGHT(sampleMeta)
        preflight = FULL_PREFLIGHT.out.result
    } else {
        CORE_PREFLIGHT(sampleMeta)
        preflight = CORE_PREFLIGHT.out.result
    }

    readsBySample = samples.map { meta, reads, assembly, rna, proteins, softmasked ->
        tuple(meta.sample_id, meta, reads)
    }
    preflightBySample = preflight.map { meta, stageDir -> tuple(meta.sample_id, stageDir) }
    qcInput = readsBySample
        .join(preflightBySample)
        .map { sampleId, meta, reads, stageDir -> tuple(meta, reads, stageDir) }
    READ_QC(qcInput)
    suppliedAssemblies = samples.map { meta, reads, assembly, rna, proteins, softmasked ->
        tuple(meta.sample_id, assembly)
    }
    assemblyInput = READ_QC.out.result
        .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
        .join(suppliedAssemblies)
        .map { sampleId, meta, readQcDir, suppliedAssembly ->
            tuple(meta, readQcDir, suppliedAssembly)
        }
    ASSEMBLE(assemblyInput)
    mappingInput = ASSEMBLE.out.result
        .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
        .join(READ_QC.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
        .map { sampleId, meta, assemblyDir, readQcDir -> tuple(meta, assemblyDir, readQcDir) }
    MAP_DEPTH_QUICKBIN(mappingInput)

    if (params.run_mode == 'full') {
        ALL_BIN_SCREEN(MAP_DEPTH_QUICKBIN.out.result)
        SSU_EXTRACT(ASSEMBLE.out.result)
        VIRAL_SCREEN(ASSEMBLE.out.result)

        routeInput = MAP_DEPTH_QUICKBIN.out.result
            .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
            .join(ALL_BIN_SCREEN.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            .join(SSU_EXTRACT.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            .join(VIRAL_SCREEN.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            .map { sampleId, meta, binningDir, screenDir, ssuDir, viralDir ->
                tuple(meta, binningDir, screenDir, ssuDir, viralDir)
            }
        ROUTE_BINS(routeInput)
        PROKARYOTE_CHARACTERIZATION(ROUTE_BINS.out.result)

        geneStages = Channel.empty()
        annotationStages = Channel.empty()
        if (!params.skip_gene_calling) {
            geneRouting = ROUTE_BINS.out.result
                .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                .join(PROKARYOTE_CHARACTERIZATION.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
                .map { sampleId, meta, routingDir, prokDir -> tuple(sampleId, meta, routingDir) }
            suppliedSoftmasked = samples.map { meta, reads, assembly, rna, proteins, softmasked ->
                tuple(meta.sample_id, softmasked)
            }
            maskingInput = geneRouting
                .join(suppliedSoftmasked)
                .map { sampleId, meta, routingDir, softmasked -> tuple(meta, routingDir, softmasked) }
            REPEAT_MASK_EUKARYOTES(maskingInput)

            rnaInputBase = ASSEMBLE.out.result
                .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                .join(REPEAT_MASK_EUKARYOTES.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            suppliedRna = samples.map { meta, reads, assembly, rna, proteins, softmasked ->
                tuple(meta.sample_id, rna)
            }
            rnaInput = rnaInputBase
                .join(suppliedRna)
                .map { sampleId, meta, assemblyDir, maskingDir, rna -> tuple(meta, assemblyDir, maskingDir, rna) }
            MAP_RNA_EVIDENCE(rnaInput)

            brakerInputBase = REPEAT_MASK_EUKARYOTES.out.result
                .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                .join(MAP_RNA_EVIDENCE.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            suppliedProteins = samples.map { meta, reads, assembly, rna, proteins, softmasked ->
                tuple(meta.sample_id, proteins)
            }
            brakerInput = brakerInputBase
                .join(suppliedProteins)
                .map { sampleId, meta, maskingDir, rnaDir, proteins -> tuple(meta, maskingDir, rnaDir, proteins) }
            BRAKER3_EUKARYOTES(brakerInput)

            prodigalInput = geneRouting
                .join(VIRAL_SCREEN.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
                .map { sampleId, meta, routingDir, viralDir -> tuple(meta, routingDir, viralDir) }
            PRODIGAL_GV_GENES(prodigalInput)

            geneStages = REPEAT_MASK_EUKARYOTES.out.result.mix(
                MAP_RNA_EVIDENCE.out.result,
                BRAKER3_EUKARYOTES.out.result,
                PRODIGAL_GV_GENES.out.result
            )
            if (!params.skip_annotation) {
                annotationInput = PRODIGAL_GV_GENES.out.result
                    .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                    .join(BRAKER3_EUKARYOTES.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
                    .map { sampleId, meta, prodigalDir, brakerDir -> tuple(meta, prodigalDir, brakerDir) }
                FUNCTIONAL_ANNOTATION(annotationInput)
                annotationStages = FUNCTIONAL_ANNOTATION.out.result
            }
        }

        stageDirs = preflight.mix(
            READ_QC.out.result,
            ASSEMBLE.out.result,
            MAP_DEPTH_QUICKBIN.out.result,
            ALL_BIN_SCREEN.out.result,
            SSU_EXTRACT.out.result,
            VIRAL_SCREEN.out.result,
            ROUTE_BINS.out.result,
            PROKARYOTE_CHARACTERIZATION.out.result,
            geneStages,
            annotationStages
        )
    } else {
        stageDirs = preflight.mix(
            READ_QC.out.result,
            ASSEMBLE.out.result,
            MAP_DEPTH_QUICKBIN.out.result
        )
    }

    groupedStages = stageDirs
        .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
        .groupTuple(by: 0)
        .map { sampleId, metas, dirs -> tuple(metas[0], dirs) }
    COLLECT_RECORDS(groupedStages)

    allRecords = COLLECT_RECORDS.out.result
        .map { meta, records -> records }
        .collect()
    BUILD_CATALOG(allRecords)
    BUILD_REPORT(BUILD_CATALOG.out.result)
}
