include { CORE_PREFLIGHT; FULL_PREFLIGHT } from './modules/preflight'
include { READ_QC; ASSEMBLE; MAP_DEPTH_QUICKBIN } from './modules/core'
include {
    CHECKEUK_SCREEN;
    CHECKM2_SCREEN;
    GVCLASS_SCREEN;
    ALL_BIN_SCREEN;
    SSU_EXTRACT;
    VIRAL_SCREEN;
    CONTIG_CLASSIFY;
    ROUTE_BINS;
    CHECKM1_QC;
    GTDBTK_CLASSIFY;
    SYMCLATRON_CLASSIFY;
    PROKARYOTE_CHARACTERIZATION;
    EUKARYOTE_QC
} from './modules/characterization'
include {
    REPEAT_MASK_EUKARYOTES;
    MAP_RNA_EVIDENCE;
    BRAKER3_EUKARYOTES;
    PRODIGAL_GV_GENES
} from './modules/gene_calling'
include { FUNCTIONAL_ANNOTATION } from './modules/annotation'
include { COLLECT_RECORDS; BUILD_CATALOG; BUILD_REPORT; BUILD_VISUALIZATION } from './modules/catalog'

workflow {
    if (!params.input) {
        error 'Required parameter --input was not provided.'
    }
    if (!(params.run_mode in ['core', 'full'])) {
        error "--run_mode must be core or full, received: ${params.run_mode}"
    }
    ['skip_gene_calling', 'skip_annotation', 'skip_repeat_masking'].each { name ->
        if (!(params[name] instanceof Boolean)) {
            error "--${name} must be a boolean, received: ${params[name]}"
        }
        def prepared = params["prepared_${name}"]
        if (prepared != null && (!(prepared instanceof Boolean) || prepared != params[name])) {
            error "--${name} differs from the prepared run; prepare a new run with the requested options."
        }
    }

    def runId = (params.run_id ?: workflow.runName).toString().replaceAll(/[^A-Za-z0-9_.-]/, '_')
    def samples = channel
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

    def sampleMeta = samples.map { meta, _reads, _assembly, _rna, _proteins, _softmasked -> meta }
    def preflight = channel.empty()
    if (params.run_mode == 'full') {
        FULL_PREFLIGHT(sampleMeta)
        preflight = FULL_PREFLIGHT.out.result
    } else {
        CORE_PREFLIGHT(sampleMeta)
        preflight = CORE_PREFLIGHT.out.result
    }

    def readsBySample = samples.map { meta, reads, _assembly, _rna, _proteins, _softmasked ->
        tuple(meta.sample_id, meta, reads)
    }
    def preflightBySample = preflight.map { meta, stageDir -> tuple(meta.sample_id, stageDir) }
    def qcInput = readsBySample
        .join(preflightBySample)
        .map { _sampleId, meta, reads, stageDir -> tuple(meta, reads, stageDir) }
    READ_QC(qcInput)
    def suppliedAssemblies = samples.map { meta, _reads, assembly, _rna, _proteins, _softmasked ->
        tuple(meta.sample_id, assembly)
    }
    def assemblyInput = READ_QC.out.result
        .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
        .join(suppliedAssemblies)
        .map { _sampleId, meta, readQcDir, suppliedAssembly ->
            tuple(meta, readQcDir, suppliedAssembly)
        }
    ASSEMBLE(assemblyInput)
    def mappingInput = ASSEMBLE.out.result
        .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
        .join(READ_QC.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
        .map { _sampleId, meta, assemblyDir, readQcDir -> tuple(meta, assemblyDir, readQcDir) }
    MAP_DEPTH_QUICKBIN(mappingInput)

    def stageDirs = preflight.mix(
        READ_QC.out.result,
        ASSEMBLE.out.result,
        MAP_DEPTH_QUICKBIN.out.result
    )

    if (params.run_mode == 'full') {
        // Mapping (20 CPUs), SSU (9) and Tiara (3) fill the 32 CPUs after assembly, so the
        // local executor can never hold mapping back behind the side stages.
        SSU_EXTRACT(ASSEMBLE.out.result)
        CONTIG_CLASSIFY(ASSEMBLE.out.result)
        GVCLASS_SCREEN(MAP_DEPTH_QUICKBIN.out.result)
        CHECKEUK_SCREEN(MAP_DEPTH_QUICKBIN.out.result)
        CHECKM2_SCREEN(MAP_DEPTH_QUICKBIN.out.result)
        // The viral screen needs only the assembly. It waits for CheckEUK so that it cannot
        // take the CPUs that mapping or GVClass need; it then uses the CPUs CheckEUK frees.
        def viralInput = ASSEMBLE.out.result
            .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
            .join(CHECKEUK_SCREEN.out.result.map { meta, _dir, _report -> tuple(meta.sample_id, true) })
            .map { _sampleId, meta, stageDir, _checkeukDone -> tuple(meta, stageDir) }
        VIRAL_SCREEN(viralInput)

        def screenInput = MAP_DEPTH_QUICKBIN.out.result
            .map { meta, binningDir -> tuple(meta.sample_id, meta, binningDir) }
            .join(CHECKEUK_SCREEN.out.result.map { meta, dir, report -> tuple(meta.sample_id, dir, report) })
            .join(CHECKM2_SCREEN.out.result.map { meta, dir, report -> tuple(meta.sample_id, dir, report) })
            .join(GVCLASS_SCREEN.out.result.map { meta, dir, report -> tuple(meta.sample_id, dir, report) })
            .map { _sampleId, meta, binningDir, checkeukDir, checkeukReport, checkm2Dir, checkm2Report, gvclassDir, gvclassReport ->
                tuple(meta, binningDir, checkeukDir, checkeukReport, checkm2Dir, checkm2Report, gvclassDir, gvclassReport)
            }
        ALL_BIN_SCREEN(screenInput)

        def routeInput = MAP_DEPTH_QUICKBIN.out.result
            .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
            .join(ALL_BIN_SCREEN.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            .join(SSU_EXTRACT.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            .join(VIRAL_SCREEN.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            .join(CONTIG_CLASSIFY.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            .map { _sampleId, meta, binningDir, screenDir, ssuDir, viralDir, tiaraDir ->
                tuple(meta, binningDir, screenDir, ssuDir, viralDir, tiaraDir)
            }
        ROUTE_BINS(routeInput)

        EUKARYOTE_QC(ROUTE_BINS.out.result)
        CHECKM1_QC(ROUTE_BINS.out.result)
        SYMCLATRON_CLASSIFY(ROUTE_BINS.out.result)
        GTDBTK_CLASSIFY(ROUTE_BINS.out.result)

        def prokaryoteInput = ROUTE_BINS.out.result
            .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
            .join(CHECKM1_QC.out.result.map { meta, dir, report -> tuple(meta.sample_id, dir, report) })
            .join(GTDBTK_CLASSIFY.out.result.map { meta, dir, report -> tuple(meta.sample_id, dir, report) })
            .join(SYMCLATRON_CLASSIFY.out.result.map { meta, dir, report -> tuple(meta.sample_id, dir, report) })
            .map { _sampleId, meta, routingDir, checkm1Dir, checkm1Report, gtdbtkDir, gtdbtkReport, symclatronDir, symclatronReport ->
                tuple(meta, routingDir, checkm1Dir, checkm1Report, gtdbtkDir, gtdbtkReport, symclatronDir, symclatronReport)
            }
        PROKARYOTE_CHARACTERIZATION(prokaryoteInput)

        def geneStages = channel.empty()
        def annotationStages = channel.empty()
        if (!params.skip_gene_calling) {
            def geneRouting = ROUTE_BINS.out.result
                .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                .join(PROKARYOTE_CHARACTERIZATION.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
                .map { sampleId, meta, routingDir, _prokaryoteDir -> tuple(sampleId, meta, routingDir) }
            def repeatStages = channel.empty()
            def eukaryoticGenomes = geneRouting.map { _sampleId, meta, routingDir -> tuple(meta, routingDir) }
            if (!params.skip_repeat_masking) {
                def suppliedSoftmasked = samples.map { meta, _reads, _assembly, _rna, _proteins, softmasked ->
                    tuple(meta.sample_id, softmasked)
                }
                def maskingInput = geneRouting
                    .join(suppliedSoftmasked)
                    .map { _sampleId, meta, routingDir, softmasked -> tuple(meta, routingDir, softmasked) }
                REPEAT_MASK_EUKARYOTES(maskingInput)
                repeatStages = REPEAT_MASK_EUKARYOTES.out.result
                eukaryoticGenomes = repeatStages
            }

            def rnaInputBase = ASSEMBLE.out.result
                .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                .join(eukaryoticGenomes.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            def suppliedRna = samples.map { meta, _reads, _assembly, rna, _proteins, _softmasked ->
                tuple(meta.sample_id, rna)
            }
            def rnaInput = rnaInputBase
                .join(suppliedRna)
                .map { _sampleId, meta, assemblyDir, maskingDir, rna -> tuple(meta, assemblyDir, maskingDir, rna) }
            MAP_RNA_EVIDENCE(rnaInput)

            def brakerInputBase = eukaryoticGenomes
                .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                .join(MAP_RNA_EVIDENCE.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
            def suppliedProteins = samples.map { meta, _reads, _assembly, _rna, proteins, _softmasked ->
                tuple(meta.sample_id, proteins)
            }
            def brakerInput = brakerInputBase
                .join(suppliedProteins)
                .map { _sampleId, meta, maskingDir, rnaDir, proteins -> tuple(meta, maskingDir, rnaDir, proteins) }
            BRAKER3_EUKARYOTES(brakerInput)

            def prodigalInput = geneRouting
                .join(VIRAL_SCREEN.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
                .map { _sampleId, meta, routingDir, viralDir -> tuple(meta, routingDir, viralDir) }
            PRODIGAL_GV_GENES(prodigalInput)

            geneStages = repeatStages.mix(
                MAP_RNA_EVIDENCE.out.result,
                BRAKER3_EUKARYOTES.out.result,
                PRODIGAL_GV_GENES.out.result
            )
            if (!params.skip_annotation) {
                def annotationInput = PRODIGAL_GV_GENES.out.result
                    .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
                    .join(BRAKER3_EUKARYOTES.out.result.map { meta, stageDir -> tuple(meta.sample_id, stageDir) })
                    .map { _sampleId, meta, prodigalDir, brakerDir -> tuple(meta, prodigalDir, brakerDir) }
                FUNCTIONAL_ANNOTATION(annotationInput)
                annotationStages = FUNCTIONAL_ANNOTATION.out.result
            }
        }

        stageDirs = stageDirs.mix(
            ALL_BIN_SCREEN.out.result,
            SSU_EXTRACT.out.result,
            VIRAL_SCREEN.out.result,
            CONTIG_CLASSIFY.out.result,
            ROUTE_BINS.out.result,
            PROKARYOTE_CHARACTERIZATION.out.result,
            EUKARYOTE_QC.out.result,
            geneStages,
            annotationStages
        )
    }

    def groupedStages = stageDirs
        .map { meta, stageDir -> tuple(meta.sample_id, meta, stageDir) }
        .groupTuple(by: 0)
        .map { _sampleId, metas, dirs -> tuple(metas[0], dirs.sort(false) { dir -> dir.name }) }
    COLLECT_RECORDS(groupedStages)

    def allRecords = COLLECT_RECORDS.out.result
        .map { _meta, records -> records }
        .collect()
    BUILD_CATALOG(allRecords)
    BUILD_REPORT(BUILD_CATALOG.out.result)
    BUILD_VISUALIZATION(BUILD_CATALOG.out.result)
}
