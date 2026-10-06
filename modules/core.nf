process READ_QC {
    tag "${meta.sample_id}"
    label 'qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(preflight_dir)

    output:
    tuple val(meta), path('01_read_qc'), emit: result

    script:
    def compressionThreads = Math.max(task.cpus - 2, 1)
    def minLength = meta.platform == 'ont' ? params.ont_min_length : params.hifi_min_length
    def filterCommand = meta.platform == 'ont' \
        ? "chopper -q ${params.ont_min_quality} -l ${minLength} --threads 1" \
        : "seqkit seq -j 1 --min-len ${minLength}"
    """
    mkdir -p 01_read_qc

    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        seqkit stats -j 1 -T "${reads}" > 01_read_qc/raw_stats.tsv

    if [[ "${reads}" == *.gz ]]; then
        pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            pigz -p 1 -dc "${reads}"
    else
        command cat "${reads}"
    fi | pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        ${filterCommand} \
      | pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        pigz -p ${compressionThreads} > 01_read_qc/filtered.fastq.gz

    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        seqkit stats -j 1 -T 01_read_qc/filtered.fastq.gz > 01_read_qc/filtered_stats.tsv
    [[ \$(awk 'NR == 2 {print \$4}' 01_read_qc/filtered_stats.tsv) -gt 0 ]]

    TOOL_VERSION=\$(pixi run --as-is --manifest-path "${params.core_manifest}" \
        --environment core ${meta.platform == 'ont' ? 'chopper --version' : 'seqkit version'} \
        | tr '\n' ' ' | sed 's/"/\\\\"/g')
    cat > 01_read_qc/stage.json <<JSON
    {
      "stage": "read_qc",
      "status": "completed",
      "tool": "${meta.platform == 'ont' ? 'chopper' : 'seqkit'}",
      "tool_version": "\${TOOL_VERSION}",
      "database_name": null,
      "database_version": null,
      "command": "${filterCommand}",
      "outputs": {
        "reads": "filtered.fastq.gz",
        "raw_stats": "raw_stats.tsv",
        "filtered_stats": "filtered_stats.tsv"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 01_read_qc
    : > 01_read_qc/filtered.fastq.gz
    printf 'file\tformat\ttype\tnum_seqs\tsum_len\nreads.fastq\tFASTQ\tDNA\t1\t4\n' > 01_read_qc/raw_stats.tsv
    cp 01_read_qc/raw_stats.tsv 01_read_qc/filtered_stats.tsv
    printf '%s\n' '{"stage":"read_qc","status":"completed","tool":"stub","tool_version":"stub","database_version":null,"command":"stub","outputs":{"reads":"filtered.fastq.gz","raw_stats":"raw_stats.tsv","filtered_stats":"filtered_stats.tsv"}}' > 01_read_qc/stage.json
    """
}

process ASSEMBLE {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(read_qc_dir), path(assembly_input)

    output:
    tuple val(meta), path('02_assembly'), emit: result

    script:
    def readModeFlag = meta.platform == 'pacbio_hifi' ? '--hifi' : '--nano-r10'
    def source = assembly_input ? assembly_input.toString() : ''
    """
    mkdir -p 02_assembly
    if [[ -n "${source}" ]]; then
        cp -L "${source}" 02_assembly/assembly_primary.fa
        TOOL_VERSION="supplied"
        COMMAND="validated supplied assembly"
    else
        pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            myloasm "${read_qc_dir}/filtered.fastq.gz" -o myloasm_output \
            -t ${task.cpus} ${readModeFlag}
        [[ -s myloasm_output/assembly_primary.fa ]]
        cp myloasm_output/assembly_primary.fa 02_assembly/assembly_primary.fa
        if [[ -s myloasm_output/final_contig_graph.gfa ]]; then
            cp myloasm_output/final_contig_graph.gfa 02_assembly/final_contig_graph.gfa
        fi
        TOOL_VERSION=\$(pixi run --as-is --manifest-path "${params.core_manifest}" \
            --environment core myloasm --version | tr '\n' ' ' | sed 's/"/\\\\"/g')
        COMMAND="myloasm reads -o output -t ${task.cpus} ${readModeFlag}"
    fi

    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        seqkit stats -T 02_assembly/assembly_primary.fa > 02_assembly/assembly_stats.tsv
    [[ \$(awk 'NR == 2 {print \$4}' 02_assembly/assembly_stats.tsv) -gt 0 ]]
    cat > 02_assembly/stage.json <<JSON
    {
      "stage": "assembly",
      "status": "completed",
      "tool": "${source ? 'supplied' : 'myloasm'}",
      "tool_version": "\${TOOL_VERSION}",
      "database_name": null,
      "database_version": null,
      "command": "\${COMMAND}",
      "outputs": {
        "assembly": "assembly_primary.fa",
        "stats": "assembly_stats.tsv"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 02_assembly
    printf '>stub_contig\nACGTACGTACGT\n' > 02_assembly/assembly_primary.fa
    printf 'file\tformat\ttype\tnum_seqs\tsum_len\nassembly_primary.fa\tFASTA\tDNA\t1\t12\n' > 02_assembly/assembly_stats.tsv
    printf '%s\n' '{"stage":"assembly","status":"completed","tool":"stub","tool_version":"stub","database_version":null,"command":"stub","outputs":{"assembly":"assembly_primary.fa","stats":"assembly_stats.tsv"}}' > 02_assembly/stage.json
    """
}

process MAP_DEPTH_QUICKBIN {
    tag "${meta.sample_id}"
    label 'mapping'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly_dir), path(read_qc_dir)

    output:
    tuple val(meta), path('03_binning'), emit: result

    script:
    def mapPreset = meta.platform == 'pacbio_hifi' ? 'map-hifi' : 'map-ont'
    def mapThreads = Math.max(task.cpus - 4, 1)
    def samtoolsThreads = Math.max(task.cpus - 1, 1)
    def heapGb = Math.max(task.memory.toGiga().intValue() - 16, 4)
    """
    [[ -s "${params.quickclade_ref}" ]]
    mkdir -p 03_binning/bins
    ASSEMBLY="${assembly_dir}/assembly_primary.fa"
    READS="${read_qc_dir}/filtered.fastq.gz"

    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        minimap2 -ax ${mapPreset} --secondary=no -t ${mapThreads} "\${ASSEMBLY}" "\${READS}" \
        2> 03_binning/minimap2.log \
      | pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        samtools view -b -F 2308 -q ${params.mapping_min_mapq} - \
      | pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        samtools sort -@ 1 -m 4G -T sort.tmp \
        -o 03_binning/primary.q20.sorted.bam -
    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        samtools index -@ ${samtoolsThreads} 03_binning/primary.q20.sorted.bam
    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        samtools quickcheck -v 03_binning/primary.q20.sorted.bam
    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        samtools flagstat -@ ${samtoolsThreads} 03_binning/primary.q20.sorted.bam \
        > 03_binning/flagstat.txt
    pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
        samtools coverage -d 0 -o 03_binning/coverage.tsv \
        03_binning/primary.q20.sorted.bam

    export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=${task.cpus}"
    pixi run --as-is --manifest-path "${params.bbtools_manifest}" --environment bbtools \
        quickbin.sh \
        in="\${ASSEMBLY}" reads=03_binning/primary.q20.sorted.bam \
        out=03_binning/bins/ covout=03_binning/quickbin_coverage.tsv \
        report=03_binning/quickbin_report.tsv chaff clade=f callssu=t normal \
        minid=${params.quickbin_min_identity} minmapq=${params.mapping_min_mapq} \
        threads=${task.cpus} -Xmx${heapGb}g
    [[ -s 03_binning/quickbin_report.tsv ]]
    [[ -s 03_binning/quickbin_coverage.tsv ]]
    BIN_COUNT=\$(find 03_binning/bins -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l)

    MYLOASM_MAP_VERSION=\$(pixi run --as-is --manifest-path "${params.core_manifest}" \
        --environment core minimap2 --version)
    BBTOOLS_VERSION=\$(pixi run --as-is --manifest-path "${params.bbtools_manifest}" \
        --environment bbtools bbversion.sh)
    [[ -n "\${BBTOOLS_VERSION}" ]]
    if [[ \${BIN_COUNT} -gt 0 ]]; then
        pixi run --as-is --manifest-path "${params.bbtools_manifest}" --environment bbtools \
            quickclade.sh \
            03_binning/bins/*.fa ref="${params.quickclade_ref}" server=f \
            out=03_binning/quickclade.tsv composition=03_binning/quickclade_composition.tsv \
            format=machine fast records=1 color=f showloading=f \
            threads=${task.cpus} -Xmx${heapGb}g
        [[ -s 03_binning/quickclade.tsv ]]
        [[ -s 03_binning/quickclade_composition.tsv ]]
        QUICKCLADE_SHA256=\$(sha256sum "${params.quickclade_ref}" | awk '{print \$1}')
        cat > 03_binning/stage.json <<JSON
    {
      "stage": "binning",
      "status": "completed",
      "tool": "minimap2;samtools;QuickBin;QuickClade",
      "tool_version": "minimap2 \${MYLOASM_MAP_VERSION}; BBTools \${BBTOOLS_VERSION}",
      "database_name": "quickclade_ref",
      "database_version": "${params.quickclade_ref_version}",
      "database_sha256": "\${QUICKCLADE_SHA256}",
      "command": "minimap2 -ax ${mapPreset}; samtools view -F 2308 -q ${params.mapping_min_mapq}; QuickBin unpaired long-read BAM contributes coverage without mate edges; quickbin minid=${params.quickbin_min_identity}; quickclade ref=local server=f",
      "outputs": {
        "bam": "primary.q20.sorted.bam",
        "bam_index": "primary.q20.sorted.bam.bai",
        "flagstat": "flagstat.txt",
        "coverage": "coverage.tsv",
        "bins": "bins",
        "quickbin_report": "quickbin_report.tsv",
        "quickbin_coverage": "quickbin_coverage.tsv",
        "quickclade": "quickclade.tsv",
        "quickclade_composition": "quickclade_composition.tsv"
      }
    }
JSON
    else
        cat > 03_binning/stage.json <<JSON
    {
      "stage": "binning",
      "status": "completed",
      "tool": "minimap2;samtools;QuickBin",
      "tool_version": "minimap2 \${MYLOASM_MAP_VERSION}; BBTools \${BBTOOLS_VERSION}",
      "database_name": null,
      "database_version": null,
      "command": "minimap2 -ax ${mapPreset}; samtools view -F 2308 -q ${params.mapping_min_mapq}; QuickBin found no accepted bins; QuickClade not run because no bins",
      "outputs": {
        "bam": "primary.q20.sorted.bam",
        "bam_index": "primary.q20.sorted.bam.bai",
        "flagstat": "flagstat.txt",
        "coverage": "coverage.tsv",
        "bins": "bins",
        "quickbin_report": "quickbin_report.tsv",
        "quickbin_coverage": "quickbin_coverage.tsv"
      }
    }
JSON
    fi
    """

    stub:
    """
    mkdir -p 03_binning/bins
    touch 03_binning/primary.q20.sorted.bam 03_binning/primary.q20.sorted.bam.bai
    printf 'stub\n' > 03_binning/flagstat.txt
    printf '#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\n' > 03_binning/coverage.tsv
    printf '>stub_contig\nACGTACGTACGT\n' > 03_binning/bins/bin_1.fa
    printf '#Num\tFile\tSize\tContigs\tGC\tDepth\tMinDepth\tMaxDepth\tTaxID\tLineage\tContig\n1\tbin_1.fa\t12\t1\t0.5\t1\t1\t1\t0\t\tstub_contig\n' > 03_binning/quickbin_report.tsv
    printf 'ID\tAvg_fold\n' > 03_binning/quickbin_coverage.tsv
    printf '#QueryName\tQ_Bases\tlineage\tConfLevel\tConfidence\nbin_1\t12\t\t0\t0\n' > 03_binning/quickclade.tsv
    printf 'query\tcomposition\n' > 03_binning/quickclade_composition.tsv
    printf '%s\n' '{"stage":"binning","status":"completed","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"bam":"primary.q20.sorted.bam","bam_index":"primary.q20.sorted.bam.bai","flagstat":"flagstat.txt","coverage":"coverage.tsv","bins":"bins","quickbin_report":"quickbin_report.tsv","quickbin_coverage":"quickbin_coverage.tsv","quickclade":"quickclade.tsv","quickclade_composition":"quickclade_composition.tsv"}}' > 03_binning/stage.json
    """
}
