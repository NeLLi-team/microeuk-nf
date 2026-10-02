process PROKARYOTE_CHARACTERIZATION {
    tag "${meta.sample_id}"
    label 'heavy_qc'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(routing_dir)

    output:
    tuple val(meta), path('08_prokaryote_qc'), emit: result

    script:
    """
    mkdir -p 08_prokaryote_qc/{checkm1,checkm2,gtdbtk,symclatron}
    BIN_COUNT=\$(find "${routing_dir}/prokaryotic" -maxdepth 1 -type f -name '*.fna' -size +0c | wc -l)
    STATUS=skipped
    REASON=no_routed_prokaryotic_bins
    if [[ \${BIN_COUNT} -gt 0 ]]; then
        STATUS=completed
        REASON=
        export CHECKM_DATA_PATH="${params.checkm1_db}"
        TMPDIR=/tmp pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment checkm1 \
            checkm lineage_wf --reduced_tree --extension fna --threads ${task.cpus} \
            --pplacer_threads 8 --tab_table --file "\$PWD/08_prokaryote_qc/checkm1/quality_report.tsv" \
            --tmpdir /tmp "${routing_dir}/prokaryotic" \
            "\$PWD/08_prokaryote_qc/checkm1/output"
        [[ -s 08_prokaryote_qc/checkm1/quality_report.tsv ]]

        export GTDBTK_DATA_PATH="${params.gtdbtk_db}"
        TMPDIR=/tmp pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment gtdbtk232 \
            gtdbtk classify_wf --genome_dir "${routing_dir}/prokaryotic" \
            --out_dir 08_prokaryote_qc/gtdbtk --extension fna \
            --cpus ${task.cpus} --pplacer_cpus 8 --tmpdir /tmp
        [[ \$(find 08_prokaryote_qc/gtdbtk -type f \
            -name 'gtdbtk.*.summary.tsv' -size +0c | wc -l) -gt 0 ]]

        pixi run --as-is --manifest-path "${params.checkm_manifest}" --environment symclatron \
            symclatron classify --genome-dir "${routing_dir}/prokaryotic" \
            --input-kind contigs --input-ext .fna --threads ${task.cpus} \
            --confidence-threshold 0.725 --output-dir 08_prokaryote_qc/symclatron
        [[ -s 08_prokaryote_qc/symclatron/symclatron_results.tsv ]]
        [[ \$(awk 'END {print NR - 1}' 08_prokaryote_qc/symclatron/symclatron_results.tsv) -eq \${BIN_COUNT} ]]
    fi

    CHECKM1_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
        --environment checkm1 python -c \
        'from pathlib import Path; import checkm; print((Path(checkm.__file__).parent / "VERSION").read_text().splitlines()[0])')
    GTDBTK_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
        --environment gtdbtk232 gtdbtk --version)
    GTDBTK_VERSION=\$(printf '%s\n' "\${GTDBTK_VERSION}" | sed 's/"/\\\\"/g')
    SYMCLATRON_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.checkm_manifest}" \
        --environment symclatron symclatron --version)
    SYMCLATRON_VERSION=\$(printf '%s\n' "\${SYMCLATRON_VERSION}" | sed 's/"/\\\\"/g')
    [[ -n "\${CHECKM1_VERSION}" && -n "\${GTDBTK_VERSION}" && -n "\${SYMCLATRON_VERSION}" ]]

    cat > 08_prokaryote_qc/stage.json <<JSON
    {
      "stage": "prokaryote_qc",
      "status": "\${STATUS}",
      "tool": "CheckM1;GTDB-Tk;Symclatron",
      "tool_version": "CheckM1 \${CHECKM1_VERSION}; GTDB-Tk \${GTDBTK_VERSION}; Symclatron \${SYMCLATRON_VERSION}",
      "database_name": "checkm1_db;gtdbtk_db;symclatron_bundled",
      "database_version": "CheckM1 ${params.checkm1_db_version}; GTDB-Tk ${params.gtdbtk_db_version}; Symclatron bundled",
      "command": "bounded per-sample batch over routed prokaryotic bins",
      "reason": "\${REASON}",
      "tools": {
        "checkm1": {"version": "\${CHECKM1_VERSION}", "database_name": "checkm1_db", "database_version": "${params.checkm1_db_version}", "output": "checkm1/quality_report.tsv"},
        "gtdbtk": {"version": "\${GTDBTK_VERSION}", "database_name": "gtdbtk_db", "database_version": "${params.gtdbtk_db_version}", "output": "gtdbtk"},
        "symclatron": {"version": "\${SYMCLATRON_VERSION}", "database_name": "symclatron_bundled", "database_version": "bundled", "output": "symclatron/symclatron_results.tsv"}
      },
      "outputs": {
        "checkm1": "checkm1/quality_report.tsv",
        "gtdbtk": "gtdbtk",
        "symclatron": "symclatron/symclatron_results.tsv"
      }
    }
JSON
    """

    stub:
    """
    mkdir -p 08_prokaryote_qc/{checkm1,gtdbtk,symclatron}
    printf 'Bin Id\tMarker lineage\n' > 08_prokaryote_qc/checkm1/quality_report.tsv
    printf 'taxon_oid\tcompleteness_UNI56\tclassification\tconfidence\tpasses_confidence_threshold\tclassification_thresholded\n' > 08_prokaryote_qc/symclatron/symclatron_results.tsv
    printf '%s\n' '{"stage":"prokaryote_qc","status":"skipped","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"checkm1":"checkm1/quality_report.tsv","gtdbtk":"gtdbtk","symclatron":"symclatron/symclatron_results.tsv"},"reason":"no_routed_prokaryotic_bins"}' > 08_prokaryote_qc/stage.json
    """
}

process REPEAT_MASK_EUKARYOTES {
    tag "${meta.sample_id}"
    label 'gene'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'
    scratch '/tmp'

    input:
    tuple val(meta), path(routing_dir), path(softmasked_input)

    output:
    tuple val(meta), path('09_repeat_masking'), emit: result

    script:
    def suppliedSoftmasked = softmasked_input ? softmasked_input.toString() : ''
    def repeatMaskerWorkers = Math.max((task.cpus as int).intdiv(4), 1)
    """
    set -o pipefail
    # Nextflow keeps .command.sh in the persistent task directory when using scratch.
    REPEATMASK_WORK_DIR=\$(dirname "\$(realpath "\${BASH_SOURCE[0]}")")
    REPEATMASK_LOG_DIR="\$PWD/09_repeat_masking/logs"
    preserve_failure_logs() {
        local status=\$?
        trap - EXIT
        if [[ \${status} -ne 0 && -d "\${REPEATMASK_LOG_DIR}" ]]; then
            local failure_dir
            if failure_dir=\$(mkdir -p "\${REPEATMASK_WORK_DIR}/repeatmask_failure_logs" \
                && mktemp -d "\${REPEATMASK_WORK_DIR}/repeatmask_failure_logs/attempt.XXXXXXXX"); then
                cp -a "\${REPEATMASK_LOG_DIR}/." "\${failure_dir}/" \
                    || printf 'Could not preserve repeat-masking logs in %s.\n' "\${failure_dir}" >&2
            else
                printf '%s\n' 'Could not create the repeat-masking failure log directory.' >&2
            fi
        fi
        exit "\${status}"
    }
    trap preserve_failure_logs EXIT
    mkdir -p 09_repeat_masking/{masked,models,logs}
    printf 'bin_id\tmask_source\trepeat_library\trepeat_family_count\tmasked_base_count\trepeatmasker_parallel_jobs\tdfam_release\n' \
        > 09_repeat_masking/logs/masking.tsv
    canonical_sha256() {
        pixi run --as-is --quiet --manifest-path "${params.core_manifest}" --environment core \
            seqkit sort -j 1 "\$1" \
          | pixi run --as-is --quiet --manifest-path "${params.core_manifest}" --environment core \
            seqkit seq -j 1 -i -u -w 0 \
          | sha256sum | awk '{print \$1}'
    }
    lowercase_base_count() {
        awk '!/^>/ {line=\$0; gsub(/[^a-z]/, "", line); count += length(line)}
             END {print count + 0}' "\$1"
    }
    reject_native_child_crashes() {
        local log=\$1
        if grep -Eiq 'Aborted|Segmentation fault|core dumped' "\${log}"; then
            printf 'Native child-process crash marker found in %s.\n' "\${log}" >&2
            exit 3
        fi
    }
    shopt -s nullglob
    BINS=("${routing_dir}"/eukaryotic/*.fna)
    STATUS=skipped
    REASON=no_routed_eukaryotic_bins
    RAN_REPEATMODELER=0
    RAN_REPEATMASKER=0
    RAN_RECON_ONLY=0
    if [[ \${#BINS[@]} -gt 0 ]]; then
        STATUS=completed
        REASON=
        for FASTA in "\${BINS[@]}"; do
            BIN=\$(basename "\${FASTA}" .fna)
            if [[ -n "${suppliedSoftmasked}" ]]; then
                [[ -s "${suppliedSoftmasked}" ]]
                pixi run --as-is --quiet --manifest-path "${params.core_manifest}" --environment core \
                    seqkit seq -j 1 -n -i "\${FASTA}" > ids.txt
                pixi run --as-is --quiet --manifest-path "${params.core_manifest}" --environment core \
                    seqkit grep -j 1 -f ids.txt "${suppliedSoftmasked}" \
                    > "09_repeat_masking/masked/\${BIN}.softmasked.fna"
                MASKED_OUT="09_repeat_masking/masked/\${BIN}.softmasked.fna"
                [[ \$(grep -c '^>' "\${MASKED_OUT}") -eq \$(wc -l < ids.txt) ]]
                [[ \$(canonical_sha256 "\${FASTA}") == \$(canonical_sha256 "\${MASKED_OUT}") ]]
                MASKED_BASES=\$(lowercase_base_count "\${MASKED_OUT}")
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "\${BIN}" supplied not_applicable 0 "\${MASKED_BASES}" 0 not_applicable \
                    >> 09_repeat_masking/logs/masking.tsv
                continue
            fi

            RAN_REPEATMODELER=1
            mkdir -p "09_repeat_masking/models/\${BIN}"
            FASTA_ABS=\$(realpath "\${FASTA}")
            pushd "09_repeat_masking/models/\${BIN}" >/dev/null
            pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" --environment repeatmasker \
                BuildDatabase -name repeatdb "\${FASTA_ABS}" \
                > "../../logs/\${BIN}.builddatabase.log" 2>&1
            REPEATMODELER_LOG="../../logs/\${BIN}.repeatmodeler.log"
            REPEATMODELER_STATUS=0
            RECON_ONLY=0
            pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" --environment repeatmasker \
                RepeatModeler -database repeatdb -threads ${task.cpus} \
                > "\${REPEATMODELER_LOG}" 2>&1 || REPEATMODELER_STATUS=\$?
            printf '%s\n' "\${REPEATMODELER_STATUS}" > "../../logs/\${BIN}.repeatmodeler.initial.exitcode"
            if [[ \${REPEATMODELER_STATUS} -ne 0 ]]; then
                reject_native_child_crashes "\${REPEATMODELER_LOG}"
                mapfile -t SCOUT_SAMPLES < <(find . -mindepth 3 -maxdepth 3 -type f \
                    -path './RM_*/round-1/sampleDB-1.fa' -size +0c -print)
                if [[ \${REPEATMODELER_STATUS} -eq 2 && \${#SCOUT_SAMPLES[@]} -eq 1 ]]; then
                    SCOUT_ROUND=\$(dirname "\$(realpath "\${SCOUT_SAMPLES[0]}")")
                    FILTERED="\${SCOUT_SAMPLES[0]}.rscons.filtered"
                    RANGES="\${SCOUT_SAMPLES[0]}.rscons-ranges.tsv"
                    FILTERED_COUNT=\$(awk '
                        /^   - RepeatScout: Running filtering stage[.][.] [1-9][0-9]* families remaining\$/ {
                            count++; families=\$6
                        }
                        END {if (count == 1) print families}
                    ' "\${REPEATMODELER_LOG}")
                    # RepeatModeler 2.0.9 skips RS families below five parsed instances,
                    # then searches the empty refined-cons.fa without checking it.
                    if [[ -z "\${FILTERED_COUNT}" || ! -s "\${FILTERED}" \
                        || ! -f "\${SCOUT_ROUND}/refined-cons.fa" || -s "\${SCOUT_ROUND}/refined-cons.fa" ]] \
                        || ! cmp -s "\${FILTERED}" "\${SCOUT_ROUND}/consensi.fa" \
                        || ! grep -Eq '^   - Refining 0 families[.][.][.] [0-9:]+ [(]hh:mm:ss[)] Elapsed Time\$' "\${REPEATMODELER_LOG}" \
                        || ! grep -Fxq "   - Redundant Families and Large Satellite Filtering.. NCBIBlastSearchEngine::search: Error...compressed subject database (\${SCOUT_ROUND}/refined-cons.fa) does not exist!" "\${REPEATMODELER_LOG}" \
                        || ! awk -F '\\t' -v expected="\${FILTERED_COUNT}" '
                            FNR == NR {
                                if (/^>/) {
                                    id=substr(\$0, 2); sub(/[[:space:]].*\$/, "", id)
                                    if (id !~ /^R=[0-9]+\$/ || id in families) bad=1
                                    families[id]=0; total++
                                }
                                next
                            }
                            {
                                if (NF != 8 || \$1 == "" || \$1 ~ /[[:space:]]/ || \$2 !~ /^-?[0-9]+\$/ \
                                    || \$3 !~ /^[0-9]+\$/ || \$4 !~ /^R=[0-9]+\$/ \
                                    || \$5 !~ /^[+-]\$/ || \$6 !~ /^[0-9]+\$/ \
                                    || \$7 !~ /^[0-9]+\$/ || \$8 !~ /^[0-9]+\$/) bad=1
                                # The native unsigned-coordinate parser ignores negative starts.
                                if (\$4 in families) seen[\$4]++
                                if (\$4 in families && \$0 ~ /^[^[:space:]]+\\t[0-9]+\\t[0-9]+/) families[\$4]++
                            }
                            END {
                                for (id in families) if (!seen[id] || families[id] >= 5) bad=1
                                exit (bad || total != expected)
                            }
                        ' "\${FILTERED}" "\${RANGES}"; then
                        cat "\${REPEATMODELER_LOG}" >&2
                        printf '%s\n' 'RepeatScout failure did not confirm empty refinement below the five-instance cutoff.' >&2
                        exit "\${REPEATMODELER_STATUS}"
                    fi
                    RECOVERY_REASON='Confirmed RepeatScout empty refinement: all filtered families have fewer than five native-parsed instances'
                elif [[ \${REPEATMODELER_STATUS} -eq 1 ]] \
                    && grep -Fxq 'build_lmer_table failed. Exit code 256' "\${REPEATMODELER_LOG}"; then
                    mapfile -t LMER_SETTINGS < <(awk '
                        /^   - RepeatScout: Running build_lmer_table [(] l = [0-9]+, min = [0-9]+ [)][.][.]\$/ {
                            gsub(",", "", \$8); print \$8, \$11
                        }' "\${REPEATMODELER_LOG}")
                    if [[ \${#SCOUT_SAMPLES[@]} -ne 1 || \${#LMER_SETTINGS[@]} -ne 1 ]]; then
                        printf '%s\n' 'Cannot identify one failed RepeatScout sample and its native l/min settings.' >&2
                        exit "\${REPEATMODELER_STATUS}"
                    fi
                    read -r LMER_LENGTH LMER_MINIMUM <<< "\${LMER_SETTINGS[0]}"
                    PROBE="../../logs/\${BIN}.repeatscout-probe"
                    cp "\${SCOUT_SAMPLES[0]}" "\${PROBE}.fa"
                    printf 'build_lmer_table -l %s -min %s -sequence %s.fa -freq %s.lfreq\n' \
                        "\${LMER_LENGTH}" "\${LMER_MINIMUM}" "\${PROBE}" "\${PROBE}" > "\${PROBE}.command.txt"
                    PROBE_STATUS=0
                    pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" --environment repeatmasker \
                        build_lmer_table -l "\${LMER_LENGTH}" -min "\${LMER_MINIMUM}" \
                        -sequence "\${PROBE}.fa" -freq "\${PROBE}.lfreq" \
                        > "\${PROBE}.stdout" 2> "\${PROBE}.stderr" || PROBE_STATUS=\$?
                    printf '%s\n' "\${PROBE_STATUS}" > "\${PROBE}.exitcode"
                    tr '\\r' '\\n' < "\${PROBE}.stderr" > "\${PROBE}.stderr.normalized"
                    reject_native_child_crashes "\${PROBE}.stdout"
                    reject_native_child_crashes "\${PROBE}.stderr"
                    if [[ \${PROBE_STATUS} -ne 1 || -e "\${PROBE}.lfreq" ]] \
                        || ! grep -Fxq 'OOPS no good lmers' "\${PROBE}.stderr.normalized"; then
                        cat "\${PROBE}.stderr" >&2
                        printf '%s\n' 'RepeatScout failure did not confirm the no-seed condition.' >&2
                        exit "\${REPEATMODELER_STATUS}"
                    fi
                    RECOVERY_REASON='Confirmed RepeatScout no-seed outcome'
                else
                    cat "\${REPEATMODELER_LOG}" >&2
                    exit "\${REPEATMODELER_STATUS}"
                fi
                mkdir failed_repeatscout
                mv "\${SCOUT_SAMPLES[0]%%/round-1/*}" failed_repeatscout/
                mv "\${REPEATMODELER_LOG}" "../../logs/\${BIN}.repeatmodeler.initial.log"
                RECON_ONLY=1
                RAN_RECON_ONLY=1
                printf '%s; retrying RepeatModeler -skipRS (RECON-only discovery).\n' "\${RECOVERY_REASON}" \
                    > "../../logs/\${BIN}.repeatmodeler.recovery.log"
                REPEATMODELER_STATUS=0
                pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" --environment repeatmasker \
                    RepeatModeler -database repeatdb -threads ${task.cpus} -skipRS \
                    > "\${REPEATMODELER_LOG}" 2>&1 || REPEATMODELER_STATUS=\$?
            fi
            printf '%s\n' "\${REPEATMODELER_STATUS}" > "../../logs/\${BIN}.repeatmodeler.exitcode"
            if [[ \${REPEATMODELER_STATUS} -ne 0 ]]; then
                cat "\${REPEATMODELER_LOG}" >&2
                exit "\${REPEATMODELER_STATUS}"
            fi
            reject_native_child_crashes "\${REPEATMODELER_LOG}"
            RAW_CONSENSUS=\$(find . -maxdepth 2 -type f -name 'consensi.fa' -size +0c -print -quit)
            CONSENSUS=\$(find . -maxdepth 2 -type f -name 'consensi.fa.classified' -size +0c -print -quit)
            if [[ -n "\${RAW_CONSENSUS}" && -n "\${CONSENSUS}" ]]; then
                if [[ ${task.cpus} -lt 4 ]]; then
                    printf '%s\n' 'RepeatMasker with RMBlast requires at least 4 allocated CPUs.' >&2
                    exit 1
                fi
                REPEAT_FAMILIES=\$(grep -c '^>' "\${CONSENSUS}")
                [[ \${REPEAT_FAMILIES} -gt 0 ]]
                ID_MAP=idmap.tsv
                MASK_QUERY=repeatmasker-query.fna
                pixi run --as-is --quiet --manifest-path "${params.core_manifest}" --environment core \
                    seqkit seq -j 1 -n -i "\${FASTA_ABS}" \
                  | awk '{printf "rm%d\\t%s\\n", NR, \$0}' > "\${ID_MAP}"
                pixi run --as-is --quiet --manifest-path "${params.core_manifest}" --environment core \
                    seqkit replace -j 1 -p '^.*\$' -r 'rm{nr}' "\${FASTA_ABS}" \
                    > "\${MASK_QUERY}"
                mkdir -p masked
                REPEATMASKER_LOG="../../logs/\${BIN}.repeatmasker.log"
                printf '%s\n' \
                    'Header contract: RepeatMasker query uses exact >rmN headers without descriptions; published de novo output restores canonical first-token identifiers without descriptions.' \
                    > "\${REPEATMASKER_LOG}"
                pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" --environment repeatmasker \
                    RepeatMasker -pa ${repeatMaskerWorkers} -lib "\${CONSENSUS}" -xsmall \
                    -dir masked "\${MASK_QUERY}" >> "\${REPEATMASKER_LOG}" 2>&1
                reject_native_child_crashes "\${REPEATMASKER_LOG}"
                mapfile -t MASKED_FILES < <(
                    find masked -maxdepth 1 -type f -name '*.masked' -size +0c -print
                )
                [[ \${#MASKED_FILES[@]} -eq 1 ]]
                awk '
                    FNR == NR {
                        separator = index(\$0, "\\t")
                        original[substr(\$0, 1, separator - 1)] = substr(\$0, separator + 1)
                        next
                    }
                    /^>/ {
                        alias = substr(\$0, 2)
                        if (!(alias in original)) {
                            printf "Unknown RepeatMasker query identifier: %s\\n", alias > "/dev/stderr"
                            exit 3
                        }
                        print ">" original[alias]
                        next
                    }
                    {print}
                ' "\${ID_MAP}" "\${MASKED_FILES[0]}" \
                    > "../../masked/\${BIN}.softmasked.fna"
                MASK_SOURCE=de_novo_repeatmasker
                REPEAT_LIBRARY=consensi.fa.classified
                MASKER_JOBS=${repeatMaskerWorkers}
                RAN_REPEATMASKER=1
            elif [[ -n "\${RAW_CONSENSUS}" || -n "\${CONSENSUS}" ]]; then
                cat "../../logs/\${BIN}.repeatmodeler.log" >&2
                printf '%s\n' \
                    'RepeatModeler produced repeat models without a nonempty classified consensus.' >&2
                exit 2
            elif grep -Fxq 'No families identified.  Perhaps the database is too small' "\${REPEATMODELER_LOG}" \
                && { grep -Fxq 'RepeatScout/RECON discovery complete: 0 families found' "\${REPEATMODELER_LOG}" \
                    || { grep -Fxq 'RepeatScout/RECON discovery complete:  families found' "\${REPEATMODELER_LOG}" \
                        && grep -Eq '^Round Time: .* Elapsed Time : 0 families discovered[.]\$' "\${REPEATMODELER_LOG}" \
                        && awk '/^ -- Input Database Coverage: [0-9]+ bp out of [0-9]+ bp [(] 100[.]00 % [)]\$/ \
                            && \$5 > 0 && \$5 == \$9 {covered=1} END {exit !covered}' "\${REPEATMODELER_LOG}"; }; }; then
                cp -L "\${FASTA_ABS}" "../../masked/\${BIN}.softmasked.fna"
                printf '%s\n' \
                    'RepeatModeler reported zero discovered families; RepeatMasker was not run.' \
                    > "../../logs/\${BIN}.repeatmasker.log"
                MASK_SOURCE=zero_repeat_families
                REPEAT_LIBRARY=none
                REPEAT_FAMILIES=0
                MASKER_JOBS=0
            else
                printf '%s\n' \
                    'RepeatModeler produced neither classified models nor the explicit zero-family marker.' >&2
                exit 2
            fi
            if [[ \${RECON_ONLY} -eq 1 ]]; then
                MASK_SOURCE="\${MASK_SOURCE}_recon_only"
            fi
            popd >/dev/null
            MASKED_OUT="09_repeat_masking/masked/\${BIN}.softmasked.fna"
            [[ \$(canonical_sha256 "\${FASTA}") == \$(canonical_sha256 "\${MASKED_OUT}") ]]
            MASKED_BASES=\$(lowercase_base_count "\${MASKED_OUT}")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "\${BIN}" "\${MASK_SOURCE}" \
                "\${REPEAT_LIBRARY}" "\${REPEAT_FAMILIES}" "\${MASKED_BASES}" \
                "\${MASKER_JOBS}" "${params.dfam_db_version}" \
                >> 09_repeat_masking/logs/masking.tsv
        done
        [[ \$(find 09_repeat_masking/masked -type f -name '*.softmasked.fna' -size +0c | wc -l) -eq \${#BINS[@]} ]]
    fi

    TOOL=protist-meta
    TOOL_VERSION="${params.workflow_version}"
    COMMAND='reuse supplied soft-masked input or skip when no eukaryotic bins are routed'
    DATABASE_NAME=null
    DATABASE_VERSION=null
    if [[ \${RAN_REPEATMODELER} -eq 1 ]]; then
        REPEATMODELER_OUTPUT=\$(
            pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
                --environment repeatmasker RepeatModeler -version
        )
        REPEATMODELER_VERSION=\$(printf '%s\n' "\${REPEATMODELER_OUTPUT}" \
            | awk '/^RepeatModeler version / {print \$3; exit}')
        [[ -n "\${REPEATMODELER_VERSION}" ]]
        TOOL=RepeatModeler
        TOOL_VERSION="RepeatModeler \${REPEATMODELER_VERSION}"
        COMMAND='RepeatModeler -threads ${task.cpus}; require classified models or the explicit native zero-family marker'
        DATABASE_NAME='"Dfam"'
        DATABASE_VERSION='"${params.dfam_db_version}"'
    fi
    if [[ \${RAN_REPEATMASKER} -eq 1 ]]; then
        REPEATMASKER_OUTPUT=\$(
            pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
                --environment repeatmasker RepeatMasker -v
        )
        REPEATMASKER_VERSION=\$(printf '%s\n' "\${REPEATMASKER_OUTPUT}" \
            | awk '/^RepeatMasker version / {print \$3; exit}')
        [[ -n "\${REPEATMASKER_VERSION}" ]]
        TOOL='RepeatModeler;RepeatMasker'
        TOOL_VERSION="RepeatModeler \${REPEATMODELER_VERSION}; RepeatMasker \${REPEATMASKER_VERSION}"
        COMMAND='RepeatModeler -threads ${task.cpus}; RepeatMasker -pa ${repeatMaskerWorkers} -xsmall; verify exact sequence identity ignoring mask case'
    fi
    if [[ \${RAN_RECON_ONLY} -eq 1 ]]; then
        COMMAND="\${COMMAND}; confirmed RepeatScout recovery: RepeatModeler -skipRS -threads ${task.cpus} (RECON-only discovery)"
    fi
    cat > 09_repeat_masking/stage.json <<JSON
    {
      "stage": "repeat_masking",
      "status": "\${STATUS}",
      "tool": "\${TOOL}",
      "tool_version": "\${TOOL_VERSION}",
      "database_name": \${DATABASE_NAME},
      "database_version": \${DATABASE_VERSION},
      "command": "\${COMMAND}",
      "reason": "\${REASON}",
      "outputs": {"masked_bins": "masked", "repeat_models": "models", "logs": "logs"}
    }
JSON
    """

    stub:
    """
    mkdir -p 09_repeat_masking/{masked,models,logs}
    printf '%s\n' '{"stage":"repeat_masking","status":"skipped","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"masked_bins":"masked","repeat_models":"models","logs":"logs"},"reason":"no_routed_eukaryotic_bins"}' > 09_repeat_masking/stage.json
    """
}

process MAP_RNA_EVIDENCE {
    tag "${meta.sample_id}"
    label 'gene'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly_dir), path(masking_dir), path(rna_reads)

    output:
    tuple val(meta), path('10_rna_alignment'), emit: result

    script:
    def rnaPreset = meta.rna_platform == 'illumina' ? 'splice:sr' \
        : meta.rna_platform == 'pacbio_isoseq' ? 'splice:hq' \
        : meta.rna_platform == 'ont_direct' ? 'splice -uf -k14' : 'splice'
    def rnaMapThreads = Math.max(task.cpus - 4, 1)
    def rnaSource = rna_reads ? rna_reads.toString() : ''
    """
    mkdir -p 10_rna_alignment/{bins,stats}
    BIN_COUNT=\$(find "${masking_dir}/masked" -maxdepth 1 -type f -name '*.softmasked.fna' -size +0c | wc -l)
    STATUS=skipped
    REASON=no_routed_eukaryotic_bins
    if [[ \${BIN_COUNT} -gt 0 && -z "${rnaSource}" ]]; then
        STATUS=pending
        REASON=no_sample_bound_rna_reads
        printf '%s\n' 'No compatible sample-bound RNA reads were supplied.' > 10_rna_alignment/pending.txt
    elif [[ \${BIN_COUNT} -gt 0 ]]; then
        STATUS=completed
        REASON=
        pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            minimap2 -ax ${rnaPreset} --secondary=no -t ${rnaMapThreads} \
            "${assembly_dir}/assembly_primary.fa" "${rnaSource}" \
          | pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            samtools view -b - \
          | pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            samtools sort -@ 1 -o all_rna.sorted.bam -
        pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            samtools flagstat -@ 4 all_rna.sorted.bam > 10_rna_alignment/stats/all.flagstat.txt
        pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            samtools view -@ 4 -b -F 2308 -q ${params.rna_min_mapq} all_rna.sorted.bam \
          | pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            samtools sort -@ 4 -o 10_rna_alignment/sample.splice.primary.unique.bam -
        pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            samtools index -@ 4 10_rna_alignment/sample.splice.primary.unique.bam
        pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
            samtools flagstat -@ 4 10_rna_alignment/sample.splice.primary.unique.bam \
            > 10_rna_alignment/stats/retained.flagstat.txt
        rm all_rna.sorted.bam

        mkdir -p tmp/faidx
        for FASTA in "${masking_dir}"/masked/*.softmasked.fna; do
            BIN=\$(basename "\${FASTA}" .softmasked.fna)
            FAI="tmp/faidx/\${BIN}.fai"
            pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
                samtools faidx --fai-idx "\${FAI}" "\${FASTA}"
            awk 'BEGIN {OFS="\t"} {print \$1, 0, \$2}' "\${FAI}" > bin.bed
            pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
                samtools view -@ 4 -b -L bin.bed \
                10_rna_alignment/sample.splice.primary.unique.bam \
                -o "10_rna_alignment/bins/\${BIN}.bam"
            pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
                samtools index -@ 4 "10_rna_alignment/bins/\${BIN}.bam"
        done
    fi

    MINIMAP2_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.core_manifest}" \
        --environment core minimap2 --version)
    MINIMAP2_VERSION=\${MINIMAP2_VERSION%%\$'\n'*}
    SAMTOOLS_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.core_manifest}" \
        --environment core samtools --version)
    SAMTOOLS_VERSION=\${SAMTOOLS_VERSION%%\$'\n'*}
    [[ -n "\${MINIMAP2_VERSION}" && -n "\${SAMTOOLS_VERSION}" ]]
    cat > 10_rna_alignment/stage.json <<JSON
    {
      "stage": "rna_alignment",
      "status": "\${STATUS}",
      "tool": "minimap2;samtools",
      "tool_version": "minimap2 \${MINIMAP2_VERSION}; \${SAMTOOLS_VERSION}",
      "database_name": null,
      "database_version": null,
      "command": "whole-assembly splice alignment; retain primary MAPQ>=${params.rna_min_mapq}; subset by bin contigs",
      "reason": "\${REASON}",
      "outputs": {"whole_assembly_bam": "sample.splice.primary.unique.bam", "bin_bams": "bins", "stats": "stats"}
    }
JSON
    """

    stub:
    """
    mkdir -p 10_rna_alignment/{bins,stats}
    printf '%s\n' '{"stage":"rna_alignment","status":"skipped","tool":"stub","tool_version":"stub","database_version":null,"command":"stub","outputs":{"bin_bams":"bins","stats":"stats"},"reason":"no_routed_eukaryotic_bins"}' > 10_rna_alignment/stage.json
    """
}

process BRAKER3_EUKARYOTES {
    tag "${meta.sample_id}"
    label 'gene'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(masking_dir), path(rna_dir), path(protein_reference)

    output:
    tuple val(meta), path('11_braker3'), emit: result

    script:
    def proteinSource = protein_reference ? protein_reference.toString() : ''
    def rnaPlatform = meta.rna_platform?.toString()?.toLowerCase() ?: ''
    """
    set -o pipefail
    mkdir -p 11_braker3/{annotations,proteins,augustus_config} braker_work/tmp
    export TMPDIR="\$PWD/braker_work/tmp"
    shopt -s nullglob
    BINS=("${masking_dir}"/masked/*.softmasked.fna)
    STATUS=skipped
    REASON=no_routed_eukaryotic_bins
    if [[ \${#BINS[@]} -gt 0 ]]; then
        HAS_PROTEINS=0
        if [[ -n "${proteinSource}" ]]; then
            [[ -s "${proteinSource}" ]]
            HAS_PROTEINS=1
        fi
        STATUS=completed
        REASON=
        MISSING_EVIDENCE=0
        AUG_SOURCE=\$(pixi run --as-is --manifest-path "${params.braker_manifest}" \
            bash -c 'for d in "\${AUGUSTUS_CONFIG_PATH:-}" "\${CONDA_PREFIX}/config" "\${CONDA_PREFIX}/share/augustus/config"; do if [[ -d "\${d}/species" ]]; then printf "%s\n" "\${d}"; exit 0; fi; done; exit 1')
        AUGUSTUS_CONFIG="\$PWD/11_braker3/augustus_config"
        for ENTRY in "\${AUG_SOURCE}"/*; do
            [[ "\${ENTRY##*/}" == species ]] && continue
            cp -a "\${ENTRY}" "\${AUGUSTUS_CONFIG}/"
        done
        mkdir -p "\${AUGUSTUS_CONFIG}/species"
        cp -a "\${AUG_SOURCE}/species/generic" "\${AUGUSTUS_CONFIG}/species/"

        for FASTA in "\${BINS[@]}"; do
            BIN=\$(basename "\${FASTA}" .softmasked.fna)
            BAM="${rna_dir}/bins/\${BIN}.bam"
            HAS_RNA=0
            if [[ -e "\${BAM}" || -e "\${BAM}.bai" ]]; then
                [[ -s "\${BAM}" && -s "\${BAM}.bai" ]]
                pixi run --as-is --manifest-path "${params.core_manifest}" --environment core \
                    samtools quickcheck -v "\${BAM}"
                ALIGNMENTS=\$(pixi run --as-is --manifest-path "${params.core_manifest}" \
                    --environment core samtools view -c "\${BAM}")
                if [[ \${ALIGNMENTS} -gt 0 ]]; then
                    HAS_RNA=1
                fi
            fi

            OUT="\$PWD/braker_work/\${BIN}"
            SPECIES="${meta.sample_id}_\${BIN}"
            ARGS=(braker.pl "--genome=\${FASTA}" "--species=\${SPECIES}" \
                "--AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG}" \
                "--workingdir=\${OUT}" "--threads=${task.cpus}" --gff3)
            if [[ \${HAS_PROTEINS} -eq 1 && \${HAS_RNA} -eq 1 ]]; then
                case "${rnaPlatform}" in
                    pacbio_isoseq|ont_cdna|ont_direct)
                        ARGS+=("--prot_seq=${proteinSource}" "--bam=\${BAM}")
                        ;;
                    illumina)
                        printf '%s\t%s\n' "\${BIN}" \
                            short_read_rna_incompatible_with_longrna_etp \
                            >> 11_braker3/pending.tsv
                        MISSING_EVIDENCE=1
                        continue
                        ;;
                    *)
                        printf '%s\t%s\n' "\${BIN}" unsupported_rna_platform \
                            >> 11_braker3/pending.tsv
                        MISSING_EVIDENCE=1
                        continue
                        ;;
                esac
            elif [[ \${HAS_PROTEINS} -eq 1 ]]; then
                ARGS+=("--prot_seq=${proteinSource}")
            elif [[ \${HAS_RNA} -eq 1 ]]; then
                case "${rnaPlatform}" in
                    illumina|pacbio_isoseq|ont_cdna|ont_direct)
                        ARGS+=("--bam=\${BAM}")
                        ;;
                    *)
                        printf '%s\t%s\n' "\${BIN}" unsupported_rna_platform \
                            >> 11_braker3/pending.tsv
                        MISSING_EVIDENCE=1
                        continue
                        ;;
                esac
            else
                printf '%s\t%s\n' "\${BIN}" no_compatible_rna_or_protein_evidence \
                    >> 11_braker3/pending.tsv
                MISSING_EVIDENCE=1
                continue
            fi
            mkdir -p "\${OUT}"
            pixi run --as-is --manifest-path "${params.braker_manifest}" \
                "\${ARGS[@]}"
            PUBLISHED="\$PWD/11_braker3/annotations/\${BIN}"
            mkdir -p "\${PUBLISHED}/species"
            for REQUIRED in braker.gtf braker.gff3 braker.aa braker.codingseq \
                hintsfile.gff braker.log what-to-cite.txt; do
                [[ -s "\${OUT}/\${REQUIRED}" ]]
                cp -L "\${OUT}/\${REQUIRED}" "\${PUBLISHED}/"
            done
            for HEADER_MAP in genome_header.map bam_header.map; do
                if [[ -f "\${OUT}/\${HEADER_MAP}" ]]; then
                    cp -L "\${OUT}/\${HEADER_MAP}" "\${PUBLISHED}/"
                fi
            done
            [[ -s "\${OUT}/species/\${SPECIES}/\${SPECIES}_parameters.cfg" ]]
            cp -aL "\${OUT}/species/\${SPECIES}" "\${PUBLISHED}/species/"
            cp -L "\${OUT}/braker.aa" "11_braker3/proteins/\${BIN}.faa"
        done
        if [[ \${MISSING_EVIDENCE} -eq 1 ]]; then
            STATUS=pending
            REASON=one_or_more_bins_lack_compatible_rna_or_protein_evidence
        fi
    fi

    BRAKER_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.braker_manifest}" \
        braker.pl --version)
    BRAKER_VERSION=\${BRAKER_VERSION%%\$'\n'*}
    if GENEMARK_HELP=\$(pixi run --as-is --quiet --manifest-path "${params.braker_manifest}" \
        gmetp.pl); then
        GENEMARK_STATUS=0
    else
        GENEMARK_STATUS=\$?
    fi
    if [[ \${GENEMARK_STATUS} -ne 1 ]]; then
        printf 'GeneMark-ETP usage probe exited %s, expected 1\n' \
            "\${GENEMARK_STATUS}" >&2
        exit 1
    fi
    if ! GENEMARK_VERSION=\$(printf '%s\n' "\${GENEMARK_HELP}" \
        | grep -E '^ETP version [0-9]+[.][0-9]+\$'); then
        GENEMARK_VERSION=
    fi
    if [[ "\${GENEMARK_VERSION}" != 'ETP version 1.02' ]]; then
        printf '%s\n' 'GeneMark-ETP usage probe did not report ETP version 1.02' >&2
        exit 1
    fi
    PROTHINT_VERSION=\$(pixi run --as-is --quiet --manifest-path "${params.braker_manifest}" \
        prothint.py --version)
    PROTHINT_VERSION=\${PROTHINT_VERSION%%\$'\n'*}
    [[ -n "\${BRAKER_VERSION}" && -n "\${PROTHINT_VERSION}" ]]
    cat > 11_braker3/stage.json <<JSON
    {
      "stage": "braker3",
      "status": "\${STATUS}",
      "tool": "BRAKER3",
      "tool_version": "\${BRAKER_VERSION}@35280d0c2f56dde026ab346c0255740c11ab80a1; \${GENEMARK_VERSION}@64a69cbc11e15d51779da021d794cbe9c6d695a7; \${PROTHINT_VERSION}",
      "database_name": "protein_reference",
      "database_version": "protein lineage ${meta.protein_lineage ?: 'not supplied'}",
      "command": "seed private AUGUSTUS config with species/generic; braker.pl with explicit private AUGUSTUS_CONFIG_PATH, default softmask handling and evidence-inferred mode; ET for RNA only, EP for protein only, patched long-read ETP for long RNA plus protein; retain native work outside publication and copy final outputs, logs, header maps and trained species",
      "reason": "\${REASON}",
      "outputs": {"annotations": "annotations", "proteins": "proteins"}
    }
JSON
    """

    stub:
    """
    mkdir -p 11_braker3/{annotations,proteins}
    printf '%s\n' '{"stage":"braker3","status":"skipped","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"annotations":"annotations","proteins":"proteins"},"reason":"no_routed_eukaryotic_bins"}' > 11_braker3/stage.json
    """
}

process PRODIGAL_GV_GENES {
    tag "${meta.sample_id}"
    label 'gene'
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(routing_dir), path(viral_dir)

    output:
    tuple val(meta), path('12_prodigal_gv'), emit: result

    script:
    """
    mkdir -p 12_prodigal_gv/{proteins,genes,gff}
    "${params.cli}" prepare-gene-inputs \
        --routing-dir "${routing_dir}" \
        --viral-dir "${viral_dir}" \
        --output-dir 12_prodigal_gv/inputs
    [[ -s 12_prodigal_gv/inputs/sources.tsv ]]
    [[ -s 12_prodigal_gv/inputs/origins.tsv ]]

    SOURCE_COUNT=\$(awk -F '\t' \
        'NR > 1 && \$6 == "included" {count++} END {print count + 0}' \
        12_prodigal_gv/inputs/sources.tsv)
    STATUS=skipped
    REASON=no_exclusive_prodigal_sources
    if [[ \${SOURCE_COUNT} -gt 0 ]]; then
        STATUS=completed
        REASON=
        while IFS=\$'\t' read -r SOURCE CATEGORY RELATIVE_FASTA; do
            FASTA="12_prodigal_gv/inputs/\${RELATIVE_FASTA}"
            [[ -s "\${FASTA}" ]]
            MODE_ARGS=(-p meta)
            if [[ "\${CATEGORY}" == prokaryotic_bin \
                && -n "${meta.genetic_code}" && "${meta.genetic_code}" != auto ]]; then
                MODE_ARGS=(-p single -g "${meta.genetic_code}")
            fi
            pixi run --as-is --manifest-path "${params.gene_manifest}" --environment prodigal-gv \
                prodigal-gv -i "\${FASTA}" -a "12_prodigal_gv/proteins/\${SOURCE}.faa" \
                -d "12_prodigal_gv/genes/\${SOURCE}.fna" -f gff \
                -o "12_prodigal_gv/gff/\${SOURCE}.gff" "\${MODE_ARGS[@]}"
            [[ -s "12_prodigal_gv/gff/\${SOURCE}.gff" ]]
            grep -q '^# Model Data:.*transl_table=' "12_prodigal_gv/gff/\${SOURCE}.gff"
            [[ -e "12_prodigal_gv/proteins/\${SOURCE}.faa" ]]
            [[ -e "12_prodigal_gv/genes/\${SOURCE}.fna" ]]
            CDS_COUNT=\$(awk '!/^#/ && NF {count++} END {print count + 0}' \
                "12_prodigal_gv/gff/\${SOURCE}.gff")
            if [[ \${CDS_COUNT} -gt 0 ]]; then
                [[ -s "12_prodigal_gv/proteins/\${SOURCE}.faa" ]]
                [[ -s "12_prodigal_gv/genes/\${SOURCE}.fna" ]]
            else
                [[ ! -s "12_prodigal_gv/proteins/\${SOURCE}.faa" ]]
                [[ ! -s "12_prodigal_gv/genes/\${SOURCE}.fna" ]]
            fi
        done < <(awk -F '\t' \
            'NR > 1 && \$6 == "included" {print \$1 "\t" \$2 "\t" \$4}' \
            12_prodigal_gv/inputs/sources.tsv)
    fi
    if ! VERSION_OUTPUT=\$(pixi run --as-is --quiet --manifest-path "${params.gene_manifest}" \
        --environment prodigal-gv prodigal-gv -v 2>&1); then
        printf '%s\n' "\${VERSION_OUTPUT}" >&2
        exit 1
    fi
    VERSION=\$(printf '%s\n' "\${VERSION_OUTPUT}" \
        | awk '/^Prodigal V/ {print; exit}' | sed 's/"/\\\\"/g')
    [[ -n "\${VERSION}" ]]
    cat > 12_prodigal_gv/stage.json <<JSON
    {
      "stage": "prodigal_gv",
      "status": "\${STATUS}",
      "tool": "Prodigal-gv",
      "tool_version": "\${VERSION}",
      "database_name": "prodigal_gv_models",
      "database_version": "bundled viral models",
      "command": "prepare exclusive routed-bin and unassigned geNomad sources; prodigal-gv meta mode with per-CDS viral codes; single mode for explicit prokaryotic genetic-code override",
      "reason": "\${REASON}",
      "outputs": {"sources": "inputs/sources.tsv", "origins": "inputs/origins.tsv", "proteins": "proteins", "genes": "genes", "gff": "gff"}
    }
JSON
    """

    stub:
    """
    mkdir -p 12_prodigal_gv/{inputs,proteins,genes,gff}
    printf 'source_id\tcategory\tbin_id\tfasta\torigin_map\tstatus\treason\n' > 12_prodigal_gv/inputs/sources.tsv
    printf 'contig_id\tparent_contig\tstart\tend\tsource_id\tbin_id\tcategory\tstatus\treason\n' > 12_prodigal_gv/inputs/origins.tsv
    printf '%s\n' '{"stage":"prodigal_gv","status":"skipped","tool":"stub","tool_version":"stub","database_version":"stub","command":"stub","outputs":{"sources":"inputs/sources.tsv","origins":"inputs/origins.tsv","proteins":"proteins","genes":"genes","gff":"gff"},"reason":"no_exclusive_prodigal_sources"}' > 12_prodigal_gv/stage.json
    """
}
