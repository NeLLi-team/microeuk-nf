# Input and dependency contract

## Sample sheet

The sample sheet is tab-separated. Copy `data/samples.example.tsv` to a file of
your choice and edit its single row. The template contains all 13 supported
fields; keep the header unchanged. The allocation launcher currently requires
exactly one data row.

`sample_id`, `platform`, and `reads` are required. A sample identifier must
start with a letter or digit; its remaining characters may also include periods,
underscores, and hyphens. File paths are resolved relative to the copied sample
sheet and must name nonempty files.

| Field | Meaning |
|---|---|
| `sample_id` | Unique identifier of 1 to 80 characters; the first must be alphanumeric. |
| `platform` | `ONT` or `PACBIO_HIFI`; PacBio CLR is unsupported. |
| `reads` | Basecalled DNA FASTQ, optionally gzip-compressed. |
| `assembly` | Optional existing assembly for a post-assembly run. |
| `rna_reads`, `rna_platform` | Optional RNA evidence; both fields are required together. |
| `protein_reference`, `protein_lineage` | Optional protein FASTA and explicit lineage label; both fields are required together. |
| `genetic_code` | Prokaryotic gene-calling override, or `auto`; it does not assign a code to every organism in a mixed sample. |
| `softmasked` | Optional masked assembly corresponding to `assembly`; supplying it requires `assembly`. |
| `basecaller`, `basecaller_model` | Optional ONT provenance. Either value may be unknown and left blank. Leave both empty for HiFi. |
| `library_prep` | Optional library preparation provenance. |

RNA platforms are `ONT_CDNA`, `ONT_DIRECT`, `PACBIO_ISOSEQ`, or `ILLUMINA`.
Caller support is checked separately. DNA and RNA input paths must differ.

ONT assembly uses Myloasm's `--nano-r10` mode for R10 SUP/HAC reads. This mode
assumes roughly 97% or higher median read accuracy. The `ONT` label alone does
not establish that accuracy. Record the caller and exact model when sequencing
metadata provides them. A model name does not establish the caller software
version. PacBio HiFi uses `--hifi`.

Read QC makes one pass over the raw file. SeqKit records the raw statistics
from the same stream that the filter reads, and records the filtered statistics
from the filter's output, both with `seqkit stats -a`. The `-a` option adds
further columns, such as N50 and the Q20 and Q30 fractions, after the standard
eight. By default, Chopper retains ONT reads at least 500 bp long with a mean Q
score of at least 10. HiFi reads are filtered with SeqKit at a minimum length of
1,000 bp. De novo assembly and DNA mapping use the filtered reads.

An explicit `genetic_code` selects Prodigal-GV single-genome mode, which accepts
codes 1–6, 9–16, and 21–25. `auto` uses metagenomic mode and retains each
prediction's reported code.

## Optional analysis stages

Full mode runs without gene calling by default. Gene calling and functional
annotation are opt-in. Repeat modeling and masking are a further opt-in that
needs gene calling. The allocation launcher and `protist-meta prepare` accept
these switches:

| Switch | Effect in full mode |
|---|---|
| `--enable-gene-calling` | Enables dedicated Prodigal-GV gene calling, RNA alignment, BRAKER3, and functional annotation. |
| `--skip-gene-calling` | Selects the default. Omits those stages and repeat masking. |
| `--skip-annotation` | With gene calling on, omits eggNOG-mapper and InterProScan functional annotation. |
| `--enable-repeat-masking` | With gene calling on, enables RepeatModeler/RepeatMasker, or uses the supplied softmasked assembly. |
| `--skip-repeat-masking` | Explicitly selects the default unmasked route. |

The defaults are `skip_gene_calling = true`, `skip_annotation = false`, and
`skip_repeat_masking = true`. `--enable-gene-calling` and `--skip-gene-calling`
are mutually exclusive, and so are the two repeat-masking switches. Functional
annotation needs gene calling. Quality assessment, contig classification
(Tiara), taxonomic classification, SSU analysis, viral screening, eukaryote QC
(BUSCO and tidk), and prokaryotic characterization always run in full mode. Gene
prediction used internally by those tools is unaffected. Core mode excludes all
optional stages regardless of these switches. When repeat masking is off, RNA
alignment and BRAKER3 use routed eukaryotic bin FASTAs directly. BRAKER3 receives
`--softmasking_off`. When full mode includes gene calling, preparation rejects
supplied softmasked input unless repeat masking is enabled.

Direct Nextflow commands use `--skip_repeat_masking`, `--skip_gene_calling` and
`--skip_annotation`, with underscores. These parameters must be booleans. A
command-line value such as `--skip_gene_calling false` arrives as a string and is
rejected; set the value in a params file or a configuration file instead, for
example `{"skip_gene_calling": false}` with `-params-file`. The hyphenated
launcher switches do not set these Nextflow parameters. Preparation records the
requested values in `inputs.json` and writes the matching Nextflow
configuration. Nextflow rejects an override that conflicts with a recorded
preparation choice. Analyses excluded by these choices are recorded as skipped,
with no products or counts; gene calling that was not requested carries the
reason `Gene calling not requested (skip_gene_calling=true).`

All supplied sample paths and evidence-pair fields are validated even when
their downstream analysis is disabled. Leave unused optional fields empty.

## Dependency registry

`conf/databases.yaml` is the path and version registry. Relative paths are
resolved from `conf/`; this includes repository manifests, sibling workspaces,
and local reference resources. The registry validates expected files and
executables but does not download or install them.

The repository supplies Pixi manifests for Nextflow, Java, Python, core tools,
BBTools, the quality and classification environments, RepeatMasker, Prodigal-GV,
eggNOG-mapper, the local BRAKER3 recipe, and the local InterProScan recipe. Full
mode also requires these targets:

| Registry entry | Required target |
|---|---|
| `checkm_manifest` | `envs/qc`: CheckM1, CheckM2, GTDB-Tk, and Symclatron Pixi environments |
| `viral_manifest` | `envs/viral`: geNomad Pixi environment (`genomad`) |
| `checkv_manifest` | `envs/viral`: CheckV Pixi environment (`checkv`) |
| `contig_manifest` | `envs/contig`: Tiara Pixi environment (`tiara`) |
| `eukqc_manifest` | `envs/eukqc`: BUSCO (`busco`) and tidk (`tidk`) Pixi environments |
| `busco_db` | BUSCO lineage directory `eukaryota_odb12.2`, with sentinel `dataset.cfg` |
| `checkeuk_app` | CheckEUK application and Pixi environment, a frozen clone |
| `gvclass_app` | GVClass application and Pixi environment, a frozen clone |
| `ssuextract_app` | SSUextract application, Pixi environment, workflow files, and database manager, a frozen clone |

The three `*_app` entries point to clones below `../protist-meta-runtimes/` and
carry a `git_commit`. Preparation and the launcher check each clone with
`git rev-parse HEAD` and `git status --porcelain`: the commit must match and the
tree must be clean. Preparation writes `params.<name>_commit` for each app, and
the stage manifests record the first 12 characters of the commit in the tool
version, for example `CheckEUK <version> (<commit>)`.

All targets and their registered reference resources are required for full mode.
Core mode requires the root and BBTools environments plus the registered
QuickClade reference. Classification uses local reference files; remote fallback
is disabled.

Gene calling off, the default, does not require `gene_manifest`,
`braker_manifest`, or `dfam_db`. Gene calling on with repeat masking off requires
`gene_manifest` and `braker_manifest` for Prodigal-GV and BRAKER3, but not
`dfam_db` or the RepeatModeler and RepeatMasker executables. Disabled functional
annotation does not require `annotation_manifest`, `interpro_manifest`,
`eggnog_db`, or `interproscan_db`. These entries are excluded from dependency
validation and preflight probes. All other enabled dependencies remain required.

The Dori setup scripts also use fixed local paths. `scripts/install-tools.slurm`
probes a fixed eggNOG data directory.
`scripts/install-custom.slurm` uses fixed InterProScan data and smoke-input
paths. The InterProScan package build requires the local source distribution,
the checked source manifests, and `.nellidb_version/version.txt`. Editing only
the registry does not redirect these setup paths.

Runtime calls use `pixi run --as-is` with installation and lock updates
disabled. Preflight rejects missing executables, missing sentinels, and
executables outside the activated Pixi environment. It also compares the
`DB_VERSION` file of the CheckEUK and GVClass databases with the registered
versions and checks that `myloasm --help` lists `--hifi` and `--clean-dir`. See
[Run and resume on
Dori](../how-to/run-on-dori.md) for the ordered setup.

## Resume identity

Resume verifies the saved sample sheet, registry, prepared configuration, and
source identity. Source identity includes every referenced Pixi manifest and
lock, and the Git commit of each `*_app` clone. Other changes to application code
are not detected when the manifest, lock, and commit remain unchanged. The
effective ONT minimum read length is always recorded and compared.

Resume requires the same optional-stage switches used for preparation, except
that the saved gene-calling choice is adopted when the command line does not
state one. A preparation without a field means that switch was false, so gene
calling, annotation, and repeat masking were all enabled. Changing the stage
selection requires a new run directory, including when using `resume-reviewed`
for a source repair.
