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

Every read-based run records raw FASTQ statistics with SeqKit before filtering.
By default, Chopper retains ONT reads at least 500 bp long with a mean Q score
of at least 10. SeqKit then records statistics for the filtered reads. De novo
assembly and DNA mapping use the filtered reads.

An explicit `genetic_code` selects Prodigal-GV single-genome mode, which accepts
codes 1–6, 9–16, and 21–25. `auto` uses metagenomic mode and retains each
prediction's reported code.

## Dependency registry

`conf/databases.yaml` is the path and version registry. Relative paths are
resolved from `conf/`; this includes repository manifests, sibling workspaces,
and local reference resources. The registry validates expected files and
executables but does not download or install them.

The repository supplies Pixi manifests for Nextflow, Java, Python, core tools,
BBTools, RepeatMasker, Prodigal-GV, eggNOG-mapper, the local BRAKER3 recipe, and
the local InterProScan recipe. Full mode also requires these separately managed
targets:

| Registry entry | Required workspace |
|---|---|
| `checkm_manifest` | CheckM1, CheckM2, GTDB-Tk, and Symclatron Pixi environments |
| `viral_manifest` | geNomad Pixi environment |
| `checkv_manifest` | CheckV Pixi environment |
| `checkeuk_app` | CheckEUK application and Pixi environment |
| `gvclass_app` | GVClass application and Pixi environment |
| `ssuextract_app` | SSUextract application, Pixi environment, workflow files, and database manager |

All six targets and their registered reference resources are required for full
mode. Core mode requires the root and BBTools environments plus the registered
QuickClade reference. Classification uses local reference files; remote fallback
is disabled.

The Dori setup scripts also use fixed local paths. `scripts/install-tools.slurm`
probes a fixed eggNOG data directory.
`scripts/install-custom.slurm` uses fixed InterProScan data and smoke-input
paths. The InterProScan package build requires the local source distribution,
the checked source manifests, and `.nellidb_version/version.txt`. Editing only
the registry does not redirect these setup paths.

Runtime calls use `pixi run --as-is` with installation and lock updates
disabled. Preflight rejects missing executables, missing sentinels, and
executables outside the activated Pixi environment. See [Run and resume on
Dori](../how-to/run-on-dori.md) for the ordered setup.

## Resume identity

Resume verifies the saved sample sheet, registry, prepared configuration, and
source identity. Source identity includes every referenced Pixi manifest and
lock. External application code changes are not detected when its manifest and
lock remain unchanged.
