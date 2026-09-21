# Input and dependency contract

The sample sheet is tab-separated. `sample_id`, `platform`, and `reads` are
required. Identifiers use letters, digits, periods, underscores, and hyphens.
Paths are resolved relative to the sample sheet; they must name nonempty files.

| Field | Meaning |
|---|---|
| `platform` | `ONT` or `PACBIO_HIFI`; PacBio CLR is unsupported. |
| `reads` | Basecalled DNA FASTQ, optionally gzip-compressed. |
| `assembly` | Optional existing assembly for a post-assembly run. |
| `rna_reads`, `rna_platform` | Optional RNA evidence; both fields are required together. |
| `protein_reference`, `protein_lineage` | Optional protein FASTA and explicit lineage label. |
| `genetic_code` | Prokaryotic gene-calling override, or `auto`; does not assign a code to every organism in a mixed sample. |
| `softmasked` | Optional masked assembly corresponding to `assembly`. |
| `basecaller`, `basecaller_model` | Optional ONT provenance. Either value may be unknown; leave it blank. Leave both empty for HiFi. |
| `library_prep` | Optional library preparation provenance. |

RNA platforms are `ONT_CDNA`, `ONT_DIRECT`, `PACBIO_ISOSEQ`, or `ILLUMINA`.
Caller support is checked separately. DNA and RNA input paths must differ.

ONT assembly uses Myloasm's `--nano-r10` mode for R10 SUP/HAC reads. This mode
assumes roughly 97% or higher median read accuracy. The `ONT` label alone does
not establish that accuracy. Record the caller and exact model when each is
present in the sequencing metadata. A model does not establish the caller's
software version. PacBio HiFi uses `--hifi`.

An explicit `genetic_code` selects Prodigal-GV single-genome mode, which accepts
codes 1–6, 9–16 and 21–25. `auto` uses metagenomic mode and retains each
prediction's reported code.

`conf/databases.yaml` is the path registry. Paths are relative to that file,
including local Pixi manifests and externally maintained reference resources.
The registry records versions and expected files. Prepare reference files before
launching an analysis. GVClass retains its native reference setup and update
behavior. The workflow accepts GVClass output only when native metadata matches
the registered database path and version. Classification uses local reference
files.

Pixi manages all application software. Nextflow and Java belong to
the main environment. Incompatible tools use separate manifests under `envs/`.
Some tools reuse existing Pixi workspaces named in the registry. Runtime calls
use `pixi run --as-is`, which neither installs packages nor updates locks.
Preflight rejects missing executables and paths outside the activated Pixi
environment.

Resume verifies the saved sample sheet, registry, prepared configuration and
source identity. The identity includes every referenced Pixi manifest and lock.
That check does not detect external application code changes when the manifest
and lock stay unchanged.
