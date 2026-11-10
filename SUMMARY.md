# Source-scoped validation

The 2026-11-10 upgrade is newer than every validation record below. It moves the
pipeline environments into the repository, freezes CheckEUK, GVClass, and
SSUextract as commit-verified clones, splits the bin screens and prokaryote QC
into parallel processes, adds Tiara contig classification and BUSCO and tidk
eukaryote QC, reads the input once in read QC, and turns gene calling off by
default. The records after the next section describe the older source revisions
that they name, not the upgraded pipeline.

## Upgrade validation, 2026-11-10

| Case | Job | Commit | Evidence |
|---|---:|---|---|
| MICRO, full mode | 27725402 | `e83ae52` | 21 completed tasks in 4 min 47 s; catalog, bundle, notebook, and HTML written |
| ONT pilot `X0341_p01`, full mode | 27722802 | `0f217d0` | 21 completed tasks in 2 h 02 min; catalog and report audit below |
| MICRO, `resume-reviewed` of job 27725402 | 27752569 | `f01bacd` | 11 cached and 10 completed tasks in 1 min 12 s |

The pilot ran with `--ont-min-length 1000` and gene calling off. Commit `e83ae52`
differs from `0f217d0` only in three CPU labels and the start order of the viral
screen. Commit `f01bacd` differs from `e83ae52` only in the tidk threshold.

The reviewed resume kept every task up to the merged bin screen and ran
`ROUTE_BINS` and all later tasks again, including CheckM1, Symclatron, and
GTDB-Tk. A source repair that needs a reviewed resume therefore repeats the
routed-bin analyses. A plain `resume` with unchanged source was not tested for
the upgraded pipeline.

Local gates at `e83ae52`: 398 unit tests pass, the Python style gate passes,
`generate_models.py --check` passes, `nextflow lint` reports no error, and
`tests/check_stage_merge.sh` passes.

The pilot catalog holds 130 bins, 1,192 contigs, 408 quality rows, 931 taxonomy
rows, 59 SSU loci, 28 viruses, 59 phenotype rows, and 175 artifacts. SQLite
integrity is `ok` with zero foreign-key violations. The executed notebook has
8 code cells and 0 errors. Routing gave 59 prokaryotic bins, 1 eukaryotic bin,
1 viral bin, and 69 unresolved bins, 6 of them conflicting. The one eukaryotic
bin is 90 kb; BUSCO found none of 125 markers in it, so the pilot does not
test BUSCO on a real eukaryotic genome.

Stage times in the pilot: read QC 42 s, assembly 1 min 45 s, mapping and
binning 1 min, geNomad and CheckV 7 min, SSUextract 9 min, CheckM2 3 min,
CheckEUK 19 min, GVClass 28 min, CheckM1 4 min, Symclatron 5 min, GTDB-Tk
1 h 15 min. GTDB-Tk spends most of that time loading its six class-level
reference trees, so its time does not scale with the bin count.

Limits: one 1% ONT subset and one four-contig data set with a supplied
assembly. No HiFi input and no gene-calling branch were executed on the
upgraded pipeline.

microeuk-nf is implemented for ONT and PacBio HiFi inputs. Validation covers
read QC, Myloasm assembly, depth mapping, QuickBin, QuickClade, all-bin quality
and taxonomy evidence, SSU extraction, viral screening, routing, repeat masking,
RNA evidence, prokaryotic and eukaryotic gene calling, functional annotation,
LinkML/Pydantic validation, atomic SQLite publication, notebook execution, HTML
export, and cached resume.

| Case | Job | Source prefix | Accepted evidence |
|---|---:|---|---|
| Current ONT | 26174581 | `4f734784d5e5` | 17 completed tasks; 17 cached tasks on resume; catalog, report, and artifact audit passed |
| Positive prokaryotic | 26165905 | `4d733a063539` | Native workflow, catalog, report, artifact inspection, and cached resume passed |
| Positive RNA-supported eukaryotic | 26165903 | `4d733a063539` | Native workflow, BRAKER3 branch, catalog, report, artifact inspection, and cached resume passed |

Full source fingerprints:

```text
ONT:      sha256:4f734784d5e5623955354126ff3a8e972f6047c126d5b1caffa951246322581e
Positive: sha256:4d733a06353939e25593a34fc4ef0ce6ae260e6bdab34391cb5167de660f9aaf
```

The recorded `sha256:` value is a workflow source fingerprint over code,
configuration, schemas, scripts, Pixi and tool-recipe files, and referenced Pixi
manifests and locks. It is distinct from a Git commit and excludes Markdown.

The positive cases are evidence for their recorded `4d733a…` revision. They are
not current-source reruns. Later changes were bounded to Pixi shell handling,
launcher return behavior, RepeatMasker identifier restoration, and selection of
RepeatModeler's final zero-family consensus.

## Current ONT acceptance

Before scientific execution, job 26174581 passed 143 tests, style checks for 32
Python files, generated-model consistency, strict documentation, Nextflow
configuration and previews, the 38-entry Pixi guard, and all 10 applied-file
hash checks. The job completed with exit `0:0` in 11:26:56 using 32 allocated
CPUs and a 384 GB request. It ran sample `X0341_p01` as
`run_823b4b2ea3996922`.

The completion audit matched all 18 catalog record counts and all 16 report
tables. SQLite integrity was `ok`, foreign-key violations were zero, and all
1,448 artifact hashes passed. Catalog, notebook, and HTML hashes were unchanged
after resume.

The accepted catalog contains:

- 150 bins and 1,196 contigs;
- 616 quality rows and 1,300 taxonomy rows;
- 57 SSU loci and 29 viruses;
- 143 phenotype rows;
- 90,992 genes and 1,680,663 annotations; and
- 1,448 artifacts.

Native quality matching covered all 150 targets for each tool. CheckEUK matched
150 QC and 36 taxonomy rows. CheckM2 matched 150 QC rows. GVClass matched one QC
and 469 taxonomy rows.

The executed notebook has 8 code cells, 3 figures, and 0 errors. The
catalog is 1,437,749,248 bytes, the notebook is 3,093,437 bytes, and the HTML
export is 2,282,848 bytes. Visual review found no clipping, missing panels,
incorrect units, or caption contradictions. Figure 3 has a minor
label-readability limitation.

Batch accounting recorded `MaxRSS` 192,589,848 K. This is batch accounting, not
a per-tool peak. The earlier 1 TB Nextflow trace value can sum shared pages and
peaks from different times, so it does not establish deduplicated physical
memory use or a memory-limit violation.

## Positive-case coverage

The positive prokaryotic catalog contains 4,323 genes, 233,766 annotations,
seven SSU loci, five viruses, nine quality rows, 45 taxonomy rows, one phenotype
row, and 58 artifacts.

The positive RNA-supported eukaryotic catalog contains 8,782 protein-coding
transcript/protein models, 504,043 annotations, 21 SSU loci, six quality rows,
20 taxonomy rows, and 702 artifacts.

## Scientific limits

- All six ONT eukaryotic analysis bins lacked compatible RNA or protein
  evidence. Their RNA and BRAKER3 stage records remain pending scientific
  states. The ONT case therefore makes no eukaryotic gene-call claim. Job
  26165903 provides the positive RNA-supported BRAKER3 evidence.
- All 57 ONT SSU assignments remain `Unclassified` with
  `no_runtime_calibration`. Reference candidates are evidence, not calibrated
  final assignments.
- GVClass domain, confidence, and model-reliability fields are exploratory
  evidence rather than final taxonomy.
- Tool-specific absent estimates remain absent. The workflow does not replace
  missing values with zero or create status-only taxonomy and quality rows.
- Route labels select downstream analyses. They do not establish final taxonomy
  or MAG acceptance.

## Reproducibility boundary

Pixi manages all declared workflow dependencies, including Nextflow and Java.
The host supplies `/bin/bash` for the initial bootstrap and Nextflow's outer
task launcher. Inner process shells and workflow tools use guarded Pixi
environments with `pixi run --as-is`; runtime does not install software or
update locks.

One bounded Slurm allocation contains local Nextflow execution. The 17
Nextflow tasks of the accepted ONT case are not 17 Slurm jobs. The central path and version registry is
`conf/databases.yaml`; the [input and dependency
contract](docs/reference/inputs.md) describes its portability limits.

Full completion receipts, reports, catalogs, and accounting records remain
local artifacts under ignored `tasks/` and `results/` paths. They are not part
of the GitHub clone. This summary preserves the accepted source revisions,
counts, and interpretation limits without representing local files as
repository links.
