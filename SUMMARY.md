# Validation status

The protist metagenome workflow is implemented and validated for ONT and PacBio
HiFi inputs. It covers read QC, Myloasm assembly, depth mapping, QuickBin,
QuickClade, all-bin quality and taxonomy evidence, SSU extraction, viral
screening, routing, repeat masking, RNA evidence, prokaryotic and eukaryotic
gene calling, functional annotation, LinkML/Pydantic validation, atomic SQLite
publication, an executed notebook, HTML export, and cached resume.

## Current production validation

The current source revision is
`sha256:4f734784d5e5623955354126ff3a8e972f6047c126d5b1caffa951246322581e`.
Before scientific execution, job 26174581 passed 143 tests, style checks for 32
Python files, generated-model consistency, strict documentation, Nextflow
configuration and previews, the 38-entry Pixi guard, and all 10 applied-file
hash checks.

Job 26174581 completed with exit 0:0 in 11:26:56 using 32 allocated CPUs and a
384 GB request. It ran sample `X0341_p01` as `run_823b4b2ea3996922`. The initial
trace has 17 COMPLETED, exit-zero tasks. The resume trace has 17 CACHED,
exit-zero tasks. Source and input identity remained unchanged.

The [completion audit](tasks/ont-pilot-retry-26174581/completion-audit.json)
has `status: passed`. It matched all 18 catalog record counts and all 16 report
tables. SQLite integrity is `ok`, foreign-key violations are zero, and all
1,448 artifact hashes pass. Catalog, notebook, and HTML hashes were unchanged
after resume.

The ONT catalog contains:

- 150 bins and 1,196 contigs;
- 616 quality rows and 1,300 taxonomy rows;
- 57 SSU loci and 29 viruses;
- 143 phenotype rows;
- 90,992 genes and 1,680,663 annotations; and
- 1,448 artifacts.

Native quality matching covered all 150 targets for each tool. CheckEUK matched
150 QC and 36 taxonomy rows. CheckM2 matched 150 QC rows. GVClass matched one QC
and 469 taxonomy rows.

The executed notebook has eight code cells, three figures, and no errors. It is
3,093,437 bytes. The HTML export is 2,282,848 bytes, and the SQLite catalog is
1,437,749,248 bytes. The [figure review](tasks/ont-pilot-retry-26174581/figure-review/findings.md)
passed with minor Figure 3 label-readability findings and no clipping, missing
panels, wrong units, or caption contradictions.

The [final accounting receipt](tasks/ont-pilot-retry-26174581/final-accounting-20260921.txt)
records batch `MaxRSS` 192,589,848 K. This is batch accounting, not a per-tool
peak. The earlier 1 TB Nextflow trace value can count shared pages and peaks
from different times more than once; it does not establish deduplicated
physical memory use or a memory-limit violation.

## Positive-case coverage

The positive prokaryotic case, job 26165905, completed its native workflow,
catalog, report, artifact inspection, and cached resume at source revision
`sha256:4d733a06353939e25593a34fc4ef0ce6ae260e6bdab34391cb5167de660f9aaf`.
Its reconciled catalog contains 4,323 genes, 233,766 annotations, seven SSU
loci, five viruses, nine quality rows, 45 taxonomy rows, one phenotype row, and
58 artifacts. See its [inspection](tasks/positive-prok-retry-26165905/inspection.json)
and [completion audit](tasks/ont-mask-ids-fix/implementation-proof-26171416/prok-completion-audit.json).

The positive RNA-supported eukaryotic case, job 26165903, completed its native
workflow, catalog, report, artifact inspection, and cached resume at the same
`4d733a…` source revision. Its reconciled catalog contains 8,782
protein-coding transcript/protein models, 504,043 annotations, 21 SSU loci, six
quality rows, 20 taxonomy rows, and 702 artifacts. See its
[inspection](tasks/positive-euk-retry-26165903/inspection.json) and
[completion audit](tasks/micro-full-26173111/euk-completion-audit.json).

These receipts remain valid for their recorded source revisions. They are not
represented as current-source reruns. The changes after `4d733a…` were reviewed
and bounded to Pixi shell handling, launcher return behavior, RepeatMasker
identifier restoration, and selection of RepeatModeler's final zero-family
consensus.

## Scientific limits

- All six ONT eukaryotic analysis bins lacked compatible RNA or protein
  evidence. Their RNA and BRAKER stage records remain pending by design, and
  the ONT case makes no eukaryotic gene-call claim. Job 26165903 supplies the
  positive RNA-supported BRAKER evidence.
- All 57 ONT SSU assignments remain `Unclassified` with
  `no_runtime_calibration`. Reference candidates are evidence, not calibrated
  final assignments.
- GVClass domain, confidence, and model-reliability fields remain exploratory
  evidence rather than final taxonomy.
- Tool-specific absent estimates remain absent. The workflow does not replace
  missing values with zero or create status-only taxonomy and quality rows.
- Route labels select downstream analyses; they are not final taxonomic
  assignments or MAG acceptance claims.

## Reproducibility and delivery

All declared application dependencies, including Nextflow and Java, are
Pixi-managed. The outer Nextflow launcher uses the documented operating-system
Bash bootstrap; workflow tools and process shells use guarded Pixi
environments. One bounded Slurm allocation contains local Nextflow execution;
the 17 Nextflow tasks are not 17 Slurm jobs.

The central database registry is `conf/databases.yaml`, with its contract in
`docs/reference/inputs.md`. The final collaborator artifacts are:

- `results/validation/ont-pilot-retry-26174581/results/catalog/protist-meta.sqlite`;
- `results/validation/ont-pilot-retry-26174581/results/report/report.executed.ipynb`;
- `results/validation/ont-pilot-retry-26174581/results/report/index.html`; and
- `tasks/ont-pilot-retry-26174581/completion-audit.json`.

Repository documentation is the reviewed delivery source. No hosted wiki
deployment is claimed. Exact commands, versions, database releases, failed
approaches, and resource evidence are recorded in the parent project
`tasks/METHODS.md`. Historical failed runs remain linked there as evidence for
the defects they exposed; they do not replace the accepted current-run counts.
