# Lab notebook

Add dated entries with the goal, command or driver, QC result, interpretation, and next check.

## 2026-09-21: final ONT acceptance

**Goal.** Validate the repaired retained-input ONT workflow through initial
execution, cached resume, catalog inspection, report publication, and the final
completion audit.

**Driver.** Job 26174581 ran sample `X0341_p01` as
`run_823b4b2ea3996922` at source revision
`sha256:4f734784d5e5623955354126ff3a8e972f6047c126d5b1caffa951246322581e`.
The production wrapper and exact commands are recorded in
`../../tasks/METHODS.md`.

**QC result.** The job completed with exit 0:0 in 11:26:56. All 17 initial
Nextflow tasks completed, and all 17 were cached on resume. Inspection passed
with SQLite integrity `ok`, zero foreign-key violations, and 1,448 artifact
hashes checked. The completion audit matched all 18 catalog counts and 16 report
tables. The catalog contains 150 bins, 1,196 contigs, 616 QC rows, 1,300
taxonomy rows, 57 SSU loci, 29 viruses, 143 phenotypes, 90,992 genes, and
1,680,663 annotations. The notebook executed eight code cells with three
figures and no errors; the HTML export and visual review passed.

**Interpretation.** The current ONT workflow passed full scientific and delivery
acceptance. All six eukaryotic analysis bins lacked compatible RNA or protein
evidence, so their RNA and BRAKER stages remain explicit pending scientific
states rather than positive gene-call evidence. The separate eukaryotic case
retains the positive RNA-supported BRAKER proof. SSU assignments remain
`Unclassified` without runtime calibration, and GVClass evidence remains
exploratory.

Repository documentation is the reviewed delivery source. No hosted wiki
deployment is implied.
