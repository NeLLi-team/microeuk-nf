# Handoff: ONT 500 bp five-sample completion

**Date:** 2026-10-02
**Branch:** main
**HEAD:** ecd68cf9b58671b938ec0e9cbd690737e22486b2

## Context & Status

X0342, X0343, X0344, X0345 and X0346 passed their recorded production audits.
Each run starts from the requested ONT reads, uses `chopper -q 10 -l 500`,
and assembles the filtered reads. The audits confirm an observed minimum read
length of 500 bp, SQLite integrity and zero foreign-key violations.

The [executed comparison notebook](https://dori-nb.newlineages.com/lab/tree/frederik/projects/08protists/protist-meta-nf/tasks/ont500-results-overview-20261002/overview.executed.ipynb)
shows all five samples as accepted. Its 13 code cells executed without errors,
and all four embedded figures passed visual review. Cloudflare Access sign-in
is required; authenticated browser rendering was not checked.

## Technical Implementation

### Work Completed

- X0343 production 26644346 completed all eight fresh finalization tasks with
  exit 0. The validated upstream stages retain their original provenance.
- Acceptance audit 26644349 passed after 21m32s. Its receipt checks original
  reads, native outputs, 16,849 artifacts, catalog relationships and 16 rendered
  report tables. The source remained frozen through the audit.
- The X0343 report contains nine executed cells, three reviewed figures and no
  errors. Its 1,831,070 Prodigal gene records and zero functional annotations
  match the requested gene-calling and annotation settings.
- Overview refresh 26647536 passed. Independent review confirmed unchanged
  input checksums, scientific selection logic and export counts. The preceding
  overview snapshot remains in the analysis directory's `pre-acceptance/`.

### Outcomes

All five runs meet the requested processing scope. The final X0343 recovery
uses the repair documented in the [recovery handoff](2026-10-02_empty-refiner-recovery.md).
Native RNA alignment and BRAKER states remain pending because compatible
evidence was not supplied. Successful execution does not confirm candidate
genomes, symbionts or viruses as biological findings.

### File Map

All paths below are relative to the repository root.

| Record | Location |
|---|---|
| Five-sample acceptance index | `tasks/ont500-finish-20261002/five-sample-acceptance.json` |
| X0342, X0345 and X0346 audits | `tasks/ont-500-20260925/final-audit-26470140/` |
| X0344 audit | `tasks/ont500-finish-20261002/final-audit/job-26641999/` |
| X0343 audit | `tasks/ont500-finish-20261002/final-audit/job-26644349/` |
| X0343 final output | `results/ont500-20260925/X0343-noannotation-20261002-refiner/` |
| Executed overview and validation | `tasks/ont500-results-overview-20261002/` |
| Exact methods and chronology | `../tasks/METHODS.md`, `../tasks/todo.md` |

## Key Decisions

| Decision | Rationale |
|---|---|
| Retain CheckEUK qualifications and host abstentions | Eukaryotic biological summaries follow Chef's CheckEUK-only constraint. Routing and other tools provide separately labeled context. |
| Preserve pending RNA/BRAKER states | Missing compatible evidence limits biological prediction, even when the workflow handles that state successfully. |
| Keep the existing notebook service | The executed files are accessible through the established Cloudflare route without a hosting change. |

## Knowledge Capture

### Lessons Learned

The audit establishes the requested scope through native output, provenance and
report reconciliation. A successful Slurm exit alone is insufficient evidence.
Earlier recovery decisions and proofs remain in the linked recovery handoff.

### Gotchas

- Generated analysis files are intentionally outside Git, matching run outputs.
  Their checksums and execution receipts are retained beside the notebook.
- The older queued job 26641846 was not cancelled because approval was absent.
  Preserve its shared finalizer and `nextflow.config`; do not treat an old
  cancellation question as approval to cancel a different job.
- Workflow acceptance does not imply complete eukaryotic genomes or confirmed
  symbiosis. The notebook retains tool-specific confidence and conflicts.

## Moving Forward

### Next Steps

No processing work remains for the five-sample goal. Use the executed overview
and its exported tables for biological review. Additional RNA/protein evidence
or functional annotation would be a separate analysis scope.

### Blockers

None for the requested five-sample processing and results delivery.
