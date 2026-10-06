# Handoff: SSU BLAST thread allocation

**Date:** 2026-10-06
**Branch:** detached worktree for delivery to origin/main
**HEAD:** f69124c103e7223251f1653b803c38f24e83e7fd

## Context & Status

The nested SSU workflow now allocates up to eight CPUs to BLAST, capped by
the parent task's available nested pool. The one-line change passed native
resource and output checks, 203 project tests, strict docs and Opus 5.5 review.

The canonical Dori checkout remains at 08829ca for production jobs 26695248
and 26695249 and source-bound audit 26695842. Keep that checkout pinned until
the runs and their audits finish. Current sample status and exact methods are
maintained in the sibling study's `tasks/todo.md` and `tasks/METHODS.md`.

## Technical Implementation

### Work Completed

`modules/characterization.nf` sets the existing `BLAST_ANNOTATE` selector to
`Math.min(8, nestedCpus)`. CMSEARCH retains four threads. Scientific parameters,
reference data and retry behavior are unchanged.

### Outcomes

Benchmark 26705890 produced byte-identical BLAST outputs for both SSU models
with eight threads. Its observed runtime reduction is relative to historical
four-thread runs, with node-load and cache differences; it is not a controlled
whole-workflow speedup estimate. Timings and comparison proof are in the study
Methods record and `tasks/x0347-x0348-20261005/ssu-threading/job-26705890/`.

Native wrapper gate 26707090 verified parent allocations of 12 and 5 CPUs.
The nested pools had 11 and 4 CPUs, BLAST used 8 and 4 threads, and CMSEARCH
used 4 threads in both cases. Both nested executors had capacity one. Summary
tables, extracted FASTA, BLAST M8 and top-hit outputs matched exactly. All
203 tests and strict docs passed with the candidate source unchanged.

The parent-12 run completed in 26707002, but its checker falsely matched a
resolver's `--cmsearch` argument. The corrected checker anchors executable
names. Gate 26707090 reused that completed native run, checked the same source
digest, and ran only the remaining parent-5 case before the final checks.
Both receipts and review records are preserved under
`tasks/x0347-x0348-20261005/ssu-threading/implementation/` in the canonical checkout.

### File Map

- `modules/characterization.nf`: one resource-allocation line.
- This handoff: validation evidence, limits and production resumption rules.

## Key Decisions

Cap BLAST at the available nested CPU pool and retain the established scientific
configuration. The short wrapper proof uses the 18-fragment RF01960 input;
the separate benchmark covers the longer RF00177 workload. This change does
not supply the missing SSU runtime taxonomic calibration.

## Knowledge Capture

### Lessons Learned

Verify generated nested configuration and native commands through the actual
outer wrapper. Source inspection alone does not prove effective thread counts.

### Gotchas

Temporary `.pixi` symlinks used the canonical frozen environment. Their paths
are recorded in `implementation/runtime-aliases.txt`; remove only those aliases
after validation consumers finish. Preserve the canonical environment and proof.

## Moving Forward

### Next Steps

1. Continue both original production jobs and inspect their source-bound audit.
2. Compare complete routing evidence against the QuickClade correction f69124c
   and replay affected downstream work with recorded source provenance.
3. Validate final catalogs and executed reports before completing the goal.

### Blockers

None. Production results and final audits remain unfinished.
