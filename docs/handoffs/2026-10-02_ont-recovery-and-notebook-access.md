# Handoff: ONT recovery and notebook access

**Date:** 2026-10-02
**Branch:** main
**HEAD:** 3f383ccf0edc5427b2b9c404ba372edba0b66964

## Context & Status

The active objective remains processing X0342–X0346 ONT reads with quality
filtering and a 500 bp minimum. X0342, X0345 and X0346 passed final audit 26470140.
X0343 and X0344 have completed characterization, including Symclatron, but lack
final catalogs and reports. Functional annotation is omitted from recovery at
the user's request; gene calling remains enabled.

## Technical Implementation

### Work Completed

- Verified all five native Symclatron tables against routed input identifiers.
  X0343 has 1,946 rows, including 36 predictions above the confidence threshold
  classified as host associated or obligate intracellular. X0344 has 783 rows
  and 20 such predictions. These are predictions, not confirmed symbioses.
- Diagnosed X0343 bin_2580: normal RepeatModeler exited zero with no families,
  explicit zero-family rounds and full 53,033 bp coverage. The wrapper wrongly
  accepted its blank numeric summary only in RECON-only recovery mode.
- Removed that single mode restriction in `modules/gene_calling.nf`, retaining
  exit, crash, model, marker and coverage checks. Native proof 26641649 completed
  0:0 in 1m28s. All 198 tests, style, schema and workflow previews passed.
  Initial review 26641656 timed out without a verdict. Bounded review 26641764
  completed 0:0 and approved the source fix and continuation drivers without
  blocking findings. Publish this fix and handoff together before recovery.
- Prepared task-local recovery under `tasks/ont500-finish-20261002/`.
  Fresh routing must reproduce original evidence bytes and route membership.
  Completed characterization is imported. X0343 reruns four gene producers;
  X0344 imports their outputs. Each case has a native MICRO gate.

### Outcomes

X0343 job 26512474 reached a 72h outer allocation timeout after the native workflow
had already failed. The old wrapper used direct `exec pixi`; the exact surviving
process state was not captured. New wrappers retain a Bash parent, following
the previously verified repair. X0344 job 26492222 failed at the 24h annotation
limit. Neither old job remains active. No cancellations were performed.

The strict documentation build failed only on links in the pre-existing,
untracked September 25 handoff. That file remains unchanged. Job 26641671
built tracked documentation plus this handoff with the unchanged strict
configuration and completed 0:0 in nine seconds.

### File Map

| Path | Purpose |
| --- | --- |
| `modules/gene_calling.nf` | Minimal zero-family correction |
| `tasks/ont500-finish-20261002/` | Review, native proof and recovery drivers |
| `../tasks/todo.md` | Authoritative progress, including Agentdash T-55414 |
| `../tasks/METHODS.md` | Commands, methods and evidence |

## Key Decisions

Preserve completed biological outputs. Do not rerun prokaryotic characterization
or the three accepted samples. Use existing annotation-skip handling and honest
recovery receipts. Preserve the RNA/BRAKER pending reasons when compatible
evidence is absent. Trust CheckEUK for eukaryotic completeness and lineage.

## Knowledge Capture

### Lessons Learned

Normal RepeatModeler can print a blank family count on a valid zero-family run.
The acceptance predicate requires explicit supporting evidence, not exit zero alone.

### Gotchas

The existing Cloudflare Jupyter service at <https://dori-nb.newlineages.com/lab>
serves `/clusterfs/jgi/scratch/science/mgs/nelli`. Public access redirects to
Cloudflare sign-in; authenticated public browser rendering was not tested.
Use `/lab/tree/` notebook links, because the local `/files/` HTML route returned 403.
No hosting configuration was changed.

X0342/X0345/X0346 notebooks at
`results/ont500-20260925/<sample>/results/report/report.executed.ipynb` each
contain eight executed code cells, three saved figures and no errors. The X0344
comparison notebook at
`tasks/ont-500-20260925/x0344-merge-assessment/assessment.executed.ipynb`
contains eight executed cells, two figures and no errors. It is a separate
analysis, not the missing final X0344 pipeline report.

## Moving Forward

### Next Steps

1. Commit and push the source fix and this handoff. Review 26641764 approves
   both and supersedes the timed-out review. Submission checks and the two
   Low findings are recorded in `tasks/ont500-finish-20261002/review-triage.md`.
2. Run the matching native MICRO gates and X0344/X0343 continuations from
   `tasks/ont500-finish-20261002/x0344-finalize/`. X0344 uses 4 CPUs and 32 GB;
   X0343 requires an explicit 32 CPU, 384 GB allocation override.
3. Audit final catalogs, notebooks, identifiers, recovery provenance and read QC.
   Combine those results with the accepted three-sample audit before completing
   the five-sample objective. Resolve the deferred owner biology-summary request.

### Blockers

No external blocker. Validation and both production recoveries remain unfinished.
