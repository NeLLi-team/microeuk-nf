# Handoff: QuickClade routing correction

**Date:** 2026-10-06
**Branch:** detached worktree for delivery to origin/main
**HEAD:** b2bb8286caff1f7f2a599d6aa43bc5ab0f1c4f25

## Context & Status

Chef requires eukaryotic-bin assessment from CheckEUK without relying on
QuickClade. All bins already receive CheckEUK, but downstream routing allowed
QuickClade alone to promote unresolved-host candidates into a prokaryotic route.
The correction is implemented and tested. Opus 5.5 review 26706827 accepted
the fixes with no blocking defects.

The canonical Dori checkout remains at 08829ca for active X0347 and X0348
production jobs 26695248 and 26695249 and their source-bound audit 26695842.
Do not pull remote main into that checkout before those runs and audits finish.
Current scientific status and exact methods have one home in the sibling study:
`../09enrichments-ont/tasks/todo.md` and `../09enrichments-ont/tasks/METHODS.md`.

## Technical Implementation

### Work Completed

- `src/protist_meta/routing.py` requires GVClass support for a prokaryotic route.
  QuickClade labels remain in the evidence table but cannot select or veto a route.
- `tests/test_routing.py` covers unresolved CheckEUK hosts and independently
  supported viral routes, including retained native labels and emitted FASTA.
- `docs/explanation/design.md` and `docs/assets/workflow.svg` describe the rule.

### Outcomes

Gate 26706802 passed 203 tests, the two-file Python style check, strict docs,
and diff checks. All four new regression cases fail against the original source.
The earlier implementation review found a documentation qualification, a stale
diagram label, and a missing linked-virus case; these were corrected before
the final gate. The follow-up review accepted those fixes; the SVG change alters
only its text label, with all geometry preserved. Proof is under
`tasks/x0347-x0348-20261005/quickclade-routing/` in the canonical Dori checkout.

### File Map

The four changed source, test and documentation files are listed above. The
study Methods record has the native bin_3 example, job receipts and full review
provenance. This handoff records how to resume the correction and production work.

## Key Decisions

Keep the existing combined eukaryotic-support rule and output schema. A high
CheckEUK completeness score alone does not assign taxonomy. Preserve unresolved
lineage and unscorable contamination rather than inventing accepted MAGs.
The correction is not applied to the running production source.

## Knowledge Capture

### Lessons Learned

All-bin screening does not prevent a later routing rule from hiding a candidate.
The regression now checks that QuickClade cannot independently determine routing.

### Gotchas

After the original-source audits, compare complete production routing evidence
with the corrected rule before deciding which downstream outputs need replay.
The live GVClass results were incomplete when this correction was prepared.

## Moving Forward

### Next Steps

1. Continue both original jobs and inspect their final source-bound audit.
2. Measure actual routing changes from complete native evidence and replay only
   affected downstream work with an explicitly recorded corrected source.
3. Preserve all-bin CheckEUK results in the final assessment and finish the
   original end-to-end goal only when both samples and final outputs pass.

### Blockers

None. Production is running; its final results and audits remain unfinished.
