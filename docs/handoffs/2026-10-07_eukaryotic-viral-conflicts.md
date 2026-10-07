# Handoff: preserve eukaryotic evidence in viral routing

**Date:** 2026-10-07
**Branch:** detached worktree for landing on main
**HEAD:** f3d894ddb861b7074b5e51ec8d583e8504979d4f

## Context & Status

GVClass `d_EUK` and linked whole-contig viral calls could send a bin to the
viral route when CheckEUK's host lineage was blank. The corrected rule retains
that mixed evidence as unresolved, with candidate class `conflicting`.
CheckEUK completeness alone does not assign taxonomy.

The X0347/X0348 production objective remains active. X0348's queued continuation
is held before execution; X0347 continues its original core stages. Their
execution checkouts and existing continuation bindings remain frozen.

## Technical Implementation

### Work Completed

- `src/protist_meta/routing.py` uses one independent-eukaryotic-evidence flag
  for the viral guard, conflict explanation and unresolved fallback.
- `tests/test_routing.py` checks the corrected route and conflict explanation
  in the existing parameterized case.
- `docs/explanation/design.md` states the evidence rule and its limits.

### Outcomes

The intended regression failed before implementation in step 26695248.43.
The Python house gate passed in .44. All 203 tests and strict documentation
passed in .45. Native X0348 preview .47 moved exactly 13 bins from viral to
unresolved/conflicting; the other 763 evidence rows stayed unchanged. All
776 identifiers and native evidence fields were preserved. Candidate counts
are 369 prokaryotic, 34 eukaryotic, 12 viral and 361 unresolved.
Independent Opus source review .46 found no blocking defect. Its tool could not
read the external Git snapshot; root verified the full diff and source hashes.
The review's low-priority SSU test gap concerns an unchanged guard condition,
as confirmed directly against the base commit.

Exact commands and receipts are retained under the canonical pipeline's
`tasks/x0347-x0348-20261005/noannotation-continuation/`, with source gate logs
in `routing-implementation/` and native validation in
`routing-preview-X0348-fixed-20261007/validation.json`.
The sibling study's `tasks/METHODS.md` owns the experiment narrative.
No configured type checker was run. Duplicate predicates and an unused
Signals field were replaced; no files or imports were deleted.

## Key Decisions

Positive eukaryotic inclusion and SSU parsing are unchanged. The prior broader
RF01960 proposal was cut after the 25-bin native audit found no RF01960 loci
without GVClass `d_EUK` among those viral candidates. QuickClade remains
contextual and cannot establish or veto a route.

Bin_3 remains unresolved, with native Pfam completeness 69.936% and unscorable
contamination. It is not promoted into the supported eukaryotic route by this
change. Without independent eukaryotic evidence, a whole-contig geNomad call
can still establish a viral candidate route. The native audit reports the
affected base fractions; remaining candidates are not validated viral genomes.

## Moving Forward

### Next Steps

1. Land the tested source correction on main and record its exact commit.
2. Bind and gate a fresh annotation-disabled continuation to that commit.
   Preserve the held job's source and task files until it is retired.
3. Complete both samples and verify native results, catalog and executed reports.
   The study's `tasks/todo.md` owns those completion checks.

### Blockers

Agentdash inbox and milestone posts intermittently return HTTP 500. Unsent
updates are recorded in the study's task record; local work continues.
