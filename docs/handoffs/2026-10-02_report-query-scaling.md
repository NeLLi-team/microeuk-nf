# Handoff: report query scaling

Date: 2026-10-02. Base revision: `1123a3570ceaeb282e7361cc5dce182cf4acad8a`.

## Status and evidence

X0344 recovery allocation 26641917 completed with exit 0:0. Its report contains
nine executed code cells, three saved figures and no errors. Final acceptance
audit 26641999 is running against the unchanged canonical checkout. X0342,
X0345 and X0346 already passed their separate final audit.

The X0344 notebook's final cell took 775.879 seconds. Its stage summary joined
both child tables before counting distinct values. The `prodigal_gv` stage has
792 inputs and 990,550 expected keys, producing 784,515,600 intermediate pairs.
The repair counts each child table through a separate indexed subquery. It
preserves the existing expressions, stage identities, columns, ordering and
empty-stage zeros. It changes no scientific output or workflow stage choice.

Plan review 26642063 approved the approach. Gate 26642134 passed 199 tests,
Python style, schema checks, strict documentation and both workflow previews.
The new regression fails on the old query with SQLite's instruction-limit
interrupt and passes on the repaired query. The first native rerender reduced
the final cell to 31.979 seconds. Its comparator incorrectly included rerun
paths and timestamps among invariant fields. The corrected verifier validates
those fields separately and retains exact comparisons for result tables and
figures. Native verification in 26642258 passed: all 16 result tables and three
figures match, stable provenance fields match, and catalog bytes are unchanged.
The full report built in 86.160 seconds; its final cell took 36.629 seconds.
Independent diff review in the same allocation returned CLEAN, with no defects,
deletions or tidiness findings. Its file access covered the production change
and regression; the task-local native comparator was outside its allowed tree.

Evidence lives in `tasks/report-scaling-20261002/`. Current production status
and exact commands remain in `../tasks/todo.md` and `../tasks/METHODS.md`.

## Decisions and next steps

Keep canonical source frozen until audit 26641999 ends. The repair is ready
for publication after the final handoff documentation check. Then start the
tested X0343 continuation without functional annotation.

The old X0343 job 26641846 remains queued. Cancellation approval is pending,
but cancellation is not a prerequisite for the corrected, already authorized
run. The old immutable script uses the configuration pattern that failed its
mandatory MICRO gate in 26641869, before production setup. Preserve canonical
`nextflow.config`, the shared continuation driver and pinned Nextflow 25.10.2
until that old job terminates or cancellation is approved. Directory guards
alone would not prevent a simultaneous-start race.

Do not complete the five-sample goal until X0343 and X0344 both pass their
production audits. All runs retain ONT read QC and the 500 bp minimum length.
