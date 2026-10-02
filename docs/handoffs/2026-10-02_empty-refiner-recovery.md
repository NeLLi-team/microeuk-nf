# Handoff: RepeatScout empty refinement recovery

Date: 2026-10-02. Base revision: `320659b`.

## Failure and repair

X0343 allocation 26642638 failed after 39 minutes. Twenty repeat-masking bins
passed before bin_487 exposed a second RepeatModeler edge case. Its four
filtered RepeatScout families had four, three, four and four usable instances,
below the native minimum of five. RepeatModeler 2.0.9 scheduled no refinement
jobs, then searched a database for an empty refined-cons.fa and exited 2.
Nextflow aborted the parallel Prodigal task. The dependent audit 26642654 was
cancelled automatically by its invalid-dependency policy.

Native reproduction 26643710 confirmed the cause with the original random
seed. A separate RECON-only run recovered one 27,566-base repeat family, so
the first failure cannot be accepted as successful zero-repeat discovery.

The module now permits the existing RECON retry for this verified condition.
It checks the native failure signature, filtered and active consensus, family
instance counts and the empty refined consensus. Negative-start ranges and
ranges for filtered-out families do not contribute to eligible counts.
Every retained family needs valid range evidence and fewer than five eligible
instances. Missing retained families, malformed evidence and child crashes
remain fatal. Initial logs and the failed directory are kept. Classified-library,
zero-family and sequence-identity checks retain their existing behavior.

## Verification

Final native Nextflow proof 26644135 passed for bin_487 and prior positive,
no-seed and normal-zero controls. The recovered bin has one classified family
and 69,750 masked bases; normalized sequence content and identifiers are intact.
Two additional recovery controls passed. Six adverse cases failed before retry
or publication. Input and source hashes remained unchanged.

Final gate 26644136 passed all 199 tests, the 32-file Python style check, schema
checks, configuration parsing, both workflow previews, strict documentation
build and source integrity. Earlier verification-script corrections and their
original receipts remain in the task evidence; they changed no production code.

Review 26643859 found a fail-closed case with observed ranges but zero usable
instances. The guard now tracks observed families separately from parsed counts.
Fix-verification review 26644039 returned CLEAN. A later one-clause change
rejects empty sequence names; root checked it directly against the native
parser, and final proof 26644135 includes that adverse case. The external
review covers the substantive guard changes before that final tightening.

Final module SHA256:
`6fbef35242fd99585dbdc126c0457242dbf917b95fdb2285b4ede879945febbd`.

## Recovery and remaining work

Use the prepared wrapper at
`tasks/repeat-empty-refiner-20261002/recovery/start.slurm` after publication and
its mandatory MICRO gate. The fresh target is
`results/ont500-20260925/X0343-noannotation-20261002-refiner/`.
It preserves completed upstream analyses and Symclatron by importing the same
validated stages. Functional annotation stays disabled. Gene calling remains
enabled. Submit a fresh success-dependent production audit.

The plan review removed an unnecessary in-place resume design: forced fresh
routing invalidates downstream task paths, and neither failed gene task has
completed work to cache. The failed output directory stays intact.

Proofs and review records are under `tasks/repeat-empty-refiner-20261002/`.
Current commands, production job IDs and final outcomes belong in
`../tasks/METHODS.md` and `../tasks/todo.md`. Freeze tracked source during the
new production run and its audit. The older queued job 26641846 remains
untouched; preserve its canonical configuration and shared finalizer until
it terminates or cancellation is approved.

X0342, X0344, X0345 and X0346 are accepted. X0343 still needs successful
production and its acceptance audit before the five-sample goal is complete.
The executed cross-sample results notebook remains available under
`tasks/ont500-results-overview-20261002/`, with X0343 marked pending acceptance.
