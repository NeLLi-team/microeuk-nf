# Handoff: Prodigal source batches

**Date:** 2026-10-07
**Branch:** detached from origin/main for validation
**HEAD:** 8f2c96b4e974a33eaeb84b8c25ef6bd027335b9d (validation base)

## Context & Status

Prodigal sources run in fixed batches capped at `task.cpus`. The process keeps
the existing scientific arguments, output paths, source ledgers and resource
requests. Each source runs in its own subshell. The parent waits for every
member of a batch, names failed sources and stops before launching another batch.

Native gate 26715500 passed. On 17 actual X0348 sources, task time was
386.154 seconds for the original serial module and 49.290 seconds for the
candidate, an 87.24% reduction. All output bytes matched. This was one
same-node comparison, not a full-workflow speedup estimate.

The original X0347 and X0348 workflows remain pinned to 08829ca. Do not pull
new source into that checkout until both original workflows and their
source-bound audits finish. The study task records own production status.

## Technical Implementation

### Work Completed

The module uses one Nextflow process and adds bounded background source calls.
The native regression runs through the actual module and managed environment.
The existing validator includes that regression without new dependencies.

### Outcomes

Gate 26715500 completed with exit 0:0 in 15 minutes 26 seconds. It passed
203 tests, Python style on 32 existing files, schema checks, strict docs,
Nextflow configuration and previews, and the native regression. Source hashes
matched before and after execution. No Python source changed.

Candidate execution at one versus two CPUs produced identical outputs for
mixed prokaryotic/viral fixtures in automatic and explicit-code modes. Checks
covered odd and exact batches, zero coding sequences, native Single and
Metagenomic modes, exact failure attribution, absent success metadata and no
later batch after failure. Event logs show command overlap and completion;
Nextflow's pipe handling means they do not independently prove parent wait calls.
Static review checked the explicit wait-all loop.

The original serial module failed the new concurrency check after successful
native execution. Its mixed-fixture outputs also matched the candidate in
automatic mode. The real-source benchmark compared original and candidate
automatic mode at 16 CPUs, with 17 sources and observed command-overlap peaks
of one and 16. No original-module explicit-code run is claimed.

Plan review 26714962 approved the bounded approach. Implementation review
26715117 found test weaknesses and had incomplete file access. Follow-up
26715502, using Claude Opus 5.5, accepted all seven repairs and reviewed the
previously unread scripts supplied in the prompt. Root resolved its remaining
syntax-check and equivalence-coverage findings with separate checks against
the saved outputs. No module defect was identified by either review.

### File Map

| File | Purpose |
|---|---|
| `modules/gene_calling.nf` | Source batches and failure propagation |
| `tests/check_prodigal_parallel.sh` | Native regression using one and two workers |
| `scripts/validate.slurm` | Existing validator with the regression and required tools |

## Key Decisions

Fixed batches keep the implementation small and preserve one Nextflow task.
They wait for the slowest source in each batch, so uneven source sizes can
leave cores idle. No worker pool, dependency or user option was added.

The benchmark selected the 17 smallest prepared sources in a 3–8 MB window;
their actual file sizes were about 3.01–3.35 MB. These are original routing
categories, not quality-qualified genomes. Serial execution ran first, so
filesystem caching may favor the candidate. The narrow size range, batch
stragglers and shared-node contention limit extrapolation. The acceptance
threshold of a 25% time reduction was set before the benchmark.

## Knowledge Capture

### Lessons Learned

The first gate, 26715116, stopped when recursive fixture copying could not
preserve directory modes. The repaired test uses `--no-preserve=mode`.
The original poly-A zero-CDS fixture produced a partial gene; a stop-rich
fixture now proves zero genes natively. `bash -n` accepts one script followed
by its arguments, so each script needs a separate syntax check.

### Evidence

Proof is retained under the canonical checkout's
`tasks/x0347-x0348-20261005/prodigal-parallel/`. `gate-26715500/` owns the
native run and measured script snapshot. `supplemental-checks-0909.txt` owns
the separate syntax, mixed-fixture equality and copied-input hash checks.
The review files and the final landing receipt are retained beside them.
Detailed commands and production history are in the study's `tasks/METHODS.md`.

## Moving Forward

### Next Steps

1. Complete both original workflows and their final source-bound audits.
2. Perform the recorded corrected-routing comparison and any required downstream
   replay, binding the selected source revision explicitly.
3. Finish and verify the final catalogs and reports before completing the goal.

### Blockers

The original workflows and final audits remain unfinished. Their live jobs
are the next source of results; no user input is needed.
