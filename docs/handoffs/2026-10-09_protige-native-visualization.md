# Handoff: Protige native visualization bundles

**Date:** 2026-10-09
**Branch:** main
**HEAD before changes:** 0951871

## Status

The approved native exporter, CLI and Nextflow integration are implemented and
verified. Root owns final integrated review and permission to commit and push.
No commit, push or deployment has been made by this worker. Detailed commands,
failed attempts and receipts are in `tasks/protige-visualization-20261009/`.

## Changes and decisions

Added `src/protist_meta/visualization.py` and its focused integration tests.
The CLI accepts `visualization export`; `BUILD_VISUALIZATION` runs after catalog
completion beside the report process. README and report reference describe the
schema. No files were deleted, no dependencies added and no catalog schema changed.

Bundles preserve native catalog rows, QC models and registered evidence tables.
Website builds assign display categories and registry acceptance. Mapping uses
digest-checked coverage, the native read_qc retained-read denominator and recorded
mapping policy. Missing coverage is nullable. Missing registered QC or taxonomy
reports from completed stages fail closed. Partial inventories do not claim a
whole-assembly residual. Sequence payloads and gene/annotation rows are excluded.
Public consumers must omit private source paths and commands.

## Verified results

- Thirteen focused exporter tests and Python style pass. The full repository
  suite passed: 216 tests in 34.99 seconds, measured 36.91 seconds with one CPU.
- Strict documentation build passed from the intended source tree in 7.37
  seconds. The initial working-tree build failed only on ten broken links in
  the unrelated untracked September 25 handoff. That file remains untouched.
- Both accepted native catalogs passed bounded local CLI export and count
  verification. X0348 took 19.53 seconds, X0347 30.67 seconds. The accepted
  marker is `local-20261010T010929Z/local-proof-complete.json` under the task
  directory. It binds source, command, execution, index and verification receipts.
- Actual Nextflow `BUILD_VISUALIZATION` for X0348 passed in the bounded local
  process proof `local-process-20261010T011725Z`: 25.16 seconds for the full
  command and checks, peak RSS 359212 KiB. The trace contains exactly one
  completed process, exit 0. Its bundle digest matches the accepted CLI output.
  Java, process shells and Python used CPU 1. All source pins passed afterward.

The first full-suite attempt timed out at 60 seconds with temporary files on
GPFS. The one approved retry used node-local temporary storage and passed.
Logs for both attempts and intended-source strict docs are in
`local-checks-20261010T011209Z` and `local-checks-20261010T011459Z`.

## Remaining work

Root must complete integrated review before landing the named source changes.
The backend worker has accepted both local bundle paths. Legacy bin-read counts
remain a separate pending job, 26795197; native proof does not validate them.

Dori job 26799930 remains a separate queued full validation attempt. Its outputs
are isolated under `job-26799930`. Older jobs were not canceled or modified.
Job 26796904 has a known global `-C` argument-order error and a ten-minute limit;
26799930 has the corrected order and a thirty-minute limit. Local proof above
is real CLI/Nextflow execution and does not claim a completed Slurm run.
