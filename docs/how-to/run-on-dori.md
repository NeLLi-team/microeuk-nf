# Run and resume on Dori

Use `scripts/run-allocation.slurm` to run one sample in one Dori allocation.
Nextflow uses its local executor inside that allocation. The full workflow has
17 processes; the core workflow has seven. Process labels set CPU and memory
requests in `conf/resources.config`. The launcher reserves 24 GB for Nextflow
and allocation overhead and gives the remaining memory to the local executor.

## Prerequisites

- Run from a Dori login node in the workflow repository. The host must provide
  `/bin/bash` for Nextflow's outer task launcher.
- Prepare one sample row using the [input contract](../reference/inputs.md).
- Install the locked Pixi environments in Slurm, including Nextflow and Java.
  Finish `scripts/bootstrap.slurm` before running `scripts/install-tools.slurm`,
  `scripts/install-custom.slurm` or `scripts/install-dfam.slurm`.
  Register the reference paths in
  `conf/databases.yaml`. Referenced application workspaces must also be installed.
- Check `SUMMARY.md` in the repository for validation limits before production use.

The launcher validates registry entries and installed executable paths. Runtime
commands use `pixi run --as-is`; they do not install missing software. Environment
installation and reference setup must finish before submission.

## Submit a sample

1. Set the repository, sample sheet and a new run directory. Replace the example
   sample sheet with your own. The run directory must not exist for a fresh run.

   ```bash
   cd /clusterfs/jgi/scratch/science/mgs/nelli/frederik/projects/08protists/protist-meta-nf
   project_root=$PWD
   samples="$project_root/data/samples.tsv"
   run_dir="$project_root/results/my-sample"
   mkdir -p tasks/logs
   ```

2. Check the shared queue and the workflow's allocation ledger before each wave.
   Count standalone allocations, helpers, retries and array elements toward the
   limit of 100 planned tasks. Submit at most 20 new tasks per wave. Allow at most
   `min(30, max(0, 80 - other_running))` running tasks for this workflow.
   Account for pending tasks in other workflows that may start. Coordinate with
   other agents and projects. When
   the allowance is zero, wait before submitting. These are operator checks;
   the launcher does not enforce shared queue limits.

   ```bash
   squeue -M perceus-00 -u fschulz --array
   ```

3. Submit the full workflow with 32 CPUs, 384 GB and up to 72 hours. Capture the
   job identifier only after `sbatch` succeeds. Record it with the script, sample
   and log paths in the run ledger.

   ```bash
   receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
     -p dori --qos=jgi_normal --chdir="$project_root" \
     --cpus-per-task=32 --mem=384G --time=72:00:00 \
     scripts/run-allocation.slurm "$samples" "$run_dir" full) || exit 1
   job_id=${receipt%%;*}
   [[ $job_id =~ ^[0-9]+$ ]] || exit 1
   printf '%s\n' "$receipt" > "tasks/allocation-${job_id}.receipt"
   ```

   For core mode, replace `full` with `core` and use at least 216 GB. Both modes
   require at least 32 CPUs. Use a longer approved QOS if the required allocation
   time exceeds the normal-QOS limit. The allocation walltime covers the entire
   graph, rather than one process.

4. Inspect logs and accounting. Routine scheduler polls must be at least
   60 seconds apart. Match the job's working directory and name before
   attributing a result to this run.

   ```bash
   sacct -M perceus-00 -j "$job_id" \
     --format=JobID,JobName,State,ExitCode,Elapsed,WorkDir -P
   tail -n 40 "tasks/logs/allocation-${job_id}.out"
   tail -n 40 "tasks/logs/allocation-${job_id}.err"
   ```

5. Inspect the run's `executions/<job-id>.<suffix>/trace.tsv` and published
   stage manifests. A completed process may record a valid scientific skip.
   Confirm identifiers, tool-specific evidence and artifact provenance before
   interpreting the final report. Successful publication writes the following
   files beneath the run directory:

   - `results/catalog/protist-meta.sqlite`
   - `results/report/report.executed.ipynb`
   - `results/report/index.html`

   The notebook must contain saved outputs without cell errors. A scheduler
   `COMPLETED` state alone does not establish biological correctness.

## Resume the same run

Retain the run's `launch/`, `work/` and `prepared/` directories. Keep the source,
sample sheet, registry and referenced Pixi manifests and locks unchanged. The
launcher rejects an identity mismatch. Start a new run directory after changing
those inputs; do not edit the saved identity record.

Repeat the queue and allocation checks, then submit the same command with the
fourth argument `resume`:

```bash
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  --cpus-per-task=32 --mem=384G --time=72:00:00 \
  scripts/run-allocation.slurm "$samples" "$run_dir" full resume) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "tasks/allocation-${job_id}.receipt"
```

Nextflow reuses eligible completed tasks from the saved session. A fresh
execution directory records the resume trace and resource reports. Check which
tasks were cached and validate the published outputs again. See the
[input reference](../reference/inputs.md) for source-identity limits.
