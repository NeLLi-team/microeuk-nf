# Run and resume one sample on Dori

`scripts/run-allocation.slurm` runs one sample inside one Dori allocation.
Nextflow uses its local executor with `process.maxForks = 1`. The full workflow
has 17 processes with default options, 16 with annotation disabled, and 12 with
gene calling disabled; the core workflow has seven. Process labels in
`conf/resources.config` set per-process CPU, memory, and time requests. The
launcher reserves 24 GB for Nextflow and allocation overhead and gives the
remaining memory to the local executor.

## Timeout retries

`READ_QC` and the nested SSUextract `BLAST_ANNOTATE` process get one retry when
Nextflow reports a process running-time timeout. The first attempt requests
eight hours and the second requests 16 hours. Other `READ_QC` errors terminate
the workflow. Nested BLAST annotation retains its existing retries for exit
code 104 and signal-style exit codes 130 through 145; other errors use its
existing `finish` strategy.

The nested workflow receives `--max_time` from the outer `SSU_EXTRACT` process.
This value caps the request for each nested process attempt. The outer process
limit must cover the eight-hour first attempt, the 16-hour retry, and all other
nested stages. The `heavy_qc` label currently gives `SSU_EXTRACT` 48 hours. The
Slurm allocation walltime does not raise either process limit.

Judge whether these limits are sufficient from task logs and measured
production runtime.

## Clone the repository

Run these commands from a Dori login node:

```bash
git clone git@github.com:nelli-team/microeuk-nf.git
cd microeuk-nf
project_root=$PWD
mkdir -p "$project_root/tasks/logs"
```

The host must provide Pixi, Slurm, and `/bin/bash`. The operating-system Bash
starts the bootstrap, which re-executes inside the repository's root Pixi Bash.
Nextflow's outer task launcher also invokes host `/bin/bash`. Nextflow, Java,
Python, the inner process shell, and all workflow tools are Pixi dependencies.

## Prepare environments and references

The repository does not provide a standalone full installation. Configure
`conf/databases.yaml` for the target clone location and filesystem. Its relative
paths are resolved from `conf/`, so moving the clone can change the referenced
targets.

Core mode needs the repository environments and the registered QuickClade
reference. Full mode also needs these external workspaces:

- the CheckM1, CheckM2, GTDB-Tk, and Symclatron environments in
  `checkm_manifest`;
- the geNomad and CheckV environments in `viral_manifest` and
  `checkv_manifest`; and
- the CheckEUK, GVClass, and SSUextract application workspaces.

Install those external workspaces and the reference resources for enabled stages
separately. The repository scripts do not provision them.

Before every setup or analysis submission wave, check the shared queue:

```bash
squeue -M perceus-00 -u fschulz --array
```

Count standalone allocations, helpers, retries, and array elements toward the
limit of 100 planned tasks. Submit at most 20 new tasks per wave. Allow at most
`min(30, max(0, 80 - other_running))` running tasks for this workflow,
including pending tasks that can start. When the allowance is zero, wait before
submitting. The setup scripts and allocation launcher do not enforce these
shared limits.

Submit the repository setup jobs in order, waiting for each job to succeed
before submitting the next:

```bash
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  scripts/bootstrap.slurm) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/bootstrap-${job_id}.receipt"
```

If the registered Dfam 4.0 directory is absent, submit
`scripts/install-dfam.slurm` next. Skip this installer when that release is
already present and valid because it refuses to replace an existing release.
See [Install Dfam](install-dfam.md) for its checks.

After Dfam is available, install the remaining repository tool environments,
then the custom BRAKER3 and InterProScan packages. Wait for the tool job to
succeed before submitting the custom-package job.

```bash
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  scripts/install-tools.slurm) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/tools-${job_id}.receipt"
```

```bash
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  scripts/install-custom.slurm) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/custom-${job_id}.receipt"
```

The last two scripts use fixed local paths. `install-tools.slurm` probes a
fixed eggNOG data directory. `install-custom.slurm` uses fixed InterProScan data
and smoke-test paths. The InterProScan recipe also requires the local source
distribution named by its source hash manifests and
`.nellidb_version/version.txt` metadata for version `5.76-107.0` and MD5
`8c9a8b153e527f8cfc7bf24ee1652d78`. Changing only
`conf/databases.yaml` does not redirect these setup checks.

The setup commands above install both optional stages. For an installation
that omits them, use the [dependency contract](../reference/inputs.md#dependency-registry)
to identify the required environments and resources. The setup scripts do not
accept the runtime skip switches.

The allocation launcher validates every enabled registry entry and executable.
Runtime commands set `PIXI_NO_INSTALL=true` and `PIXI_FROZEN=true` and use
`pixi run --as-is`. They do not install missing software or update locks.

## Prepare one sample

Copy the 13-field template to a file of your choice and edit its single example
row. The allocation launcher currently accepts exactly one sample row. Input
paths are resolved relative to the copied sample sheet.

```bash
samples=/absolute/path/to/run-inputs/samples.tsv
mkdir -p "$(dirname "$samples")"
cp "$project_root/data/samples.example.tsv" "$samples"
# Edit "$samples" before submission.
```

See the [input and dependency contract](../reference/inputs.md) for field
meanings and evidence-pair requirements.

## Submit a fresh run

1. Choose a new run directory. It must not exist before a fresh run.

   ```bash
   run_dir="$project_root/results/my-sample"
   ```

2. Repeat the shared queue and task-budget check from the setup section. Do not
   submit while the workflow's running-task allowance is zero.

   ```bash
   squeue -M perceus-00 -u fschulz --array
   ```

3. Submit full mode with 32 CPUs, 384 GB, and up to 72 hours. Capture the job
   identifier only after `sbatch` succeeds.

   ```bash
   receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
     -p dori --qos=jgi_normal --chdir="$project_root" \
     --cpus-per-task=32 --mem=384G --time=72:00:00 \
     scripts/run-allocation.slurm "$samples" "$run_dir" full) || exit 1
   job_id=${receipt%%;*}
   [[ $job_id =~ ^[0-9]+$ ]] || exit 1
   printf '%s\n' "$receipt" > "$project_root/tasks/allocation-${job_id}.receipt"
   ```

   For core mode, replace `full` with `core` and request at least 216 GB.
   Both modes require at least 32 CPUs. Use `--qos=jgi_long` when the required
   walltime exceeds the verified normal-QOS limit. The allocation walltime
   covers the full process graph.

4. Inspect logs and accounting. Keep routine scheduler polls at least 60 seconds
   apart. Match the job's working directory and name before attributing the
   result to this run.

   ```bash
   sacct -M perceus-00 -j "$job_id" \
     --format=JobID,JobName,State,ExitCode,Elapsed,WorkDir -P
   tail -n 40 "$project_root/tasks/logs/allocation-${job_id}.out"
   tail -n 40 "$project_root/tasks/logs/allocation-${job_id}.err"
   ```

5. Inspect `executions/<job-id>.<suffix>/trace.tsv` and the stage
   manifests below the run directory. A completed process can record a valid
   scientific skip. Confirm identifiers, evidence, and artifact provenance.
   A successful run writes:

   - `results/catalog/protist-meta.sqlite`
   - `results/report/report.executed.ipynb`
   - `results/report/index.html`

   The notebook must contain saved outputs without cell errors. Slurm
   `COMPLETED` alone does not establish scientific correctness. Review the
   [source-scoped validation
   summary](https://github.com/nelli-team/microeuk-nf/blob/main/SUMMARY.md)
   before interpreting production results.

## Choose gene calling and annotation

Append `--skip-annotation` to the launcher arguments to retain gene models
without functional annotation. Append `--skip-gene-calling` to omit both.
The rest of full mode, including bin quality assessment and classification,
still runs. Keep the same allocation minimums because these analyses remain
enabled.

For example, after preparing the sample sheet and checking the queue, submit
a fresh run without either optional stage:

```bash
run_dir="$project_root/results/my-sample-no-genes"
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  --cpus-per-task=32 --mem=384G --time=72:00:00 \
  scripts/run-allocation.slurm "$samples" "$run_dir" full \
  --skip-gene-calling) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/allocation-${job_id}.receipt"
```

The catalog, executed notebook, and HTML report are still produced. Their stage
tables identify the disabled analyses. No gene models or functional annotation
results are reported for those stages.

## Resume the same run

Retain the run's `launch/`, `work/`, and `prepared/` directories. Keep the
source, sample sheet, registry, and referenced Pixi manifests and locks
unchanged. The launcher rejects an identity mismatch. Changed inputs require a
new run directory. A reviewed source repair can use the separate procedure below.

Repeat any skip switches from the fresh run after `resume`. For example, a run
started with `--skip-annotation` requires `full resume --skip-annotation`.
Changing the switches requires a new run directory. The launcher rejects a
scope change before starting Nextflow.

Repeat the queue checks, then submit the same run with the fourth argument
`resume`:

```bash
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  --cpus-per-task=32 --mem=384G --time=72:00:00 \
  scripts/run-allocation.slurm "$samples" "$run_dir" full resume) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/allocation-${job_id}.receipt"
```

Nextflow reuses eligible tasks from the saved session. A new execution directory
records the resume trace and resource reports. Check cached task states and
validate the published outputs again.

## Resume after a reviewed source repair

Use this procedure only after reviewing the exact source changes and confirming
which completed scientific stages remain valid. It does not determine whether
arbitrary code or dependency updates are compatible with cached results.

Keep the original skip switches after the two source digests. A reviewed source
repair does not permit changing the run's stage selection.

The launcher reruns `ROUTE_BINS`, `COLLECT_RECORDS`, `BUILD_CATALOG` and
`BUILD_REPORT`. Nextflow decides cache eligibility for other tasks from its
normal task hashes. A changed task script can therefore repeat scientific work.
In particular, changing `FULL_PREFLIGHT` can also invalidate downstream tasks
because its output is an input to read processing. Do not assume a source repair
will retain the scientific cache.
Changes to shared Python helpers or external tools need their own cache review;
the four forced tasks alone do not cover such changes.

1. Wait for every run using this checkout to stop before updating its source.
   Retain the same checkout path, input files and saved run directories.
   Keep `prepared/` and its identity file unchanged.

2. After reviewing the repair, record the original and repaired content digests.
   The repaired checkout and referenced dependency manifests must remain fixed
   during execution.

   ```bash
   original_digest=$(awk -F '\t' '$1 == "source_revision" {print $2}' \
     "$run_dir/prepared/resume-identity.tsv")
   repaired_digest=$(pixi run --as-is --manifest-path "$project_root/pixi.toml" \
     python -c 'import sys; from pathlib import Path; from protist_meta.inputs import _source_digest; print(_source_digest(Path(sys.argv[1])))' \
     "$project_root")
   printf 'Original: %s\nRepaired: %s\n' "$original_digest" "$repaired_digest"
   ```

3. Check the shared queue, then submit the explicit source transition.

   ```bash
   squeue -M perceus-00 -u fschulz --array
   receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
     -p dori --qos=jgi_normal --chdir="$project_root" \
     --cpus-per-task=32 --mem=384G --time=72:00:00 \
     scripts/run-allocation.slurm "$samples" "$run_dir" full resume-reviewed \
     "$original_digest" "$repaired_digest") || exit 1
   job_id=${receipt%%;*}
   [[ $job_id =~ ^[0-9]+$ ]] || exit 1
   printf '%s\n' "$receipt" > "$project_root/tasks/allocation-${job_id}.receipt"
   ```

4. Verify the new execution trace and outputs. All input identity checks still
   apply. The execution directory contains `recovery.json`, which records both
   source digests, the original preparation identity, the allocation configuration
   digest and the forced tasks. The report copies that receipt and displays its
   contents and SHA-256 in the executed notebook and HTML.

If this recovery attempt fails, repeat `resume-reviewed` with the same two
digests. Each attempt writes its own receipt. Plain `resume` still compares
against the original preparation source and rejects the repaired source.

The catalog retains the original preparation source identity so that existing
task metadata and caches remain comparable. The displayed recovery receipt
identifies the repaired execution source. Interpret these together; the original
catalog identity alone does not describe the code used for recovered outputs.

## Check the timeout retry policy

After installing the repository environments, submit the standalone retry
check in a small allocation. It tests timeout retries and ordinary exit-code
failures with both installed Nextflow engines. Run the 32-CPU MICRO workflow
gate in a separate allocation.

```bash
squeue -M perceus-00 -u fschulz --array
proof_dir="$project_root/tasks/resource-retry-$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$project_root/tasks/logs"
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  --cpus-per-task=2 --mem=16G --time=00:10:00 \
  --output="$project_root/tasks/logs/resource-retry-%j.out" \
  --error="$project_root/tasks/logs/resource-retry-%j.err" \
  --wrap="bash tests/check_resource_retries.sh '$proof_dir'") || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt"
```

After the allocation completes, `$proof_dir/status.txt` contains `PASS`.
