# Run and resume one sample on Dori

`scripts/run-allocation.slurm` runs one sample inside one Dori allocation.
Nextflow uses its local executor. `process.maxForks = 1` allows one running
task per process, so tasks of different processes overlap while the executor
has free CPUs and memory. The launcher sets the executor limits to the
allocation's CPU count and to its memory minus a 24 GB reserve for Nextflow and
allocation overhead.

## Process graph

The number of processes depends on the run mode and the optional stages:

| Run | Processes |
|---|---:|
| Core mode | 8 |
| Full mode with default options (gene calling off) | 21 |
| Full mode with `--enable-gene-calling` | 25 |
| Full mode with `--enable-gene-calling --skip-annotation` | 24 |
| Full mode with `--enable-gene-calling --enable-repeat-masking` | 26 |
| Full mode with `--enable-gene-calling --enable-repeat-masking --skip-annotation` | 25 |

Core mode runs a preflight check, `READ_QC`, `ASSEMBLE`, `MAP_DEPTH_QUICKBIN`,
and four collection processes (`COLLECT_RECORDS`, `BUILD_CATALOG`,
`BUILD_REPORT`, `BUILD_VISUALIZATION`). Full mode replaces the core preflight
with `FULL_PREFLIGHT` and adds 13 processes:

- Assembly consumers: `SSU_EXTRACT`, `CONTIG_CLASSIFY` (Tiara), and
  `VIRAL_SCREEN` (geNomad, then CheckV).
- Bin screens: `GVCLASS_SCREEN`, `CHECKEUK_SCREEN`, and `CHECKM2_SCREEN`, which
  run as separate processes. The merge process `ALL_BIN_SCREEN` builds the
  published `04_all_bin_screen` directory from their outputs.
- `ROUTE_BINS`, which needs the merged screen, SSU, viral, and Tiara results.
- Routed-bin analyses: `EUKARYOTE_QC` (BUSCO, then tidk), `CHECKM1_QC`,
  `GTDBTK_CLASSIFY`, and `SYMCLATRON_CLASSIFY`. The merge process
  `PROKARYOTE_CHARACTERIZATION` builds the published `08_prokaryote_qc`
  directory from the last three.

`--enable-gene-calling` adds `MAP_RNA_EVIDENCE`, `BRAKER3_EUKARYOTES`,
`PRODIGAL_GV_GENES`, and `FUNCTIONAL_ANNOTATION`. `--skip-annotation` removes
`FUNCTIONAL_ANNOTATION`. `--enable-repeat-masking` adds
`REPEAT_MASK_EUKARYOTES`.

Tool processes that feed a merge process do not publish. Each saves its own
`.command.sh` as `<tool>/command.sh`, and the merge process records that file in
`stage.json`.

`READ_QC` reads the input once. It decompresses the file, tees the stream into
`seqkit stats -a` for the raw table, filters with Chopper (ONT) or SeqKit
(HiFi), tees the filtered stream into a second `seqkit stats -a`, and compresses
with pigz. A failed statistics job fails the task.

## Resource labels

Process labels in `conf/resources.config` set the CPU, memory, and time request
of each process.

| Label | Processes | CPUs | Memory | Time |
|---|---|---:|---:|---:|
| `qc` | `CORE_PREFLIGHT`, `FULL_PREFLIGHT`, `ROUTE_BINS` | 4 | 16 GB | 2 h |
| `read_qc` | `READ_QC` | 16 | 32 GB | 6 h, then 12 h on retry |
| `assembly` | `ASSEMBLE` | 32 | 256 GB | 24 h |
| `mapping` | `MAP_DEPTH_QUICKBIN` | 16 | 96 GB | 12 h |
| `ssu` | `SSU_EXTRACT` | 8 | 48 GB | 48 h |
| `contig` | `CONTIG_CLASSIFY` | 8 | 64 GB | 12 h |
| `viral` | `VIRAL_SCREEN` | 12 | 96 GB | 24 h |
| `gvclass` | `GVCLASS_SCREEN` | 20 | 96 GB | 48 h |
| `checkeuk` | `CHECKEUK_SCREEN` | 12 | 96 GB | 48 h |
| `checkm2` | `CHECKM2_SCREEN` | 8 | 64 GB | 48 h |
| `checkm1` | `CHECKM1_QC` | 16 | 96 GB | 48 h |
| `gtdbtk` | `GTDBTK_CLASSIFY` | 32 | 256 GB | 48 h |
| `symclatron` | `SYMCLATRON_CLASSIFY` | 16 | 64 GB | 48 h |
| `eukqc` | `EUKARYOTE_QC` | 16 | 96 GB | 24 h |
| `gene` | `REPEAT_MASK_EUKARYOTES`, `MAP_RNA_EVIDENCE`, `BRAKER3_EUKARYOTES`, `PRODIGAL_GV_GENES` | 16 | 160 GB | 48 h |
| `annotation` | `FUNCTIONAL_ANNOTATION` | 16 | 128 GB | 24 h |
| `report` | `ALL_BIN_SCREEN`, `PROKARYOTE_CHARACTERIZATION`, and the four collection processes | 4 | 16 GB | 2 h |

Independent stages overlap inside the 32-CPU allocation. After assembly,
`MAP_DEPTH_QUICKBIN` (16 CPUs), `SSU_EXTRACT` (8), and `CONTIG_CLASSIFY` (8)
together fill 32 CPUs, and `VIRAL_SCREEN` (12) and the bin screens start as CPUs
free up. The three bin screens request 12, 8, and 20 CPUs, so they cannot all
run at once. `main.nf` submits them in the order GVClass, CheckEUK, CheckM2.
`ASSEMBLE` and `GTDBTK_CLASSIFY` each request all 32 CPUs and 256 GB, so each
runs alone. The actual overlap depends on how long each task runs; read it from
`timeline.html` and `trace.tsv` in the execution directory.

Full mode completes the prokaryotic merge process before it starts optional
gene calling for the same sample. The native outputs of CheckM1, GTDB-Tk, and
Symclatron remain available if gene calling later fails. `EUKARYOTE_QC` does not
gate gene calling. The dependency is per sample; a failure elsewhere in a
multi-sample Nextflow run can still stop the workflow.

## Timeout retries

`READ_QC` and the nested SSUextract `BLAST_ANNOTATE` process get one retry when
Nextflow reports a process running-time timeout. The first `READ_QC` attempt
requests six hours and the second requests 12 hours. The first `BLAST_ANNOTATE`
attempt requests eight hours and the second requests 16 hours. Other `READ_QC`
errors terminate the workflow. Nested BLAST annotation retains its existing
retries for exit code 104 and signal-style exit codes 130 through 145; other
errors use its existing `finish` strategy.

The nested workflow receives `--max_time` from the outer `SSU_EXTRACT` process.
This value caps the request for each nested process attempt. The outer process
limit must cover the eight-hour first attempt, the 16-hour retry, and all other
nested stages. The `ssu` label currently gives `SSU_EXTRACT` 48 hours. The
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
reference. Full mode also needs the environments that the pipeline owns, the
three frozen application clones, and the reference resources of the enabled
stages.

The pipeline owns these Pixi workspaces below `envs/`:

| Workspace | Environments | Registry entry |
|---|---|---|
| `envs/qc` | `checkm1`, `checkm2`, `gtdbtk`, `symclatron` | `checkm_manifest` |
| `envs/viral` | `genomad`, `checkv` | `viral_manifest`, `checkv_manifest` |
| `envs/contig` | `tiara` | `contig_manifest` |
| `envs/eukqc` | `busco`, `tidk` | `eukqc_manifest` |

CheckEUK, GVClass, and SSUextract are frozen Git clones below
`../protist-meta-runtimes/`, next to the repository clone:
`checkeuk-0.13.0-86ab78a`, `gvclass-2.1.1-1517cc0`, and
`ssuextract-1.2.2-4dcf85d`. Each registry entry carries a `git_commit`.
`protist-meta prepare` and the launcher fail unless the clone is at that commit
with a clean working tree.

Install the four workspaces once, then each clone's environment:

```bash
for name in qc viral contig eukqc; do
  pixi install --all --manifest-path "$project_root/envs/$name/pixi.toml" || exit 1
done
for clone in checkeuk-0.13.0-86ab78a gvclass-2.1.1-1517cc0 ssuextract-1.2.2-4dcf85d; do
  (cd "$project_root/../protist-meta-runtimes/$clone" && pixi install --frozen) || exit 1
done
```

After installing `envs/qc`, run `pixi run --manifest-path envs/qc/pixi.toml --environment symclatron symclatron setup` once to download the Symclatron data bundle (1.1 GB) into that environment; preflight checks for it and fails early when it is missing.

**Do not install, reinstall, or update any environment while a job uses it.**
A running job crashed with a bus error on 2026-11-10 after an environment
changed under it. Wait until every run that uses the environment has finished,
then reinstall. This includes the repository setup jobs below, which reinstall
the BBTools, gene-calling, and annotation environments.

Install the reference resources of the enabled stages separately. The repository
scripts do not provision them.

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
`scripts/install-dfam.slurm` next. `install-tools.slurm` binds this directory
into the `repeatmasker` environment, so it needs the directory to exist. A run
with the default options does not use Dfam; only gene calling with repeat
masking does. Skip this installer when that release is already present and valid
because it refuses to replace an existing release. See
[Install Dfam](install-dfam.md) for its checks.

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

The setup commands above also install the gene-calling and annotation stages.
Gene calling is off by default and repeat masking is opt-in. For an installation that omits them, use the
[dependency contract](../reference/inputs.md#dependency-registry) to identify the
required environments and resources. The setup scripts do not accept the
runtime skip switches.

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

   For core mode, replace `full` with `core` and request at least 280 GB.
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

## Published results and the work directory

Every process publishes with `publishDir mode: 'link'`. The files below
`results/` are hard links to files in `work/`. Deleting `work/` therefore does
not free the space of the published files, and `results/` and `work/` must stay
on one filesystem. The launcher places both below the run directory. Do not edit
a published file in place: the edit also changes the cached task output in
`work/`. Copy the file first.

## Set the ONT minimum read length

Append `--ont-min-length 1000` to the allocation launcher arguments to retain
ONT reads of at least 1,000 bases. Chopper applies the existing Phred quality
threshold of 10. Without this option, `nextflow.config` supplies the length
default of 500 bases. The option does not change PacBio HiFi filtering.

The launcher always records the effective length in
`prepared/resume-identity.tsv`, including the default of 500. On `resume` or
`resume-reviewed` it compares the saved value with the requested value, or with
the default when the option is omitted, and rejects a difference before starting
Nextflow. Repeat the same option, or omit it when the saved value is the
default. Use a new run directory to change the filter. A run prepared before the
identity recorded this key keeps working when the option is omitted.

## Choose gene calling and annotation

Gene calling is off by default. In full mode the default run omits repeat
masking, RNA alignment, BRAKER3, dedicated Prodigal-GV gene calling, and
functional annotation. The rest of full mode, including bin quality assessment
and classification, still runs.

Append `--enable-gene-calling` to the launcher arguments to switch gene calling
and functional annotation on. Add `--skip-annotation` to retain gene models
without functional annotation. `--skip-gene-calling` selects the default
explicitly. Supply at most one of `--enable-gene-calling` and
`--skip-gene-calling`. Functional annotation needs gene calling. Keep the same
allocation minimums for every combination.

For example, after preparing the sample sheet and checking the queue, submit a
fresh run with gene calling and annotation:

```bash
run_dir="$project_root/results/my-sample-genes"
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  --cpus-per-task=32 --mem=384G --time=72:00:00 \
  scripts/run-allocation.slurm "$samples" "$run_dir" full \
  --enable-gene-calling) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/allocation-${job_id}.receipt"
```

The catalog, executed notebook, and HTML report are produced in every case.
Their stage tables identify the disabled analyses with the reason `Gene calling
not requested (skip_gene_calling=true).`; no gene models or functional
annotation results are reported for those stages.

## Resume the same run

Retain the run's `launch/`, `work/`, and `prepared/` directories. Keep the
source, sample sheet, registry, and referenced Pixi manifests and locks
unchanged. The launcher rejects an identity mismatch. Changed inputs require a
new run directory. A reviewed source repair can use the separate procedure below.

Resume keeps the saved gene-calling choice. The launcher reads it from
`prepared/inputs.json` unless the command line supplies `--enable-gene-calling`
or `--skip-gene-calling`, and it rejects an explicit option that differs from
the saved choice. A run prepared before this option existed counts as gene
calling enabled. Repeat `--skip-annotation` after `resume` when the fresh run
used it: a run started with `--enable-gene-calling --skip-annotation` requires
`full resume --skip-annotation`.

Repeat masking is off by default. Opt in with `--enable-repeat-masking` on
preparation. Resumes retain the prepared repeat choice unless an explicit repeat
switch is supplied. Runs prepared before this option existed count as masking
enabled. An ordinary resume rejects an explicit change to the saved choice.
The repeat-choice restriction applies only in full mode with gene calling
enabled; core mode and runs with gene calling off exclude repeat computation.
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

Keep the original annotation skip switch after the two source digests. The
saved gene-calling choice carries over as in a plain `resume`. To explicitly disable repeat masking for a run prepared with masking
enabled, add `--skip-repeat-masking` to `resume-reviewed`. This is the only
stage-choice change allowed during reviewed recovery. The launcher keeps every
prepared file unchanged, writes the effective repeat choice to the execution
configuration, and records the prepared and effective choices in `recovery.json`.
Repeat this explicit switch on each reviewed recovery attempt. Do not use
`--enable-repeat-masking` and `--skip-repeat-masking` together.

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

After installing the repository environments and the SSUextract clone, submit
the standalone retry check in a small allocation. It tests timeout retries and
ordinary exit-code failures for `READ_QC` (6 h, then 12 h) with the repository's
Nextflow engine and for `BLAST_ANNOTATE` (8 h, then 16 h) with the SSUextract
clone's engine. Run the 32-CPU MICRO workflow gate in a separate allocation.

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
