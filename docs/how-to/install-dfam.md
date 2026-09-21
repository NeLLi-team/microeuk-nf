# Install Dfam 4.0 for repeat classification

The full workflow requires the FamDB root and curated-consensus partition from
Dfam 4.0. `scripts/install-dfam.py` uses fixed release URLs, checks the
publisher MD5 and byte count for each compressed file, and records SHA-256
checksums for the installed files.

Complete `scripts/bootstrap.slurm` and the [shared Dori queue
check](run-on-dori.md#prepare-environments-and-references) first. From the
repository root, submit the Dfam job on Dori:

```bash
project_root=$PWD
mkdir -p "$project_root/tasks/logs"
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  scripts/install-dfam.slurm) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/dfam-${job_id}.receipt"
```

The Slurm script installs the locked `repeatmasker` Pixi environment before it
installs the database. The installer refuses to replace an existing Dfam release
directory. If the registered release is already present and valid, skip this
job.

The validated database and checksum receipts are written to the `dfam_db` path
in `conf/databases.yaml`. The script then binds that directory in the installed
FamDB `famdb.conf`.

Run `scripts/install-tools.slurm` after the Dfam job succeeds and after every
later Pixi reinstall of the `repeatmasker` environment. Reinstalling the
environment recreates `famdb.conf`, so the binding must be written again.

```bash
receipt=$(sbatch --parsable -M perceus-00 -A grp-org-sc-mgs \
  -p dori --qos=jgi_normal --chdir="$project_root" \
  scripts/install-tools.slurm) || exit 1
job_id=${receipt%%;*}
[[ $job_id =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$receipt" > "$project_root/tasks/tools-${job_id}.receipt"
```

After the job completes, check the binding and release metadata with the
configured Dfam path:

```bash
dfam_db=/absolute/path/from/conf/databases.yaml
bash scripts/configure-famdb.sh check envs/gene/pixi.toml "$dfam_db"
PIXI_NO_INSTALL=true PIXI_FROZEN=true \
pixi run --as-is --quiet --manifest-path envs/gene/pixi.toml \
  --environment repeatmasker famdb.py info
```

The metadata must report database `Dfam`, version `4.0`, and FamDB creation
format `3.0.0`. Full-workflow preflight also exports `repeat_peps` and
`fasta_all` and requires at least one record from each command.
