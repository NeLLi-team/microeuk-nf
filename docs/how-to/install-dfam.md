# Install Dfam 4.0 for repeat classification

The full workflow requires the FamDB root and curated-consensus partition from
Dfam 4.0. The installer reads the fixed release URLs in
`scripts/install-dfam.py`. It checks the publisher MD5 and byte count for each
compressed file, then records SHA-256 checksums for the installed files.

First install the root Pixi environment with `scripts/bootstrap.slurm`.
Submit the Dfam installer from the project root:

```bash
mkdir -p tasks/logs
sbatch -M perceus-00 -A grp-org-sc-mgs -p dori --qos=jgi_normal \
  --kill-on-invalid-dep=yes --chdir="$PWD" scripts/install-dfam.slurm
```

The installer refuses to replace an existing release directory. It writes the
validated files and checksum receipts to the `dfam_db` path in
`conf/databases.yaml`, then binds that directory in the installed FamDB
`famdb.conf`.

Run the tool installer after every Pixi reinstall of the `repeatmasker`
environment. Reinstalling the environment recreates `famdb.conf`, so the setup
step must run again.

```bash
sbatch -M perceus-00 -A grp-org-sc-mgs -p dori --qos=jgi_normal \
  --kill-on-invalid-dep=yes --chdir="$PWD" scripts/install-tools.slurm
```

Check the binding and release metadata after the job completes:

```bash
bash scripts/configure-famdb.sh check \
  envs/gene/pixi.toml \
  /clusterfs/jgi/scratch/science/mgs/nelli/databases/dfam/4.0/famdb
pixi run --as-is --quiet --manifest-path envs/gene/pixi.toml \
  --environment repeatmasker famdb.py info
```

The metadata must report database `Dfam`, version `4.0`, and FamDB creation
format `3.0.0`. Full workflow preflight also exports `repeat_peps` and
`fasta_all` and requires at least one record from each command.
