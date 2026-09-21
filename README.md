# protist-meta-nf

Nextflow workflow under development for ONT and PacBio HiFi protist enrichment
metagenomes. The requested analysis includes assembly, binning, domain-specific
quality assessment, taxonomy, RNA-supported gene calling, functional annotation,
and a validated SQLite catalog with an executed notebook and HTML report.

Pixi manages Nextflow, Java, and the scientific tool environments.
Reference databases are declared in `conf/databases.yaml`; sample paths are
resolved relative to the input TSV. The workflow does not modify raw reads.

See [validation status](SUMMARY.md) before running an analysis,
[run on Dori](docs/how-to/run-on-dori.md) for submission and resume commands,
and [design](docs/explanation/design.md) for the scientific contract.

- `main.nf` and `modules/`: Nextflow process graph and tool invocations.
- `pixi.toml`, `pixi.lock`, and `envs/`: runtime and tool environments.
- `schema/catalog.yaml`: field meanings, units, types, and record relationships.
- `src/protist_meta/`: input validation, result ingestion, SQLite, and reports.
- `data/samples.example.tsv`: input column contract; replace example paths.
- `docs/`: workflow reference and operator documentation.

Run outputs and working notes are excluded from Git. The source repository has
no configured remote or deployment workflow.
