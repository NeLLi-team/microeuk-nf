# Report API

`protist_meta.report.build_report(database: Path, output_dir: Path) -> Path`
creates an executed Jupyter notebook and an HTML report from one published
protist catalog.

## Prerequisites

Complete the repository bootstrap before calling the API. Run the code through
the existing locked root Pixi environment; report execution does not install
packages or create a user-level Jupyter kernel.

The input must be a catalog built by `protist_meta.catalog.build_catalog`. The
caller must reserve the compute allocation used for notebook execution.

## Example

```python
from pathlib import Path

from protist_meta.report import build_report

notebook = build_report(
    Path("results/catalog/protist-meta.sqlite"),
    Path("results/report"),
)
```

The call returns the path to `report.executed.ipynb`. `index.html` is written in
the same output directory.

## Behavior

1. Retains the database path relative to the report directory, resolves the
   original source for provenance, and creates the output directory.
2. Generates a notebook with fixed SQL queries for the catalog schema.
3. Selects the `python3` kernel only from the active Pixi environment.
4. Executes every cell with `nbclient`. A cell error stops the build.
5. Writes the executed notebook and exports HTML with figures embedded.

The function does not install or change a user-level Jupyter kernel.

## Rerun a delivered notebook

Keep the delivered catalog and report in this sibling layout:

```text
bundle/
├── catalog/
│   └── protist-meta.sqlite
└── report/
    ├── report.executed.ipynb
    └── index.html
```

Run nbconvert from the report directory through the project's locked Pixi
environment:

```bash
PROJECT_ROOT=/path/to/microeuk-nf
BUNDLE_ROOT=/path/to/bundle
cd "$BUNDLE_ROOT/report"
PIXI_NO_INSTALL=true PIXI_FROZEN=true \
pixi run --as-is --quiet --manifest-path "$PROJECT_ROOT/pixi.toml" \
  jupyter nbconvert --execute --to notebook --inplace \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=600 \
  report.executed.ipynb
```

The locked environment supplies the `python3` kernelspec, which points to its
Pixi Python interpreter. Running from `report/` preserves
the notebook's `../catalog/protist-meta.sqlite` reference. Execution stops if
the sibling database SHA256 differs from the source digest recorded when the
report was built. The command updates the notebook; it does not replace the
saved HTML export.

## Data access

The notebook opens the supplied SQLite file with `mode=ro` and enables SQLite
query-only mode. It does not scan workflow directories or open artifact paths.
The report records provenance from the database and execution environment:

- SHA-256 digest, relative execution path, and resolved source path of the
  SQLite file
- catalog schema version
- workflow name, workflow version, and source revision
- sample and run identifiers
- UTC execution time
- kernel name, Python executable, environment prefix, and package versions
- stage tool versions, database versions, expected-key counts, and input counts
- artifact paths, sizes, and SHA-256 digests

## Report contents

The notebook contains complete tables for reads, assemblies, bins, quality
assessments, taxonomy assertions, same-rank taxonomy conflict candidates, SSU
loci, viral calls, phenotypes, stage states, methods, and artifacts. Tables are
HTML-escaped and placed in scrollable containers.

Protein-coding models and functional-annotation hits are summarized with SQL
counts. Each translated BRAKER isoform has a separate model record. The report
lists the complete feature artifacts without loading those files into memory.

Up to three figures are included when source records exist:

- retained and rejected reads by sample and run
- cataloged assembly sequence by sample and run
- bin counts by evidence-routing class

Plots with sample/run categories show at most 30 categories. The corresponding
tables remain complete.

## Interpretation

`NA` means that a tool did not report a value. It does not mean zero. Quality
estimates retain their tool identity and unit. CheckEUK completeness estimates
and marker counts remain in separate columns.

Candidate classes route later analysis; they are not accepted-MAG calls. A
same-rank taxonomy conflict is a review candidate, not a resolved taxonomic
decision. A skipped, pending, failed, unrun, or empty BRAKER branch does not
establish annotation completeness or biological gene absence.

A core-only fixture or smoke run proves database and report wiring. It does not
prove that scientific analysis branches ran or that cataloged genomes are
complete.

## Errors

`FileNotFoundError` is raised when the database path does not exist. `ValueError`
is raised when the path is not a regular file, when the delivered catalog
SHA256 differs from the report source, or when `catalog_metadata` does not
contain exactly one row. SQLite query failures and notebook cell errors are
propagated to the caller. Output files are written only after notebook
execution succeeds.
