# Protist metagenome workflow

This project builds a Nextflow workflow for mixed protist, bacterial, archaeal,
and viral communities sequenced with ONT or PacBio HiFi. Its outputs include
independent quality and taxonomy evidence, a validated SQLite catalog, and a
notebook with saved outputs and an HTML export.

The core, real HiFi, full negative control, positive prokaryotic, positive
RNA-supported eukaryotic, and repaired ONT cases passed their recorded
validation checks. The final ONT run completed all 17 processes and cached all
17 on resume. SQLite integrity, foreign keys, 1,448 artifact hashes, native
quality evidence, the executed notebook, and the HTML report passed their
comparison checks. The repository's `SUMMARY.md` records the source boundary
and scientific limits for each case.

- [Design and scientific scope](explanation/design.md)
- [Run and resume on Dori](how-to/run-on-dori.md)
- [Input and dependency reference](reference/inputs.md)
- [Report reference](reference/report.md)

Wiki builds do not execute notebooks.
