# microeuk-nf

microeuk-nf analyzes ONT or PacBio HiFi metagenomes that contain microbial
eukaryotes, bacteria, archaea, and viruses. It publishes independent quality and
taxonomy evidence, a validated SQLite catalog, an executed notebook, and an HTML
report.

<figure markdown="span">
  ![Workflow schematic from long-read input through assembly, bin characterization, evidence-based routing, gene calling, annotation, and catalog and report publication](assets/workflow.svg){ width="1200" }
  <figcaption>Candidate routes guide downstream analysis; they are not final taxonomy or accepted MAG labels.</figcaption>
</figure>

The current ONT acceptance case completed 17 tasks and cached all 17 on resume.
Its catalog and report audits matched 18 record counts and 16 report tables, and
all 1,448 artifact hashes passed. Positive prokaryotic and RNA-supported
eukaryotic cases cover additional branches at their recorded source revisions. The
[source-scoped validation summary](https://github.com/nelli-team/microeuk-nf/blob/main/SUMMARY.md)
records the exact source revisions, counts, and interpretation limits.

- [Scientific design](explanation/design.md)
- [Run and resume on Dori](how-to/run-on-dori.md)
- [Input and dependency contract](reference/inputs.md)
- [Report API](reference/report.md)

Generated run outputs and full validation receipts are local artifacts excluded
from Git. This documentation is built locally; no hosted documentation
deployment is claimed.
