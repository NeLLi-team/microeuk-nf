"""Build an executed, self-contained report from a protist catalog database."""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import sys
from pathlib import Path

import nbformat
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
from nbclient import NotebookClient
from nbconvert import HTMLExporter
from nbformat.notebooknode import NotebookNode

__all__ = ["build_report"]

KERNEL_NAME = "python3"
KERNEL_DISPLAY_NAME = "Python 3 (Pixi)"
EXECUTED_NOTEBOOK_NAME = "report.executed.ipynb"
HTML_REPORT_NAME = "index.html"

INTRO = """
# Protist metagenome catalog report

**Purpose.** Summarize one published SQLite catalog without reading workflow
output directories or artifact files.

**Data source.** The supplied database is opened in SQLite read-only mode. Its
SHA-256 digest, schema and workflow versions, run identities, methods, and
artifact provenance are recorded below.

**Analysis.** The report covers read filtering, assemblies, bin classes,
tool-specific quality estimates, taxonomy evidence and same-rank conflict
candidates, SSU loci, viral calls, phenotypes, protein-coding models, annotation
hits, stage states, and provenance. SQL summaries keep large feature tables out
of the notebook.

Tables report cataloged rows. Figures appear only when source records exist.
Counts establish catalog coverage, not biological completeness. A core-only
fixture or smoke run validates report execution and database wiring only.
""".strip()

SETUP_CODE = '''
import hashlib
import html
import importlib.metadata
import platform
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import HTML, Markdown, display
from matplotlib_inline.backend_inline import set_matplotlib_formats

set_matplotlib_formats("png")

DATABASE_PATH = Path(__DATABASE_LITERAL__)
SOURCE_DATABASE_PATH = Path(__SOURCE_DATABASE_LITERAL__)
SOURCE_DATABASE_SHA256 = __SOURCE_DATABASE_SHA256_LITERAL__


def query(sql, parameters=()):
    return pd.read_sql_query(sql, connection, params=parameters)


def show(title, sql, empty_message="No cataloged records."):
    frame = query(sql)
    display(HTML(f"<h3>{html.escape(title)}</h3>"))
    if frame.empty:
        message = html.escape(empty_message)
        display(HTML(f"<p><em>{message}</em></p>"))
        return frame
    table = frame.to_html(index=False, escape=True, na_rep="NA")
    display(
        HTML(
            '<div style="overflow:auto;max-height:32rem;'
            'border:1px solid #d0d7de">' + table + "</div>"
        )
    )
    return frame


def database_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if database_sha256(DATABASE_PATH) != SOURCE_DATABASE_SHA256:
    raise ValueError("delivered catalog SHA256 does not match the report source")
DATABASE_URI = DATABASE_PATH.resolve(strict=True).as_uri() + "?mode=ro"
connection = sqlite3.connect(DATABASE_URI, uri=True)
connection.execute("PRAGMA query_only = ON")


metadata = query(
    """
    SELECT schema_version, workflow_name, workflow_version, source_revision
    FROM catalog_metadata
    """
)
if len(metadata) != 1:
    raise ValueError("catalog_metadata must contain exactly one provenance row")

execution_utc = datetime.now(timezone.utc).isoformat()
provenance = pd.DataFrame(
    [
        ("database_path", str(DATABASE_PATH)),
        ("source_database_path", str(SOURCE_DATABASE_PATH)),
        ("database_sha256", SOURCE_DATABASE_SHA256),
        ("schema_version", metadata.at[0, "schema_version"]),
        ("workflow_name", metadata.at[0, "workflow_name"]),
        ("workflow_version", metadata.at[0, "workflow_version"]),
        ("source_revision", metadata.at[0, "source_revision"]),
        ("execution_utc", execution_utc),
        ("kernel_name", "python3"),
        ("python", platform.python_version()),
        ("python_executable", sys.executable),
        ("environment_prefix", sys.prefix),
        ("sqlite", sqlite3.sqlite_version),
        ("pandas", importlib.metadata.version("pandas")),
        ("matplotlib", importlib.metadata.version("matplotlib")),
        ("nbclient", importlib.metadata.version("nbclient")),
    ],
    columns=["field", "value"],
)
display(HTML(provenance.to_html(index=False, escape=True)))

report_tables = [
    "samples", "runs", "artifacts", "stages", "stage_inputs",
    "stage_expected_keys", "read_stats", "assemblies", "contigs", "bins",
    "memberships", "qc", "taxonomy", "ssu", "viruses", "phenotypes",
    "genes", "annotations",
]
record_counts = {
    table: int(query(f'SELECT COUNT(*) AS n FROM "{table}"').at[0, "n"])
    for table in report_tables
}
counts = pd.DataFrame(record_counts.items(), columns=["table", "records"])
display(
    {
        "text/html": counts.to_html(index=False, escape=True),
        "application/json": {"record_counts": record_counts},
    },
    raw=True,
)

scope_text = " ".join(metadata.astype(str).iloc[0]).lower()
if "fixture" in scope_text or "smoke" in scope_text:
    display(
        Markdown(
            "**Scope: fixture/smoke report.** These rows exercise catalog and report "
            "wiring. They do not validate scientific branches or biological "
            "completeness."
        )
    )
'''.strip()

READ_ASSEMBLY = """
## Reads and assemblies

Read counts and base counts are separate. Assembly length and N50 are in base
pairs; GC is a fraction from zero to one. `NA` means the producer did not report
a value, not zero.
""".strip()

READ_ASSEMBLY_CODE = '''
runs = show(
    "Runs",
    "SELECT sample_id, run_id, platform, rna_technology, basecaller, "
    "basecaller_model, library_prep "
    "FROM runs ORDER BY sample_id, run_id",
    "No workflow runs are cataloged.",
)
read_stats = show(
    "Read filtering",
    """
    SELECT sample_id, run_id, record_id, total_reads, total_bases,
           retained_reads, retained_bases, rejected_reads, rejected_bases
    FROM read_stats ORDER BY sample_id, run_id, record_id
    """,
    "No read statistics are cataloged.",
)
assemblies = show(
    "Assembly summaries",
    """
    SELECT sample_id, run_id, assembly_id, assembly_length_bp, contig_count,
           n50_bp, gc_fraction
    FROM assemblies ORDER BY sample_id, run_id, assembly_id
    """,
    "No assembly summaries are cataloged.",
)
'''.strip()

READ_FIGURE_CODE = """
if read_stats.empty:
    display(Markdown("*No read records are available for this figure.*"))
else:
    read_plot = (
        read_stats.groupby(["sample_id", "run_id"], as_index=False)[
            ["retained_reads", "rejected_reads"]
        ]
        .sum()
        .assign(total=lambda x: x["retained_reads"] + x["rejected_reads"])
        .nlargest(30, "total")
        .sort_values(["sample_id", "run_id"])
    )
    labels = read_plot["sample_id"] + "\\n" + read_plot["run_id"]
    height = max(3, 0.42 * len(labels) + 1.2)
    fig, ax = plt.subplots(figsize=(8, height))
    ax.barh(labels, read_plot["retained_reads"], color="#0072B2", label="Retained")
    ax.barh(
        labels, read_plot["rejected_reads"],
        left=read_plot["retained_reads"], color="#E69F00", label="Rejected",
    )
    ax.set_xlabel("Reads")
    ax.set_ylabel("Sample / run")
    ax.invert_yaxis()
    ax.legend(frameon=False, ncol=2, loc="lower right", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    display(fig)
    plt.close(fig)
""".strip()

READ_CAPTION = """
**Figure 1. Read outcomes.** Retained and rejected reads are summed by sample
and run. At most the 30 runs with the most classified reads are plotted; the
table remains complete.
""".strip()

ASSEMBLY_FIGURE_CODE = """
if assemblies.empty:
    display(Markdown("*No assembly records are available for this figure.*"))
else:
    assembly_plot = (
        assemblies.groupby(["sample_id", "run_id"], as_index=False)[
            "assembly_length_bp"
        ]
        .sum().nlargest(30, "assembly_length_bp")
        .sort_values(["sample_id", "run_id"])
    )
    labels = assembly_plot["sample_id"] + "\\n" + assembly_plot["run_id"]
    height = max(3, 0.42 * len(labels) + 1.2)
    fig, ax = plt.subplots(figsize=(8, height))
    ax.barh(labels, assembly_plot["assembly_length_bp"] / 1_000_000,
            color="#0072B2")
    ax.set_xlabel("Assembly sequence (Mbp)")
    ax.set_ylabel("Sample / run")
    ax.invert_yaxis()
    fig.tight_layout()
    display(fig)
    plt.close(fig)
""".strip()

ASSEMBLY_CAPTION = """
**Figure 2. Cataloged assembly sequence.** Assembled base pairs are summed by
sample and run and displayed in megabase pairs. At most the 30 longest
sample/run assemblies are plotted.
""".strip()

BIN_QC = """
## Bins and quality evidence

Candidate classes route later analysis; they are not accepted-MAG calls. Each
QC row retains its tool. Generic completeness and contamination are percentages
with tool-specific definitions. CheckEUK percentages, marker counts out of 74
or the reported denominator, and distinct families out of 1,085 stay separate.
CheckM1 and CheckM2 provide prokaryotic quality evidence, including when they run
on eukaryotic candidates. A successful run does not validate the scores outside
the prokaryotic domain. For GVClass, query support describes the current row,
while model reliability describes validation of the model. Estimates labeled
`unvalidated_legacy` are exploratory and do not certify the estimate.
""".strip()

BIN_QC_CODE = '''
bins = show(
    "Bins",
    """
    SELECT sample_id, run_id, bin_id, assembly_id, candidate_class
    FROM bins ORDER BY sample_id, run_id, bin_id
    """,
    "No bins are cataloged.",
)
qc = show(
    "Tool-specific quality assessments",
    """
    SELECT q.sample_id, q.run_id, q.record_id, q.target_kind, q.target_id,
           s.tool_name, s.tool_version, q.completeness_percent,
           q.contamination_percent, q.completeness_basis,
           q.contamination_basis,
           q.completeness_query_support, q.completeness_model_reliability,
           q.contamination_query_support, q.contamination_model_reliability,
           q.checkeuk_pfam_ridge_completeness_percent,
           q.checkeuk_bom_completeness_percent, q.checkeuk_inventory74_present,
           q.checkeuk_expected_markers_present,
           q.checkeuk_expected_markers_total, q.checkeuk_pfam1085_distinct
    FROM qc q JOIN stages s USING (sample_id, run_id, stage_id)
    ORDER BY q.sample_id, q.run_id, q.target_kind, q.target_id, s.tool_name
    """,
    "No quality assessments are cataloged.",
)
'''.strip()

BIN_FIGURE_CODE = """
if bins.empty:
    display(Markdown("*No bin records are available for this figure.*"))
else:
    class_counts = bins.groupby("candidate_class", as_index=False).size()
    class_counts = class_counts.sort_values("size", ascending=False)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(class_counts["candidate_class"], class_counts["size"], color="#009E73")
    ax.set_ylabel("Bins")
    ax.tick_params(axis="x", labelrotation=45)
    fig.tight_layout()
    display(fig)
    plt.close(fig)
""".strip()

BIN_CAPTION = """
**Figure 3. Evidence-routing classes.** Counts show every cataloged bin in its
recorded candidate class. They do not state genome quality or MAG acceptance.
""".strip()

EVIDENCE = """
## Taxonomy, SSU, viral, and phenotype evidence

Taxonomy assertions remain tool-specific. Conflict candidates have more than
one label for the same target and rank; this does not resolve synonyms or
reference-release differences. For GVClass, `lineage_confidence` qualifies the
full native lineage and is not a per-rank probability. SSU and viral coordinates
are one-based and inclusive. SSU final assignments remain separate from raw
candidate reference evidence and its assignment method. Phenotypes retain class,
score, threshold, and applicability.
""".strip()

EVIDENCE_CODE = '''
taxonomy = show(
    "All taxonomy assertions",
    """
    SELECT t.sample_id, t.run_id, t.record_id, t.target_kind, t.target_id,
           s.tool_name, s.tool_version, s.database_name, s.database_version,
           t.rank, t.taxon_name, t.taxon_id, t.score, t.lineage_confidence
    FROM taxonomy t JOIN stages s USING (sample_id, run_id, stage_id)
    ORDER BY t.sample_id, t.run_id, t.target_kind, t.target_id, t.rank
    """,
    "No taxonomy evidence is cataloged.",
)
taxonomy_conflicts = show(
    "Same-rank taxonomy conflict candidates",
    """
    SELECT t.sample_id, t.run_id, t.target_kind, t.target_id, t.rank,
           COUNT(DISTINCT t.taxon_name) AS taxon_labels,
           GROUP_CONCAT(DISTINCT s.tool_name || ': ' || t.taxon_name) AS evidence
    FROM taxonomy t JOIN stages s USING (sample_id, run_id, stage_id)
    GROUP BY t.sample_id, t.run_id, t.target_kind, t.target_id, t.rank
    HAVING COUNT(DISTINCT t.taxon_name) > 1
    ORDER BY t.sample_id, t.run_id, t.target_kind, t.target_id, t.rank
    """,
    "No same-rank label conflicts were found.",
)
ssu = show(
    "SSU loci: final assignments and raw reference candidates",
    """
    SELECT sample_id, run_id, record_id, contig_id, ssu_type, start, end,
           strand, hit_accession, taxon_name AS final_taxon_name,
           taxonomy_assignment_method,
           reference_identifiers, reference_source, reference_versions,
           reference_taxonomy,
           identity_percent, query_coverage_percent
    FROM ssu ORDER BY sample_id, run_id, contig_id, start, record_id
    """,
    "No SSU loci are cataloged.",
)
viruses = show(
    "Viral contig calls",
    """
    SELECT v.sample_id, v.run_id, v.record_id, v.contig_id, v.start, v.end,
           v.strand, s.tool_name, s.tool_version, v.score, v.viral_taxon
    FROM viruses v JOIN stages s USING (sample_id, run_id, stage_id)
    ORDER BY v.sample_id, v.run_id, v.contig_id, v.start, s.tool_name
    """,
    "No viral calls are cataloged.",
)
phenotypes = show(
    "Phenotype predictions",
    """
    SELECT p.sample_id, p.run_id, p.record_id, p.bin_id, s.tool_name,
           s.tool_version, p.raw_class, p.raw_score, p.threshold,
           p.passes_threshold, p.thresholded_class, p.is_applicable,
           p.applicability_reason
    FROM phenotypes p JOIN stages s USING (sample_id, run_id, stage_id)
    ORDER BY p.sample_id, p.run_id, p.bin_id, s.tool_name
    """,
    "No phenotype predictions are cataloged.",
)
'''.strip()

PROVENANCE = """
## Coding models, stage states, methods, and artifacts

SQLite supplies counts of coding models and annotation hits. Each translated
BRAKER isoform has a separate model record. The report links to complete feature
artifacts. Each stage records its status, reason, command, software and database
versions, input count, and expected-key count.
""".strip()

PROVENANCE_CODE = '''
gene_counts = show(
    "Protein-coding model counts",
    """
    SELECT g.sample_id, g.run_id, g.bin_id, g.stage_id, s.name AS stage_name,
           s.status, s.status_reason, s.tool_name, s.tool_version,
           COUNT(*) AS coding_models
    FROM genes g JOIN stages s USING (sample_id, run_id, stage_id)
    GROUP BY g.sample_id, g.run_id, g.bin_id, g.stage_id, s.name, s.status,
             s.status_reason, s.tool_name, s.tool_version
    ORDER BY g.sample_id, g.run_id, g.bin_id, g.stage_id
    """,
    "No protein-coding models are cataloged.",
)
annotation_counts = show(
    "Functional-annotation counts",
    """
    SELECT sample_id, run_id, database_name, database_version,
           COUNT(*) AS annotation_hits,
           COUNT(DISTINCT gene_id) AS annotated_models
    FROM annotations
    GROUP BY sample_id, run_id, database_name, database_version
    ORDER BY sample_id, run_id, database_name, database_version
    """,
    "No functional annotations are cataloged.",
)
braker_states = show(
    "BRAKER states",
    """
    SELECT s.sample_id, s.run_id, s.stage_id, s.name, s.status,
           s.status_reason, s.tool_name, s.tool_version,
           COUNT(g.gene_id) AS cataloged_models
    FROM stages s LEFT JOIN genes g USING (sample_id, run_id, stage_id)
    WHERE lower(s.name) LIKE '%braker%' OR lower(s.tool_name) LIKE '%braker%'
    GROUP BY s.sample_id, s.run_id, s.stage_id, s.name, s.status,
             s.status_reason, s.tool_name, s.tool_version
    ORDER BY s.sample_id, s.run_id, s.stage_id
    """,
    "No BRAKER stage is cataloged.",
)
display(
    Markdown(
        "A skipped, pending, failed, unrun, or empty BRAKER branch does not "
        "establish annotation completeness or biological gene absence."
    )
)
feature_artifacts = show(
    "Complete feature artifacts",
    """
    SELECT sample_id, run_id, artifact_id, kind, path, sha256, size_bytes,
           producer_stage_id
    FROM artifacts
    WHERE kind IN ('gene_annotation', 'protein_fasta', 'functional_annotation')
    ORDER BY sample_id, run_id, kind, artifact_id
    """,
    "No gene or functional-annotation artifacts are cataloged.",
)
stages = show(
    "Stage states and methods",
    """
    SELECT s.sample_id, s.run_id, s.stage_id, s.name, s.status,
           s.status_reason, s.notes, s.tool_name, s.tool_version, s.command,
           s.database_name, s.database_version, s.database_sha256,
           (SELECT COUNT(DISTINCT i.position)
            FROM stage_inputs i
            WHERE i.sample_id = s.sample_id AND i.run_id = s.run_id
              AND i.stage_id = s.stage_id) AS input_artifacts,
           (SELECT COUNT(DISTINCT e.collection || ':' || e.record_id)
            FROM stage_expected_keys e
            WHERE e.sample_id = s.sample_id AND e.run_id = s.run_id
              AND e.stage_id = s.stage_id) AS expected_records
    FROM stages s
    ORDER BY s.sample_id, s.run_id, s.stage_id
    """,
    "No workflow stages are cataloged.",
)
artifacts = show(
    "Artifact provenance",
    """
    SELECT sample_id, run_id, artifact_id, kind, path, sha256, size_bytes,
           producer_stage_id
    FROM artifacts ORDER BY sample_id, run_id, kind, artifact_id
    """,
    "No artifacts are cataloged.",
)
connection.close()
display(Markdown(f"Report execution completed at `{execution_utc}`."))
'''.strip()

INTERPRETATION = """
## Interpretation boundary

This report does not infer that a skipped, pending, failed, or empty stage
succeeded, and it does not treat a union of classifier labels as an accepted
MAG. Scientific branch coverage requires verified tool outputs and cataloged
provenance.
""".strip()


def build_report(
    database: Path, output_dir: Path, *, execution_provenance: Path | None = None
) -> Path:
    """Execute a read-only catalog report and export its HTML sibling.

    The ``python3`` kernelspec is loaded only from ``sys.prefix``. User-level
    Jupyter state is unchanged.

    Args:
        database: Published protist catalog SQLite file.
        output_dir: Destination for the executed notebook and HTML report.
        execution_provenance: Optional reviewed-source recovery receipt.

    Returns:
        Path to ``report.executed.ipynb``.

    Raises:
        FileNotFoundError: If ``database`` does not exist.
        ValueError: If the database or recovery receipt has invalid provenance.
        TypeError: If a recovery receipt JSON object has an invalid type.
    """
    database_reference = Path(os.path.relpath(database, start=output_dir))
    source_database = database.resolve(strict=True)
    if not source_database.is_file():
        raise ValueError(f"catalog database is not a regular file: {source_database}")
    receipt = (
        _read_execution_provenance(execution_provenance, source_database)
        if execution_provenance is not None
        else None
    )
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    notebook = _build_notebook(
        database_reference,
        source_database,
        _sha256(source_database),
    )
    if receipt is not None:
        notebook.cells.insert(2, nbformat.v4.new_code_cell(_receipt_code(receipt)))
        (output_dir / "recovery.json").write_bytes(receipt)
    _execute_notebook(notebook, output_dir)
    notebook_path = output_dir / EXECUTED_NOTEBOOK_NAME
    nbformat.write(notebook, notebook_path)
    _export_html(notebook, output_dir / HTML_REPORT_NAME)
    return notebook_path


def _read_execution_provenance(path: Path, database: Path) -> bytes:
    receipt = path.read_bytes()
    provenance = json.loads(receipt)
    if not isinstance(provenance, dict):
        raise TypeError("execution provenance must be a JSON object")
    if (
        type(provenance.get("schema_version")) is not int
        or provenance.get("schema_version") != 1
    ):
        raise ValueError("execution provenance schema_version must be 1")
    if provenance.get("recovery_mode") != "reviewed-source":
        raise ValueError("execution provenance recovery_mode must be reviewed-source")
    for field in (
        "sample_id",
        "run_id",
        "original_source_revision",
        "execution_source_revision",
        "execution_directory",
        "started_at",
        "preparation_identity_sha256",
        "allocation_config_sha256",
    ):
        if not isinstance(provenance.get(field), str) or not provenance[field]:
            raise ValueError(f"execution provenance {field} must be a nonempty string")
    if not Path(provenance["execution_directory"]).is_absolute():
        raise ValueError("execution provenance execution_directory must be absolute")
    if not isinstance(provenance.get("preparation_identity"), dict):
        raise TypeError("execution provenance preparation_identity must be an object")
    if provenance.get("forced_processes") != [
        "ROUTE_BINS",
        "COLLECT_RECORDS",
        "BUILD_CATALOG",
        "BUILD_REPORT",
    ]:
        raise ValueError("execution provenance forced_processes do not match recovery")
    _check_catalog_identity(
        database,
        provenance["original_source_revision"],
        provenance["sample_id"],
        provenance["run_id"],
    )
    return receipt


def _check_catalog_identity(
    database: Path, original_source_revision: str, sample_id: str, run_id: str
) -> None:
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as connection:
        sources = connection.execute(
            "SELECT source_revision FROM catalog_metadata"
        ).fetchall()
        if sources != [(original_source_revision,)]:
            raise ValueError("execution provenance original_source_revision mismatch")
        run = connection.execute(
            "SELECT 1 FROM runs WHERE sample_id = ? AND run_id = ?",
            (sample_id, run_id),
        ).fetchone()
        if run is None:
            raise ValueError(
                "execution provenance sample_id/run_id absent from catalog"
            )


def _receipt_code(receipt: bytes) -> str:
    digest = hashlib.sha256(receipt).hexdigest()
    return (
        "import json\n\n"
        f"recovery = json.loads({receipt.decode('utf-8')!r})\n"
        f"recovery_sha256 = {digest!r}\n"
        'display(HTML("<h2>Reviewed-source recovery</h2>"))\n'
        "for label, value in (\n"
        '    ("Original preparation source", recovery["original_source_revision"]),\n'
        '    ("Resumed execution source", recovery["execution_source_revision"]),\n'
        '    ("Recovery receipt SHA256", recovery_sha256),\n'
        "):\n"
        '    display(HTML("<p>" + html.escape(label + ": " + value) + "</p>"))\n'
        'display(HTML("<pre>" + html.escape(json.dumps(recovery, indent=2))'
        ' + "</pre>"))'
    )


def _build_notebook(
    database: Path,
    source_database: Path,
    source_database_sha256: str,
) -> NotebookNode:
    setup = SETUP_CODE
    replacements = {
        "__DATABASE_LITERAL__": json.dumps(str(database)),
        "__SOURCE_DATABASE_LITERAL__": json.dumps(str(source_database)),
        "__SOURCE_DATABASE_SHA256_LITERAL__": json.dumps(source_database_sha256),
    }
    for placeholder, value in replacements.items():
        setup = setup.replace(placeholder, value)
    specs = [
        ("markdown", INTRO),
        ("code", setup),
        ("markdown", READ_ASSEMBLY),
        ("code", READ_ASSEMBLY_CODE),
        ("code", READ_FIGURE_CODE),
        ("markdown", READ_CAPTION),
        ("code", ASSEMBLY_FIGURE_CODE),
        ("markdown", ASSEMBLY_CAPTION),
        ("markdown", BIN_QC),
        ("code", BIN_QC_CODE),
        ("code", BIN_FIGURE_CODE),
        ("markdown", BIN_CAPTION),
        ("markdown", EVIDENCE),
        ("code", EVIDENCE_CODE),
        ("markdown", PROVENANCE),
        ("code", PROVENANCE_CODE),
        ("markdown", INTERPRETATION),
    ]
    cells = [
        nbformat.v4.new_code_cell(source)
        if kind == "code"
        else nbformat.v4.new_markdown_cell(source)
        for kind, source in specs
    ]
    return nbformat.v4.new_notebook(
        cells=cells,
        metadata={
            "kernelspec": {
                "display_name": KERNEL_DISPLAY_NAME,
                "language": "python",
                "name": KERNEL_NAME,
            },
            "language_info": {
                "name": "python",
                "version": f"{sys.version_info.major}.{sys.version_info.minor}",
            },
        },
    )


def _execute_notebook(notebook: NotebookNode, output_dir: Path) -> None:
    kernels_dir = Path(sys.prefix) / "share" / "jupyter" / "kernels"
    spec_manager = KernelSpecManager(kernel_dirs=[str(kernels_dir)])
    kernel_manager = KernelManager(
        kernel_name=KERNEL_NAME,
        kernel_spec_manager=spec_manager,
    )
    client = NotebookClient(
        notebook,
        km=kernel_manager,
        kernel_name=KERNEL_NAME,
        timeout=600,
        allow_errors=False,
    )
    client.execute(cwd=str(output_dir))


def _export_html(notebook: NotebookNode, output_path: Path) -> None:
    exporter = HTMLExporter()
    exporter.template_name = "lab"
    exporter.embed_images = True
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True
    exporter.mathjax_url = ""
    exporter.require_js_url = ""
    body, _ = exporter.from_notebook_node(notebook)
    output_path.write_text(body, encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
