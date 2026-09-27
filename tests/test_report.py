from __future__ import annotations

import hashlib
import json
import re
import sqlite3
import sys
from pathlib import Path

import nbformat
import pytest
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
from nbclient import NotebookClient
from nbformat.notebooknode import NotebookNode

from protist_meta.catalog import build_catalog
from protist_meta.report import build_report

FIXTURE = Path(__file__).parent / "fixtures" / "catalog.json"
EXPECTED_COUNTS = {
    "samples": 1,
    "runs": 1,
    "artifacts": 14,
    "stages": 12,
    "stage_inputs": 13,
    "stage_expected_keys": 12,
    "read_stats": 1,
    "assemblies": 1,
    "contigs": 1,
    "bins": 1,
    "memberships": 1,
    "qc": 1,
    "taxonomy": 1,
    "ssu": 1,
    "viruses": 1,
    "phenotypes": 1,
    "genes": 1,
    "annotations": 1,
}


def test_build_report_executes_fixture_and_exports_html(tmp_path: Path) -> None:
    database = build_catalog(FIXTURE, tmp_path / "catalog.sqlite")
    database_digest = _sha256(database)

    with sqlite3.connect(database) as connection:
        ssu_evidence = connection.execute(
            """
            SELECT reference_identifiers, reference_source, reference_versions,
                   reference_taxonomy, taxonomy_assignment_method, taxon_name
            FROM ssu WHERE record_id = 'ssu1'
            """
        ).fetchone()

    receipt_path = _write_receipt(tmp_path / "execution.json")
    receipt_bytes = receipt_path.read_bytes()
    receipt_digest = _sha256(receipt_path)
    notebook_path = build_report(
        database, tmp_path / "report", execution_provenance=receipt_path
    )

    notebook = _read_notebook(notebook_path)
    html = (notebook_path.parent / "index.html").read_text(encoding="utf-8")
    assert notebook_path.name == "report.executed.ipynb"
    assert notebook_path.stat().st_size > 0
    assert len(html) > 10_000
    assert _error_outputs(notebook) == []
    assert _record_counts(notebook) == EXPECTED_COUNTS
    assert len(_png_outputs(notebook)) == 3
    assert "fixture-revision" in html
    assert (notebook_path.parent / "recovery.json").read_bytes() == receipt_bytes
    notebook_html = "\n".join(_html_outputs(notebook))
    for rendered in (html, notebook_html):
        assert "Original preparation source: fixture-revision" in rendered
        assert "Resumed execution source: reviewed-revision" in rendered
        assert f"Recovery receipt SHA256: {receipt_digest}" in rendered
        assert "&lt;/pre&gt;&lt;script&gt;" in rendered
        assert "</pre><script>" not in rendered
    assert "Bacillota" in html
    assert ssu_evidence == (
        "SILVA:REF_000001.1",
        "SILVA",
        "SILVA:138.2",
        (
            '[{"assignment_method":"native","taxonomy":"Bacteria;'
            'Candidateobacteria","taxonomy_source":"SILVA"}]'
        ),
        "no_runtime_calibration",
        "Unclassified",
    )
    assert "SSU loci: final assignments and raw reference candidates" in html
    assert "SILVA:REF_000001.1" in html
    assert "<td>SILVA</td>" in html
    assert "SILVA:138.2" in html
    assert "Candidateobacteria" in html
    assert "no_runtime_calibration" in html
    assert "does not establish annotation completeness" in html
    assert "exploratory" in html
    assert "per-rank probability" in html
    for heading in (
        "completeness_query_support",
        "completeness_model_reliability",
        "contamination_query_support",
        "contamination_model_reliability",
        "lineage_confidence",
    ):
        assert heading in html
    assert (
        re.search(r'<(?:script|link)[^>]+(?:src|href)=["\x27]https?://', html) is None
    )
    assert _sha256(database) == database_digest


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("original_source_revision", "wrong", "original_source_revision mismatch"),
        ("sample_id", "wrong", "sample_id/run_id absent"),
        ("run_id", "wrong", "sample_id/run_id absent"),
        ("schema_version", True, "schema_version must be 1"),
        ("preparation_identity", [], "preparation_identity must be an object"),
        ("execution_directory", "relative", "execution_directory must be absolute"),
    ],
)
def test_report_rejects_invalid_provenance(
    tmp_path: Path, field: str, value: object, message: str
) -> None:
    database = build_catalog(FIXTURE, tmp_path / "catalog.sqlite")
    original_database = database.read_bytes()
    receipt_path = _write_receipt(tmp_path / "execution.json")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt[field] = value
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")

    with pytest.raises((ValueError, TypeError), match=message):
        build_report(database, tmp_path / "report", execution_provenance=receipt_path)

    assert not (tmp_path / "report").exists()
    assert database.read_bytes() == original_database


def test_relocated_report_reruns_with_sibling_catalog(tmp_path: Path) -> None:
    source_bundle = tmp_path / "source-bundle"
    database = build_catalog(
        FIXTURE,
        source_bundle / "catalog" / "protist-meta.sqlite",
    )
    source_database = database.resolve()
    source_digest = _sha256(source_database)
    report_work = tmp_path / "report-work"
    report_work.mkdir()
    (report_work / "catalog").symlink_to(database.parent, target_is_directory=True)
    staged_database = report_work / "catalog" / database.name
    assert staged_database.resolve() == source_database
    notebook_path = build_report(
        staged_database,
        report_work / "report",
    )
    notebook = _read_notebook(notebook_path)
    setup = notebook.cells[1].source
    assert not (notebook_path.parent / "recovery.json").exists()
    assert "Reviewed-source recovery" not in "\n".join(_html_outputs(notebook))
    assert 'DATABASE_PATH = Path("../catalog/protist-meta.sqlite")' in setup
    assert str(source_database) in setup
    assert source_digest in setup

    notebook_path.parent.rename(source_bundle / "report")
    relocated_bundle = tmp_path / "relocated-bundle"
    source_bundle.rename(relocated_bundle)
    assert not source_database.exists()
    relocated_notebook_path = relocated_bundle / "report" / notebook_path.name
    relocated_notebook = _read_notebook(relocated_notebook_path)

    kernels_dir = Path(sys.prefix) / "share" / "jupyter" / "kernels"
    spec_manager = KernelSpecManager(kernel_dirs=[str(kernels_dir)])
    kernel_spec = spec_manager.get_kernel_spec("python3")
    assert Path(kernel_spec.argv[0]).resolve() == Path(sys.executable).resolve()
    kernel_manager = KernelManager(
        kernel_name="python3",
        kernel_spec_manager=spec_manager,
    )
    client = NotebookClient(
        relocated_notebook,
        km=kernel_manager,
        kernel_name="python3",
        timeout=600,
        allow_errors=False,
    )
    client.execute(cwd=str(relocated_notebook_path.parent))

    assert _error_outputs(relocated_notebook) == []
    rendered_html = "\n".join(_html_outputs(relocated_notebook))
    assert "../catalog/protist-meta.sqlite" in rendered_html
    assert str(source_database) in rendered_html
    assert source_digest in rendered_html


def _write_receipt(path: Path) -> Path:
    receipt = {
        "schema_version": 1,
        "recovery_mode": "reviewed-source",
        "sample_id": "S1",
        "run_id": "R1",
        "original_source_revision": "fixture-revision",
        "execution_source_revision": "reviewed-revision",
        "execution_directory": "/reviewed/source",
        "started_at": "2026-09-27T00:00:00+00:00",
        "forced_processes": [
            "ROUTE_BINS",
            "COLLECT_RECORDS",
            "BUILD_CATALOG",
            "BUILD_REPORT",
        ],
        "preparation_identity": {
            "input": "</pre><script>alert('unsafe')</script>\n\"'\\",
            "checked": True,
        },
        "preparation_identity_sha256": "b" * 64,
        "allocation_config_sha256": "c" * 64,
    }
    path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return path


def _read_notebook(path: Path) -> NotebookNode:
    with path.open(encoding="utf-8") as handle:
        return nbformat.read(handle, as_version=4)


def _outputs(notebook: NotebookNode) -> list[NotebookNode]:
    return [
        output
        for cell in notebook.cells
        if cell.cell_type == "code"
        for output in cell.get("outputs", [])
    ]


def _error_outputs(notebook: NotebookNode) -> list[NotebookNode]:
    return [output for output in _outputs(notebook) if output.output_type == "error"]


def _png_outputs(notebook: NotebookNode) -> list[NotebookNode]:
    return [
        output for output in _outputs(notebook) if "image/png" in output.get("data", {})
    ]


def _html_outputs(notebook: NotebookNode) -> list[str]:
    return [
        output["data"]["text/html"]
        for output in _outputs(notebook)
        if "text/html" in output.get("data", {})
    ]


def _record_counts(notebook: NotebookNode) -> dict[str, int]:
    bundles = [
        output["data"]["application/json"]["record_counts"]
        for output in _outputs(notebook)
        if "application/json" in output.get("data", {})
    ]
    assert len(bundles) == 1
    return bundles[0]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
