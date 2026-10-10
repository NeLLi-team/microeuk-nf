"""Portable native bundles preserve counts, policy and source evidence."""

import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

from protist_meta.catalog import build_catalog
from protist_meta.cli import main
from protist_meta.visualization import export_visualization

FIXTURE = Path(__file__).parent / "fixtures" / "catalog.json"


@pytest.fixture
def catalog(tmp_path: Path) -> Path:
    """Create a catalog with digest-bound native mapping and quality evidence."""
    records = tmp_path / "records.json"
    records.write_text(
        FIXTURE.read_text(encoding="utf-8").replace('"filter"', '"read_qc"'),
        encoding="utf-8",
    )
    database = build_catalog(records, tmp_path / "catalog.sqlite")
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text(
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tmeandepth\n"
        "ctg1\t1\t1000\t6\t900\t2.5\n",
        encoding="utf-8",
    )
    report = tmp_path / "checkeuk.tsv"
    report.write_text(
        "genome\tstatus\tlineage\tcontamination_status\n"
        "bin1\tunscorable:no_host_consensus\t\tunassessed\n",
        encoding="utf-8",
    )
    with sqlite3.connect(database) as connection:
        for identifier in ("read-qc-report", "taxonomy-report"):
            table = tmp_path / f"{identifier}.tsv"
            table.write_text("record\tvalue\nfixture\t1\n", encoding="utf-8")
            connection.execute(
                "UPDATE artifacts SET path=?, sha256=?, size_bytes=? "
                "WHERE artifact_id=?",
                (
                    str(table),
                    hashlib.sha256(table.read_bytes()).hexdigest(),
                    table.stat().st_size,
                    identifier,
                ),
            )
        connection.executemany(
            "UPDATE artifacts SET path=?, sha256=?, size_bytes=? WHERE artifact_id=?",
            [
                (
                    str(coverage),
                    hashlib.sha256(coverage.read_bytes()).hexdigest(),
                    coverage.stat().st_size,
                    "coverage-table",
                ),
                (
                    str(report),
                    hashlib.sha256(report.read_bytes()).hexdigest(),
                    report.stat().st_size,
                    "qc-report",
                ),
            ],
        )
        connection.execute(
            "UPDATE stages SET command=? WHERE stage_id='binning'",
            (
                (
                    "minimap2 -ax map-ont --secondary=no; samtools view -F 2308 -q 20; "
                    "samtools coverage -d 0"
                ),
            ),
        )
    return database


def test_cli_exports_native_evidence_and_conserved_read_counts(
    catalog: Path, tmp_path: Path
) -> None:
    output = tmp_path / "visualization"

    assert (
        main(
            [
                "visualization",
                "export",
                "--catalog",
                str(catalog),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    index = json.loads((output / "index.json").read_text(encoding="utf-8"))
    bundle = json.loads((output / index["bundles"][0]["path"]).read_text())

    assert bundle["schema_version"] == 1
    assert bundle["sample_id"] == "S1"
    assert bundle["assembly"]["inventory_status"] == "complete"
    assert bundle["mapping"]["input_reads"] == 8
    assert bundle["mapping"]["read_unit"] == "unpaired reads"
    assert bundle["mapping"]["assembly_mapped_reads"] == 6
    assert bundle["mapping"]["bin_mapped_reads"] == 6
    assert bundle["mapping"]["unbinned_mapped_reads"] == 0
    assert bundle["mapping"]["without_qualifying_alignment_reads"] == 2
    assert bundle["mapping"]["minimum_mapq"] == 20
    assert bundle["mapping"]["excluded_sam_flags"] == 3844
    assert bundle["mapping"]["per_bin"] == [
        {"bin_id": "bin1", "numreads": 6, "covbases": 900}
    ]
    assert bundle["evidence"]["checkeuk"]["rows"][0]["status"] == (
        "unscorable:no_host_consensus"
    )
    assert bundle["evidence"]["gvclass"]["status"] == "unassessed"
    assert bundle["qc"][0]["tool_name"] == "CheckEUK"
    assert "genes" not in bundle
    assert "annotations" not in bundle
    assert export_visualization(catalog, output)[0].name == index["bundles"][0]["path"]


def test_export_rejects_changed_registered_coverage(
    catalog: Path, tmp_path: Path
) -> None:
    (tmp_path / "coverage.tsv").write_text("changed", encoding="utf-8")

    with pytest.raises(ValueError, match="SHA256"):
        export_visualization(catalog, tmp_path / "out")


def test_missing_mapping_preserves_null_counts(catalog: Path, tmp_path: Path) -> None:
    (tmp_path / "coverage.tsv").unlink()

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["mapping"]["status"] == "unavailable"
    assert bundle["mapping"]["input_reads"] == 8
    assert bundle["mapping"]["read_unit"] == "unpaired reads"
    assert bundle["mapping"]["assembly_mapped_reads"] is None
    assert bundle["mapping"]["bin_mapped_reads"] is None


def test_zero_bins_keep_mapped_reads_unbinned(catalog: Path, tmp_path: Path) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute("DELETE FROM memberships")
        connection.execute("DELETE FROM bins")

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["bins"] == []
    assert bundle["mapping"]["bin_mapped_reads"] == 0
    assert bundle["mapping"]["unbinned_mapped_reads"] == 6


def test_export_rejects_reads_above_denominator(catalog: Path, tmp_path: Path) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute("UPDATE read_stats SET retained_reads=5, rejected_reads=5")

    with pytest.raises(ValueError, match="denominator"):
        export_visualization(catalog, tmp_path / "out")


def test_partial_inventory_does_not_claim_complete_assembly_recovery(
    catalog: Path, tmp_path: Path
) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute(
            "UPDATE assemblies SET contig_count=2, assembly_length_bp=2000"
        )

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["assembly"]["inventory_status"] == "partial"
    assert bundle["mapping"]["status"] == "partial"
    assert bundle["mapping"]["without_qualifying_alignment_reads"] is None


def test_unknown_mapping_policy_does_not_invent_primary_counts(
    catalog: Path, tmp_path: Path
) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute(
            "UPDATE stages SET command='unknown' WHERE stage_id='binning'"
        )

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["mapping"]["status"] == "unavailable"
    assert bundle["mapping"]["minimum_mapq"] is None
    assert bundle["mapping"]["assembly_mapped_reads"] is None


def test_run_without_assembly_exports_inventory(catalog: Path, tmp_path: Path) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute("DELETE FROM assemblies")
        connection.execute("DELETE FROM contigs")
        connection.execute("DELETE FROM memberships")
        connection.execute("DELETE FROM bins")

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["assembly"] is None
    assert bundle["read_stats"][0]["retained_reads"] == 8
    assert bundle["mapping"]["status"] == "unavailable"


@pytest.mark.parametrize("artifact_name", ["qc-report", "taxonomy-report"])
def test_completed_stage_requires_registered_evidence(
    catalog: Path, tmp_path: Path, artifact_name: str
) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute(
            "UPDATE artifacts SET path=? WHERE artifact_id=?",
            (str(tmp_path / "missing.tsv"), artifact_name),
        )

    with pytest.raises(FileNotFoundError, match=artifact_name):
        export_visualization(catalog, tmp_path / "out")


def test_denominator_uses_read_qc_when_other_read_statistics_exist(
    catalog: Path, tmp_path: Path
) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute(
            "INSERT INTO read_stats SELECT sample_id, run_id, 'other-stats', "
            "'assemble', 'assembly-fasta', total_reads, total_bases, "
            "retained_reads, retained_bases, rejected_reads, rejected_bases "
            "FROM read_stats"
        )

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["mapping"]["input_reads"] == 8
    assert bundle["mapping"]["denominator_source"] == "reads1"
    assert bundle["mapping"]["status"] == "available"


def test_unavailable_evidence_preserves_failed_stage_status(
    catalog: Path, tmp_path: Path
) -> None:
    (tmp_path / "checkeuk.tsv").unlink()
    with sqlite3.connect(catalog) as connection:
        connection.execute(
            "UPDATE stages SET status='failed', status_reason='tool_failed' "
            "WHERE stage_id='checkeuk'"
        )

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["evidence"]["checkeuk"]["status"] == "unavailable"
    assert bundle["evidence"]["checkeuk"]["stage_status"] == "failed"
    assert bundle["evidence"]["checkeuk"]["status_reason"] == (
        "registered_artifact_unavailable"
    )


def test_non_read_qc_statistics_do_not_supply_mapping_denominator(
    catalog: Path, tmp_path: Path
) -> None:
    with sqlite3.connect(catalog) as connection:
        connection.execute(
            "UPDATE read_stats SET stage_id='assemble', artifact_id='assembly-fasta'"
        )

    paths = export_visualization(catalog, tmp_path / "out")
    bundle = json.loads(paths[0].read_text(encoding="utf-8"))

    assert bundle["mapping"]["input_reads"] is None
    assert bundle["mapping"]["assembly_mapped_reads"] == 6
    assert bundle["mapping"]["status"] == "partial"
