"""End-to-end regression test for core result collection."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sqlite3
from pathlib import Path

from protist_meta.catalog import build_catalog, validate_bundle
from protist_meta.collect import collect_records


def _write_tsv(
    path: Path,
    fields: tuple[str, ...],
    rows: list[dict[str, str]],
) -> None:
    """Write one small external-tool table."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_stage(path: Path, stage: str, outputs: dict[str, str]) -> None:
    """Write an explicit completed-stage replay manifest."""
    manifest = {
        "stage": stage,
        "status": "completed",
        "tool": "fixture",
        "tool_version": "test",
        "command": "portable regression fixture",
        "outputs": outputs,
    }
    (path / "stage.json").write_text(json.dumps(manifest) + "\n", encoding="utf-8")


def test_collects_preflight_receipts_without_duplicate_manifest(tmp_path: Path) -> None:
    """Preflight reference receipts are hashed once alongside its manifest."""
    reads = tmp_path / "reads.fastq"
    reads.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    preflight = tmp_path / "preflight"
    preflight.mkdir()
    info = preflight / "famdb-info.txt"
    exports = preflight / "famdb-exports.tsv"
    info.write_text("Database : Dfam\nVersion : 4.0\n", encoding="utf-8")
    exports.write_text("export\trecords\nfasta_all\t1\n", encoding="utf-8")
    _write_stage(
        preflight,
        "preflight",
        {
            "report": "stage.json",
            "famdb_info": info.name,
            "famdb_exports": exports.name,
        },
    )
    sample = tmp_path / "sample.json"
    sample.write_text(
        json.dumps(
            {
                "sample_id": "sample1",
                "run_id": "run1",
                "platform": "ONT",
                "reads": str(reads),
                "workflow_version": "0.1.0",
                "source_revision": "test-source",
            }
        ),
        encoding="utf-8",
    )

    bundle_path = collect_records(sample, [preflight], tmp_path / "bundle.json")
    bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
    artifacts = [
        artifact
        for artifact in bundle["artifacts"]
        if artifact["producer_stage_id"] == "preflight"
    ]

    assert len(artifacts) == 3
    assert [Path(artifact["path"]).name for artifact in artifacts].count(
        "stage.json"
    ) == 1
    receipts = {
        Path(artifact["path"]).name: artifact
        for artifact in artifacts
        if Path(artifact["path"]).name != "stage.json"
    }
    assert set(receipts) == {"famdb-info.txt", "famdb-exports.tsv"}
    for path in (info, exports):
        assert (
            receipts[path.name]["sha256"]
            == hashlib.sha256(path.read_bytes()).hexdigest()
        )


def test_collect_core_outputs_build_valid_catalog(tmp_path: Path) -> None:
    """Core native files retain counts, links, and exact expected keys."""
    reads = tmp_path / "reads.fastq"
    reads.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    stages = tmp_path / "stages"
    read_qc = stages / "01_read_qc"
    assembly = stages / "02_assembly"
    binning = stages / "03_binning"
    bins = binning / "bins"
    for directory in (read_qc, assembly, bins):
        directory.mkdir(parents=True)

    stats_fields = ("file", "format", "type", "num_seqs", "sum_len")
    stats = [
        {
            "file": "reads.fastq",
            "format": "FASTQ",
            "type": "DNA",
            "num_seqs": "1",
            "sum_len": "4",
        }
    ]
    _write_tsv(read_qc / "raw_stats.tsv", stats_fields, stats)
    _write_tsv(read_qc / "filtered_stats.tsv", stats_fields, stats)
    with gzip.open(read_qc / "filtered.fastq.gz", "wt", encoding="utf-8") as handle:
        handle.write(reads.read_text(encoding="utf-8"))
    _write_stage(
        read_qc,
        "read_qc",
        {
            "reads": "filtered.fastq.gz",
            "raw_stats": "raw_stats.tsv",
            "filtered_stats": "filtered_stats.tsv",
        },
    )

    fasta = ">ctg1\nACGT\n"
    (assembly / "assembly_primary.fa").write_text(fasta, encoding="utf-8")
    _write_stage(assembly, "assembly", {"assembly": "assembly_primary.fa"})

    (bins / "bin_0.fa").write_text(fasta, encoding="utf-8")
    _write_tsv(
        binning / "coverage.tsv",
        (
            "#rname",
            "startpos",
            "endpos",
            "numreads",
            "covbases",
            "coverage",
            "meandepth",
        ),
        [
            {
                "#rname": "ctg1",
                "startpos": "1",
                "endpos": "4",
                "numreads": "1",
                "covbases": "4",
                "coverage": "100",
                "meandepth": "2",
            }
        ],
    )
    _write_tsv(
        binning / "quickclade.tsv",
        ("#QueryName", "R_TaxID", "lineage"),
        [
            {
                "#QueryName": "bin_0.fa",
                "R_TaxID": "123",
                "lineage": "d__Eukaryota;k__Fungi",
            }
        ],
    )
    _write_stage(
        binning,
        "binning",
        {
            "coverage": "coverage.tsv",
            "bins": "bins",
            "quickclade": "quickclade.tsv",
        },
    )

    sample = tmp_path / "sample.json"
    sample.write_text(
        json.dumps(
            {
                "sample_id": "sample1",
                "run_id": "run1",
                "platform": "ONT",
                "reads": str(reads),
                "genetic_code": "auto",
                "workflow_version": "0.1.0",
                "source_revision": "test-source",
                "run_mode": "core",
            }
        ),
        encoding="utf-8",
    )
    bundle_path = tmp_path / "bundle.json"
    database_path = tmp_path / "catalog.sqlite"

    collect_records(sample, sorted(stages.iterdir()), bundle_path)
    bundle = validate_bundle(bundle_path)
    build_catalog(bundle_path, database_path)

    assert len(bundle.contigs) == 1
    assert bundle.contigs[0].mean_depth_x == 2.0
    assert len(bundle.memberships) == 1
    assert bundle.memberships[0].contig_id == "ctg1"
    completed = {
        stage.stage_id: stage for stage in bundle.stages if stage.status == "completed"
    }
    assert set(completed) == {"read_qc", "assembly", "binning"}
    assert completed["read_qc"].command == "portable regression fixture"
    assert set(completed["binning"].expected_result_keys) == {
        "bins:bin_0",
        "memberships:bin_0:ctg1",
        "taxonomy:quickclade:bin_0",
    }
    with sqlite3.connect(database_path) as connection:
        assert connection.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute("SELECT COUNT(*) FROM contigs").fetchone() == (1,)
