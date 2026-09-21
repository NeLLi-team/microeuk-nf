"""Tests for exclusive gene-calling input preparation and coordinate lift-over."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from protist_meta.gene_inputs import (
    lift_prodigal_gene,
    load_prodigal_origins,
    prepare_gene_inputs,
)
from protist_meta.sequences import read_fasta


def _routing(root: Path, bins: dict[str, dict[str, str]]) -> Path:
    for route in ("prokaryotic", "viral", "eukaryotic", "unresolved"):
        (root / route).mkdir(parents=True)
    for route, route_bins in bins.items():
        for bin_id, fasta in route_bins.items():
            (root / route / f"{bin_id}.fna").write_text(fasta, encoding="utf-8")
    return root


def _viral(
    root: Path,
    rows: list[dict[str, str]],
    fasta: str,
) -> Path:
    summary_dir = root / "genomad/input_summary"
    summary_dir.mkdir(parents=True)
    with (summary_dir / "input_virus_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("seq_name", "length", "topology", "coordinates"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    (summary_dir / "input_virus.fna").write_text(fasta, encoding="utf-8")
    return root


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_prepare_excludes_called_parents_and_lifts_unresolved_provirus(
    tmp_path: Path,
) -> None:
    routing = _routing(
        tmp_path / "routing",
        {
            "prokaryotic": {"prok1": ">called_prok\n" + "A" * 300 + "\n"},
            "viral": {"viral1": ">called_viral\n" + "C" * 200 + "\n"},
            "eukaryotic": {"euk1": ">called_euk\n" + "G" * 250 + "\n"},
            "unresolved": {"uncertain1": ">uncertain\n" + "T" * 300 + "\n"},
        },
    )
    viral = _viral(
        tmp_path / "viral",
        [
            {
                "seq_name": "called_prok|provirus_101_200",
                "length": "100",
                "topology": "Provirus",
                "coordinates": "101-200",
            },
            {
                "seq_name": "called_euk|provirus_51_150",
                "length": "100",
                "topology": "Provirus",
                "coordinates": "51-150",
            },
            {
                "seq_name": "uncertain|provirus_51_150",
                "length": "100",
                "topology": "Provirus",
                "coordinates": "51-150",
            },
            {
                "seq_name": "unbinned_full",
                "length": "120",
                "topology": "No terminal repeats",
                "coordinates": "NA",
            },
        ],
        (
            ">called_prok|provirus_101_200\n"
            + "A" * 100
            + "\n>called_euk|provirus_51_150\n"
            + "C" * 100
            + "\n>uncertain|provirus_51_150\n"
            + "G" * 100
            + "\n>unbinned_full\n"
            + "T" * 120
            + "\n"
        ),
    )

    sources_path = prepare_gene_inputs(routing, viral, tmp_path / "prepared")
    sources = _rows(sources_path)
    origins = _rows(tmp_path / "prepared/origins.tsv")

    assert [row["category"] for row in sources] == [
        "prokaryotic_bin",
        "viral_bin",
        "genomad_unbinned",
    ]
    assert len([row for row in sources if row["bin_id"] == "prok1"]) == 1
    exclusions = {
        row["parent_contig"]: row for row in origins if row["status"] == "excluded"
    }
    assert exclusions["called_prok"]["bin_id"] == "prok1"
    assert exclusions["called_euk"]["bin_id"] == "euk1"
    assert {row["reason"] for row in exclusions.values()} == {
        "parent_assigned_to_called_bin"
    }

    genomad_source = sources[-1]
    staged = list(read_fasta(tmp_path / "prepared" / genomad_source["fasta"]))
    assert len(staged) == 2
    assert all("|" not in record.identifier for record in staged)
    included = {
        row["parent_contig"]: row
        for row in origins
        if row["source_id"] == "genomad-unbinned" and row["status"] == "included"
    }
    assert included["uncertain"]["bin_id"] == "uncertain1"
    assert included["unbinned_full"]["start"] == "1"
    assert included["unbinned_full"]["end"] == "120"

    lifted = lift_prodigal_gene(
        {
            "gene_id": "prodigal-gv:genomad-unbinned:native_1",
            "contig_id": included["uncertain"]["contig_id"],
            "bin_id": None,
            "start": 4,
            "end": 90,
            "strand": "plus",
            "phase": 0,
            "genetic_code": 11,
        },
        load_prodigal_origins(tmp_path / "prepared" / genomad_source["origin_map"]),
    )
    assert lifted["contig_id"] == "uncertain"
    assert lifted["bin_id"] == "uncertain1"
    assert lifted["start"] == 54
    assert lifted["end"] == 140


def test_lift_prodigal_genes_reuses_parsed_origins(tmp_path: Path) -> None:
    origin_path = tmp_path / "origins.tsv"
    origin_path.write_text(
        "contig_id\tparent_contig\tstart\tend\tsource_id\tbin_id\tcategory\t"
        "status\treason\n"
        "fragment\tparent\t101\t400\tsource\tbin1\tgenomad_provirus\tincluded\t\n"
        "whole\twhole\t1\t600\tsource\t\tgenomad_full\tincluded\t\n",
        encoding="utf-8",
    )
    origins = load_prodigal_origins(origin_path)
    origin_path.unlink()

    first = lift_prodigal_gene(
        {"contig_id": "fragment", "start": 4, "end": 90}, origins
    )
    second = lift_prodigal_gene(
        {"contig_id": "fragment", "start": 100, "end": 300}, origins
    )
    whole = lift_prodigal_gene({"contig_id": "whole", "start": 1, "end": 600}, origins)

    assert first == {"contig_id": "parent", "start": 104, "end": 190, "bin_id": "bin1"}
    assert second == {
        "contig_id": "parent",
        "start": 200,
        "end": 400,
        "bin_id": "bin1",
    }
    assert whole == {"contig_id": "whole", "start": 1, "end": 600, "bin_id": None}


def test_load_prodigal_origins_rejects_duplicate_contigs(tmp_path: Path) -> None:
    origin_path = tmp_path / "origins.tsv"
    origin_path.write_text(
        "contig_id\tparent_contig\tstart\tend\tsource_id\tbin_id\tcategory\t"
        "status\treason\n"
        "fragment\tparent\t101\t400\tsource\tbin1\tgenomad_provirus\tincluded\t\n"
        "fragment\tother\t1\t300\tsource\tbin2\tgenomad_provirus\tincluded\t\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="repeats identifier: fragment"):
        load_prodigal_origins(origin_path)


@pytest.mark.parametrize("header", [">", "> \t"])
def test_read_fasta_rejects_empty_identifier(tmp_path: Path, header: str) -> None:
    fasta = tmp_path / "invalid.fna"
    fasta.write_text(f"{header}\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="FASTA header lacks an identifier"):
        list(read_fasta(fasta))


def test_prepare_keeps_nonoverlapping_proviruses_with_unique_safe_ids(
    tmp_path: Path,
) -> None:
    routing = _routing(tmp_path / "routing", {})
    viral = _viral(
        tmp_path / "viral",
        [
            {
                "seq_name": "parent|provirus_1_50",
                "length": "50",
                "topology": "Provirus",
                "coordinates": "1-50",
            },
            {
                "seq_name": "parent|provirus_101_150",
                "length": "50",
                "topology": "Provirus",
                "coordinates": "101-150",
            },
        ],
        ">parent|provirus_1_50\n"
        + "A" * 50
        + "\n>parent|provirus_101_150\n"
        + "C" * 50
        + "\n",
    )

    sources = _rows(prepare_gene_inputs(routing, viral, tmp_path / "prepared"))
    records = list(read_fasta(tmp_path / "prepared" / sources[0]["fasta"]))

    assert len(records) == 2
    assert len({record.identifier for record in records}) == 2
    assert all(_identifier_safe(record.identifier) for record in records)


def test_prepare_rejects_duplicate_genomad_identifiers(tmp_path: Path) -> None:
    routing = _routing(tmp_path / "routing", {})
    viral = _viral(
        tmp_path / "viral",
        [
            {
                "seq_name": "duplicate",
                "length": "4",
                "topology": "No terminal repeats",
                "coordinates": "NA",
            }
        ],
        ">duplicate\nAAAA\n>duplicate\nCCCC\n",
    )

    with pytest.raises(ValueError, match="repeats a sequence identifier"):
        prepare_gene_inputs(routing, viral, tmp_path / "prepared")


def _identifier_safe(value: str) -> bool:
    return value.replace("_", "").replace("-", "").isalnum()
