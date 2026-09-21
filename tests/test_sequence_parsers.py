"""Tests for native SSUextract, geNomad, and CheckV table parsing."""

from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path

import pytest

from protist_meta.sequence_parsers import (
    ParsedSequenceRows,
    parse_checkv,
    parse_genomad,
    parse_ssu,
)

SSU_FIELDS = (
    "name",
    "sample",
    "model",
    "length",
    "coordinates",
    "strand",
    "sequence_type",
    "contig_name",
    "blast_sseqid",
    "blast_pident",
    "blast_length",
    "blast_bitscore",
    "is_assembled",
    "reference_source",
    "taxonomy",
    "taxonomy_source",
    "taxonomy_assignment_method",
    "reference_identifiers",
    "reference_versions",
    "reference_taxonomy",
    "reference_taxonomy_source",
    "blast_query_coverage",
)
GENOMAD_FIELDS = (
    "seq_name",
    "length",
    "topology",
    "coordinates",
    "n_genes",
    "genetic_code",
    "virus_score",
    "fdr",
    "n_hallmarks",
    "marker_enrichment",
    "taxonomy",
)
CHECKV_FIELDS = (
    "contig_id",
    "contig_length",
    "provirus",
    "proviral_length",
    "gene_count",
    "viral_genes",
    "host_genes",
    "checkv_quality",
    "miuvig_quality",
    "completeness",
    "completeness_method",
    "contamination",
    "kmer_freq",
    "warnings",
)


def _write_tsv(
    path: Path,
    fields: tuple[str, ...],
    rows: list[dict[str, str]],
) -> Path:
    """Write a small native-shaped TSV fixture."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_parse_ssu_normalizes_coordinates_and_blast_units(tmp_path: Path) -> None:
    """SSUextract coordinates stay inclusive and query coverage becomes percent."""
    report = _write_tsv(
        tmp_path / "cmsearch_summary.tsv",
        SSU_FIELDS,
        [
            {
                "name": "sample|RF01960|ctg_1|101-200|strand_-|simple",
                "sample": "sample",
                "model": "RF01960",
                "length": "100",
                "coordinates": "101-200",
                "strand": "-",
                "sequence_type": "simple",
                "contig_name": "ctg_1",
                "blast_sseqid": "REF_PR2_AB123.1",
                "blast_pident": "99.5",
                "blast_length": "75",
                "blast_bitscore": "140",
                "is_assembled": "False",
                "reference_source": "PR2",
                "taxonomy": "Unclassified",
                "taxonomy_source": "FINAL_SOURCE_SENTINEL",
                "taxonomy_assignment_method": "no_runtime_calibration",
                "reference_identifiers": "PR2:AB123.1",
                "reference_versions": "PR2:5.0",
                "reference_taxonomy": (
                    '[{"assignment_method":"native","taxonomy":"Eukaryota;TSAR;'
                    'Stramenopiles","taxonomy_source":"PR2"}]'
                ),
                "reference_taxonomy_source": "SIBLING_SOURCE_SENTINEL",
                "blast_query_coverage": "0.75",
            },
            {
                "name": "sample|RF00177|ctg_2|1-80|strand_+|simple",
                "sample": "sample",
                "model": "RF00177",
                "length": "80",
                "coordinates": "1-80",
                "strand": "+",
                "sequence_type": "simple",
                "contig_name": "ctg_2",
                "taxonomy": "Unclassified",
            },
        ],
    )

    parsed = parse_ssu(report)

    assert parsed.target_ids == (
        "sample|RF01960|ctg_1|101-200|strand_-|simple",
        "sample|RF00177|ctg_2|1-80|strand_+|simple",
    )
    assert parsed.ssu == [
        {
            "record_id": "ssuextract:ssu:ctg_1:101-200:RF01960",
            "contig_id": "ctg_1",
            "ssu_type": "rrna_18s",
            "start": 101,
            "end": 200,
            "strand": "minus",
            "hit_accession": "REF_PR2_AB123.1",
            "reference_identifiers": "PR2:AB123.1",
            "reference_source": "PR2",
            "reference_versions": "PR2:5.0",
            "reference_taxonomy": (
                '[{"assignment_method":"native","taxonomy":"Eukaryota;TSAR;'
                'Stramenopiles","taxonomy_source":"PR2"}]'
            ),
            "taxon_name": "Unclassified",
            "taxonomy_assignment_method": "no_runtime_calibration",
            "identity_percent": 99.5,
            "query_coverage_percent": 75.0,
        },
        {
            "record_id": "ssuextract:ssu:ctg_2:1-80:RF00177",
            "contig_id": "ctg_2",
            "ssu_type": "rrna_16s",
            "start": 1,
            "end": 80,
            "strand": "plus",
            "taxon_name": "Unclassified",
        },
    ]


def test_parse_genomad_maps_whole_virus_and_exact_provirus_parent(
    tmp_path: Path,
) -> None:
    """A provirus on ctg_10 cannot be assigned to the ctg_1 substring."""
    report = _write_tsv(
        tmp_path / "input_virus_summary.tsv",
        GENOMAD_FIELDS,
        [
            {
                "seq_name": "ctg_1",
                "length": "1000",
                "topology": "No terminal repeats",
                "coordinates": "NA",
                "n_genes": "0",
                "genetic_code": "11",
                "virus_score": "0",
                "fdr": "NA",
                "n_hallmarks": "0",
                "marker_enrichment": "0",
                "taxonomy": "Unclassified",
            },
            {
                "seq_name": "ctg_10|provirus_101_200",
                "length": "100",
                "topology": "Provirus",
                "coordinates": "101-200",
                "n_genes": "4",
                "genetic_code": "11",
                "virus_score": "0.9",
                "fdr": "0.01",
                "n_hallmarks": "2",
                "marker_enrichment": "3.5",
                "taxonomy": (
                    "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;"
                    "Caudoviricetes;;Demerecviridae"
                ),
            },
        ],
    )

    parsed = parse_genomad(report, {"ctg_1": 1000, "ctg_10": 2000})

    assert parsed.target_ids == ("ctg_1", "ctg_10|provirus_101_200")
    assert parsed.viruses == [
        {
            "record_id": "genomad:virus:ctg_1",
            "contig_id": "ctg_1",
            "start": 1,
            "end": 1000,
            "strand": "unknown",
            "score": 0.0,
            "viral_taxon": "Unclassified",
        },
        {
            "record_id": "genomad:virus:ctg_10:101-200",
            "contig_id": "ctg_10",
            "start": 101,
            "end": 200,
            "strand": "unknown",
            "score": 0.9,
            "viral_taxon": (
                "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;"
                "Caudoviricetes;;Demerecviridae"
            ),
        },
    ]
    family = [row for row in parsed.taxonomy if row["rank"] == "family"]
    assert family == [
        {
            "record_id": ("genomad:taxonomy:genomad:virus:ctg_10:101-200:family"),
            "target_kind": "virus",
            "target_id": "genomad:virus:ctg_10:101-200",
            "rank": "family",
            "taxon_name": "Demerecviridae",
        }
    ]
    assert any("genetic code" in gap for gap in parsed.evidence_gaps)


def test_parse_genomad_accepts_finite_negative_marker_enrichment(
    tmp_path: Path,
) -> None:
    """Signed native marker enrichment remains valid for a called virus."""
    report = _write_tsv(
        tmp_path / "input_virus_summary.tsv",
        GENOMAD_FIELDS,
        [
            {
                "seq_name": "ctg_1",
                "length": "1000",
                "topology": "No terminal repeats",
                "coordinates": "NA",
                "n_genes": "4",
                "genetic_code": "11",
                "virus_score": "0.91",
                "fdr": "NA",
                "n_hallmarks": "0",
                "marker_enrichment": "-3.7554",
                "taxonomy": "Unclassified",
            }
        ],
    )

    parsed = parse_genomad(report, {"ctg_1": 1000})

    assert parsed.target_ids == ("ctg_1",)
    assert parsed.viruses[0]["record_id"] == "genomad:virus:ctg_1"


def test_parse_genomad_rejects_nonfinite_marker_enrichment(tmp_path: Path) -> None:
    """Relaxing the sign constraint does not admit nonfinite evidence."""
    report = _write_tsv(
        tmp_path / "input_virus_summary.tsv",
        GENOMAD_FIELDS,
        [
            {
                "seq_name": "ctg_1",
                "length": "1000",
                "topology": "No terminal repeats",
                "coordinates": "NA",
                "n_genes": "4",
                "genetic_code": "11",
                "virus_score": "0.91",
                "fdr": "NA",
                "n_hallmarks": "0",
                "marker_enrichment": "-inf",
                "taxonomy": "Unclassified",
            }
        ],
    )

    with pytest.raises(ValueError, match="marker_enrichment is not finite"):
        parse_genomad(report, {"ctg_1": 1000})


def test_parse_checkv_preserves_zero_and_null_predictions(tmp_path: Path) -> None:
    """CheckV zero estimates remain numbers while NA completeness stays absent."""
    report = _write_tsv(
        tmp_path / "quality_summary.tsv",
        CHECKV_FIELDS,
        [
            {
                "contig_id": "ctg_1",
                "contig_length": "1000",
                "provirus": "No",
                "proviral_length": "NA",
                "gene_count": "3",
                "viral_genes": "1",
                "host_genes": "0",
                "checkv_quality": "Low-quality",
                "miuvig_quality": "Genome-fragment",
                "completeness": "0",
                "completeness_method": "HMM-based (lower-bound)",
                "contamination": "0",
                "kmer_freq": "1",
                "warnings": "",
            },
            {
                "contig_id": "ctg_10|provirus_101_200",
                "contig_length": "100",
                "provirus": "Yes",
                "proviral_length": "95",
                "gene_count": "2",
                "viral_genes": "1",
                "host_genes": "1",
                "checkv_quality": "Not-determined",
                "miuvig_quality": "Genome-fragment",
                "completeness": "NA",
                "completeness_method": "NA",
                "contamination": "5.7",
                "kmer_freq": "1",
                "warnings": "",
            },
        ],
    )
    virus_ids = {
        "ctg_1": "genomad:virus:ctg_1",
        "ctg_10|provirus_101_200": "genomad:virus:ctg_10:101-200",
    }

    parsed = parse_checkv(report, virus_ids)

    assert parsed.target_ids == ("ctg_1", "ctg_10|provirus_101_200")
    assert parsed.qc == [
        {
            "record_id": "checkv:qc:genomad:virus:ctg_1",
            "target_kind": "virus",
            "target_id": "genomad:virus:ctg_1",
            "completeness_percent": 0.0,
            "completeness_basis": "HMM-based (lower-bound)",
            "contamination_percent": 0.0,
            "contamination_basis": "CheckV predicted host-region fraction",
        },
        {
            "record_id": "checkv:qc:genomad:virus:ctg_10:101-200",
            "target_kind": "virus",
            "target_id": "genomad:virus:ctg_10:101-200",
            "contamination_percent": 5.7,
            "contamination_basis": "CheckV predicted host-region fraction",
        },
    ]


@pytest.mark.parametrize(
    ("parser", "fields", "arguments"),
    [
        (parse_ssu, SSU_FIELDS, ()),
        (parse_genomad, GENOMAD_FIELDS, ({},)),
        (parse_checkv, CHECKV_FIELDS, ({},)),
    ],
)
def test_native_header_only_reports_are_valid(
    tmp_path: Path,
    parser: Callable[..., ParsedSequenceRows],
    fields: tuple[str, ...],
    arguments: tuple[object, ...],
) -> None:
    """Each tool can complete with a native header and no result rows."""
    report = _write_tsv(tmp_path / "empty.tsv", fields, [])

    parsed = parser(report, *arguments)

    assert parsed.target_ids == ()
    assert parsed.ssu == []
    assert parsed.viruses == []
    assert parsed.qc == []


def test_genomad_rejects_coordinate_suffix_disagreement(tmp_path: Path) -> None:
    """A provirus coordinate column cannot override its exact native ID."""
    report = _write_tsv(
        tmp_path / "input_virus_summary.tsv",
        GENOMAD_FIELDS,
        [
            {
                "seq_name": "ctg_1|provirus_101_200",
                "length": "100",
                "topology": "Provirus",
                "coordinates": "101-201",
                "n_genes": "1",
                "genetic_code": "11",
                "virus_score": "0.9",
                "fdr": "NA",
                "n_hallmarks": "1",
                "marker_enrichment": "2",
                "taxonomy": "Unclassified",
            }
        ],
    )

    with pytest.raises(ValueError, match="provirus coordinates disagree"):
        parse_genomad(report, {"ctg_1": 1000})


def test_genomad_rejects_duplicate_native_ids(tmp_path: Path) -> None:
    """Two geNomad rows cannot create the same viral record."""
    row = {
        "seq_name": "ctg_1",
        "length": "1000",
        "topology": "DTR",
        "coordinates": "NA",
        "n_genes": "2",
        "genetic_code": "11",
        "virus_score": "0.9",
        "fdr": "NA",
        "n_hallmarks": "1",
        "marker_enrichment": "2",
        "taxonomy": "Unclassified",
    }
    report = _write_tsv(
        tmp_path / "input_virus_summary.tsv",
        GENOMAD_FIELDS,
        [row, row],
    )

    with pytest.raises(ValueError, match="duplicate geNomad target"):
        parse_genomad(report, {"ctg_1": 1000})


def test_checkv_requires_exact_genomad_target_coverage(tmp_path: Path) -> None:
    """A CheckV report cannot omit a geNomad virus supplied in its mapping."""
    report = _write_tsv(tmp_path / "quality_summary.tsv", CHECKV_FIELDS, [])

    with pytest.raises(ValueError, match="target sets differ"):
        parse_checkv(report, {"ctg_1": "genomad:virus:ctg_1"})


def test_missing_native_columns_are_rejected(tmp_path: Path) -> None:
    """A nonempty table with a partial native header fails at the boundary."""
    report = _write_tsv(
        tmp_path / "quality_summary.tsv",
        ("contig_id", "completeness"),
        [{"contig_id": "ctg_1", "completeness": "90"}],
    )

    with pytest.raises(ValueError, match="missing required columns"):
        parse_checkv(report, {"ctg_1": "genomad:virus:ctg_1"})
