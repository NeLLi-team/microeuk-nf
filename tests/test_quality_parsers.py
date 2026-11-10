"""Tests for native quality, taxonomy, and phenotype table parsing."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from protist_meta.quality_parsers import (
    GtdbtkCall,
    SymclatronEvidence,
    parse_quality,
    read_gtdbtk_calls,
)

CHECKEUK_FIELDS = (
    "genome",
    "status",
    "error_type",
    "error_message",
    "warning",
    "input_type",
    "completeness",
    "completeness_basis",
    "completeness_bom",
    "contamination",
    "contamination_type",
    "contamination_domains",
    "HGT_signal",
    "n_contigs_cont",
    "n_bp_cont",
    "cont_contig_contamination_pct",
    "contamination_status",
    "lineage",
    "completeness_note",
    "n_contigs",
    "assembly_size",
    "n50",
    "gc_pct",
    "clade_used",
    "pfam_status",
    "pfam_model",
    "pfam_union_families",
    "pfam_observed_families",
    "pfam_present_fraction",
    "pfam_completeness",
    "pfam_minus_bom_pp",
    "pfam_bom_agreement",
    "pfam_raw_domains",
    "pfam_retained_loci",
    "pfam_excluded_or_echo_domains",
    "pfam_translation_geometry",
    "pfam_warnings",
)
CHECKM1_FIELDS = (
    "Bin Id",
    "Marker lineage",
    "# genomes",
    "# markers",
    "# marker sets",
    "0",
    "1",
    "2",
    "3",
    "4",
    "5+",
    "Completeness",
    "Contamination",
    "Strain heterogeneity",
)
CHECKM2_FIELDS = (
    "Name",
    "Completeness",
    "Contamination",
    "Completeness_Model_Used",
    "Translation_Table_Used",
    "Coding_Density",
    "Contig_N50",
    "Average_Gene_Length",
    "Genome_Size",
    "GC_Content",
    "Total_Coding_Sequences",
    "Total_Contigs",
    "Max_Contig_Length",
    "Additional_Notes",
)
GVCLASS_FIELDS = (
    "query",
    "quality_models",
    "taxonomy_majority",
    "species_tree_nn_taxonomy",
    "taxonomy_confidence",
    "estimated_completeness",
    "completeness_query_support",
    "estimated_completeness_strategy",
    "completeness_model_group",
    "completeness_model_reliability",
    "estimated_contamination",
    "contamination_query_support",
    "estimated_contamination_strategy",
    "contamination_model_reliability",
    "contamination_type",
)
GTDBTK_FIELDS = (
    "user_genome",
    "classification",
    "closest_genome_reference",
    "closest_genome_reference_radius",
    "closest_genome_taxonomy",
    "closest_genome_ani",
    "closest_genome_af",
    "closest_placement_reference",
    "closest_placement_radius",
    "closest_placement_taxonomy",
    "closest_placement_ani",
    "closest_placement_af",
    "pplacer_taxonomy",
    "classification_method",
    "note",
    "other_related_references(genome_id,species_name,radius,ANI,AF)",
    "msa_percent",
    "translation_table",
    "red_value",
    "warnings",
)
SYMCLATRON_FIELDS = (
    "taxon_oid",
    "completeness_UNI56",
    "classification",
    "confidence",
    "passes_confidence_threshold",
    "classification_thresholded",
)
QUICKCLADE_FIELDS = (
    "#QueryName",
    "Q_GC",
    "Q_Bases",
    "Q_Contigs",
    "RefName",
    "R_TaxID",
    "R_GC",
    "R_Bases",
    "R_Contigs",
    "R_Level",
    "GCdif",
    "STRdif",
    "HHdif",
    "CAGAdif",
    "k3dif",
    "k4dif",
    "k5dif",
    "lineage",
    "ConfLevel",
    "Confidence",
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


def test_checkeuk_keeps_distinct_completeness_evidence(tmp_path: Path) -> None:
    """CheckEUK Pfam, BOM, inventory, and expected counts remain distinct."""
    summary = _write_tsv(
        tmp_path / "checkeuk_report.tsv",
        CHECKEUK_FIELDS,
        [
            {
                "genome": "/work/bin.01.fasta",
                "status": "ok",
                "completeness": "0",
                "completeness_basis": "pfam_union1085_ridge",
                "completeness_bom": "",
                "contamination": "0",
                "lineage": "sg__SAR;p__Ochrophyta",
                "completeness_note": "markers 0/74 (BOM reference)",
                "pfam_status": "scored",
                "pfam_union_families": "1085",
                "pfam_observed_families": "0",
                "pfam_completeness": "0",
            }
        ],
    )
    detail = _write_tsv(
        tmp_path / "checkeuk_report_detail.tsv",
        ("genome", "n_present_core", "n_present_expected", "n_expected"),
        [
            {
                "genome": "bin.01.fasta",
                "n_present_core": "0",
                "n_present_expected": "0",
                "n_expected": "74",
            }
        ],
    )

    parsed = parse_quality("CheckEUK", summary, detail)

    assert parsed.target_ids == ("bin.01",)
    assert parsed.qc == [
        {
            "record_id": "checkeuk:qc:bin.01",
            "target_kind": "bin",
            "target_id": "bin.01",
            "completeness_percent": 0.0,
            "completeness_basis": "pfam_union1085_ridge",
            "contamination_percent": 0.0,
            "contamination_basis": "CheckEUK report contamination",
            "checkeuk_pfam_ridge_completeness_percent": 0.0,
            "checkeuk_inventory74_present": 0,
            "checkeuk_expected_markers_present": 0,
            "checkeuk_expected_markers_total": 74,
            "checkeuk_pfam1085_distinct": 0,
        }
    ]
    assert [row["taxon_name"] for row in parsed.taxonomy] == ["SAR", "Ochrophyta"]
    assert "checkeuk_bom_completeness_percent" not in parsed.qc[0]
    assert all("unavailable without" not in gap for gap in parsed.evidence_gaps)


def test_checkeuk_filtered_target_remains_covered_without_qc(
    tmp_path: Path,
) -> None:
    """A native filtered row records why CheckEUK did not score the target."""
    summary = _write_tsv(
        tmp_path / "checkeuk_report.tsv",
        CHECKEUK_FIELDS,
        [
            {
                "genome": "small-bin.fa",
                "status": "filtered",
                "completeness": "",
                "completeness_basis": "",
                "completeness_bom": "",
                "contamination": "",
                "lineage": "",
                "completeness_note": "all 2 contigs below --min-contig 3000",
                "pfam_status": "not_run_filtered",
                "pfam_union_families": "",
                "pfam_observed_families": "",
                "pfam_completeness": "",
            }
        ],
    )
    detail = _write_tsv(
        tmp_path / "checkeuk_report_detail.tsv",
        ("genome", "n_present_core", "n_present_expected", "n_expected"),
        [
            {
                "genome": "small-bin.fa",
                "n_present_core": "",
                "n_present_expected": "",
                "n_expected": "",
            }
        ],
    )

    parsed = parse_quality("CheckEUK", summary, detail)

    assert parsed.target_ids == ("small-bin",)
    assert parsed.qc == []
    assert parsed.taxonomy == []
    assert any(
        "small-bin" in gap and "all 2 contigs below --min-contig 3000" in gap
        for gap in parsed.evidence_gaps
    )


@pytest.mark.parametrize("contamination", [0.0, 214.81, 300.0])
def test_checkm_native_reports_become_tool_specific_qc(
    tmp_path: Path, contamination: float
) -> None:
    """CheckM1 and CheckM2 keep their native model or lineage basis."""
    checkm1 = _write_tsv(
        tmp_path / "checkm1.tsv",
        CHECKM1_FIELDS,
        [
            {
                "Bin Id": "bin-a.fna",
                "Marker lineage": "Eukaryota",
                "Completeness": "0",
                "Contamination": str(contamination),
            }
        ],
    )
    checkm2 = _write_tsv(
        tmp_path / "quality_report.tsv",
        CHECKM2_FIELDS,
        [
            {
                "Name": "bin.2.fa",
                "Completeness": "91.25",
                "Contamination": "1.5",
                "Completeness_Model_Used": "Neural Network (Specific Model)",
            }
        ],
    )

    parsed1 = parse_quality("checkm1", checkm1)
    parsed2 = parse_quality("checkm2", checkm2)

    assert parsed1.target_ids == ("bin-a",)
    assert parsed1.qc[0]["completeness_percent"] == 0.0
    assert parsed1.qc[0]["contamination_percent"] == contamination
    assert parsed1.qc[0]["completeness_basis"] == "Eukaryota"
    assert parsed2.target_ids == ("bin.2",)
    assert parsed2.qc[0]["record_id"] == "checkm2:qc:bin.2"
    assert parsed2.qc[0]["completeness_basis"] == "Neural Network (Specific Model)"


def test_gvclass_retains_unknown_target_in_coverage(tmp_path: Path) -> None:
    """A GVClass row without estimates still contributes its exact target ID."""
    report = _write_tsv(
        tmp_path / "gvclass_summary.tsv",
        GVCLASS_FIELDS,
        [
            {
                "query": "/results/bin_3.fasta",
                "taxonomy_majority": "d_NCLDV;p_Nucleocytoviricota",
                "taxonomy_confidence": "low_support",
                "estimated_completeness": "nd",
                "completeness_query_support": "unassigned",
                "estimated_completeness_strategy": "gvog8 regression",
                "completeness_model_group": "NCLDV",
                "completeness_model_reliability": "unavailable",
                "estimated_contamination": "1.99",
                "contamination_query_support": "supported",
                "estimated_contamination_strategy": "marker duplication",
                "contamination_model_reliability": "unvalidated_legacy",
                "contamination_type": "none",
            },
            {
                "query": "bin_4.fa",
                "taxonomy_majority": "nd",
                "taxonomy_confidence": "nd",
                "estimated_completeness": "nd",
                "estimated_contamination": "nd",
            },
        ],
    )

    parsed = parse_quality("GVClass", report)

    assert parsed.target_ids == ("bin_3", "bin_4")
    assert len(parsed.qc) == 1
    assert parsed.qc[0]["contamination_percent"] == 1.99
    assert parsed.qc[0]["completeness_query_support"] == "unassigned"
    assert parsed.qc[0]["completeness_model_reliability"] == "unavailable"
    assert parsed.qc[0]["contamination_query_support"] == "supported"
    assert parsed.qc[0]["contamination_model_reliability"] == "unvalidated_legacy"
    assert [row["taxon_name"] for row in parsed.taxonomy] == [
        "NCLDV",
        "Nucleocytoviricota",
    ]
    assert {row["lineage_confidence"] for row in parsed.taxonomy} == {"low_support"}
    assert all(row["target_id"] != "bin_4" for row in parsed.qc)
    assert all(row["target_id"] != "bin_4" for row in parsed.taxonomy)
    assert any("retained native report" in gap for gap in parsed.evidence_gaps)


def test_gtdbtk_parses_each_rank_and_flags_native_evidence_gap(tmp_path: Path) -> None:
    """GTDB-Tk taxonomy survives related-reference fields above the CSV default."""
    report = _write_tsv(
        tmp_path / "gtdbtk.bac120.summary.tsv",
        GTDBTK_FIELDS,
        [
            {
                "user_genome": "bin_5",
                "classification": (
                    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria"
                ),
                "classification_method": "ANI/Placement",
                "closest_genome_ani": "98.2",
                "closest_genome_af": "0.91",
                "other_related_references(genome_id,species_name,radius,ANI,AF)": (
                    "GCF_000005845.2,s__Escherichia coli,95,98.2,0.91;" * 4096
                ),
                "red_value": "0.973",
            }
        ],
    )

    original_limit = csv.field_size_limit(131072)
    try:
        parsed = parse_quality("GTDB-Tk", report)
    finally:
        csv.field_size_limit(original_limit)

    assert parsed.target_ids == ("bin_5",)
    assert [row["rank"] for row in parsed.taxonomy] == [
        "domain",
        "phylum",
        "class",
    ]
    assert parsed.taxonomy[0]["record_id"] == "gtdbtk:taxonomy:bin_5:domain"
    assert any("ANI" in gap and "RED" in gap for gap in parsed.evidence_gaps)


BACTERIA = "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria"
APPLICABLE = {
    "eukaryotic_bins": frozenset(),
    "viral_bins": frozenset(),
    "gtdbtk": {"bin_6": GtdbtkCall(BACTERIA, "N/A")},
}


def _symclatron_report(tmp_path: Path, uni56: str = "80.5") -> Path:
    return _write_tsv(
        tmp_path / "symclatron_results.tsv",
        SYMCLATRON_FIELDS,
        [
            {
                "taxon_oid": "bin_6.faa",
                "completeness_UNI56": uni56,
                "classification": "Symbiont;Host-associated",
                "confidence": "0.9",
                "passes_confidence_threshold": "True",
                "classification_thresholded": "Symbiont;Host-associated",
            }
        ],
    )


def _symclatron_phenotype(
    tmp_path: Path,
    uni56: str = "80.5",
    **changes: object,
) -> dict[str, object]:
    evidence = SymclatronEvidence(**{**APPLICABLE, **changes})
    report = _symclatron_report(tmp_path, uni56)
    parsed = parse_quality("Symclatron", report, symclatron_evidence=evidence)
    assert len(parsed.phenotypes) == 1
    return parsed.phenotypes[0]


def test_symclatron_is_applicable_when_every_check_passes(tmp_path: Path) -> None:
    """A bacterial, non-viral, non-eukaryotic, complete bin is interpretable."""
    report = _symclatron_report(tmp_path)
    evidence = SymclatronEvidence(**APPLICABLE)

    parsed = parse_quality("Symclatron", report, symclatron_evidence=evidence)

    assert parsed.target_ids == ("bin_6",)
    assert parsed.qc[0]["completeness_percent"] == 80.5
    assert parsed.qc[0]["completeness_basis"] == "UNI56 marker inventory (56 markers)"
    assert parsed.phenotypes == [
        {
            "record_id": "symclatron:phenotype:bin_6",
            "bin_id": "bin_6",
            "raw_class": "Symbiont;Host-associated",
            "raw_score": 0.9,
            "threshold": 0.725,
            "passes_threshold": True,
            "thresholded_class": "Symbiont;Host-associated",
            "is_applicable": True,
            "applicability_reason": None,
        }
    ]
    assert any("GTDB-Tk warnings" in gap for gap in parsed.evidence_gaps)


@pytest.mark.parametrize(
    ("uni56", "changes", "reason"),
    [
        pytest.param(
            "80.5",
            {"eukaryotic_bins": frozenset({"bin_6"})},
            "CheckEUK assigns eukaryotic taxonomy",
            id="checkeuk-lineage",
        ),
        pytest.param(
            "80.5",
            {"viral_bins": frozenset({"bin_6"})},
            "GVClass majority is viral",
            id="gvclass-viral-majority",
        ),
        pytest.param(
            "80.5",
            {"gtdbtk": {"bin_6": GtdbtkCall(BACTERIA, "Domain Questionable")}},
            "GTDB-Tk flags the domain as questionable",
            id="gtdbtk-domain-questionable",
        ),
        pytest.param(
            "80.5",
            {"gtdbtk": {"bin_6": GtdbtkCall("Unclassified Bacteria", "N/A")}},
            "No GTDB-Tk bacterial or archaeal classification",
            id="gtdbtk-unclassified",
        ),
        pytest.param(
            "80.5",
            {"gtdbtk": {}},
            "No GTDB-Tk bacterial or archaeal classification",
            id="gtdbtk-missing",
        ),
        pytest.param(
            "49.99",
            {},
            "UNI56 completeness below 50 percent",
            id="uni56-below-floor",
        ),
    ],
)
def test_symclatron_reports_the_failed_check(
    tmp_path: Path,
    uni56: str,
    changes: dict[str, object],
    reason: str,
) -> None:
    """Each failed check makes the call inapplicable with its own reason."""
    phenotype = _symclatron_phenotype(tmp_path, uni56, **changes)

    assert phenotype["is_applicable"] is False
    assert phenotype["applicability_reason"] == reason
    assert phenotype["raw_class"] == "Symbiont;Host-associated"
    assert phenotype["raw_score"] == 0.9
    assert phenotype["passes_threshold"] is True


@pytest.mark.parametrize(
    "classification",
    ["d__Bacteria;p__Pseudomonadota", "d__Archaea;p__Asgardarchaeota"],
)
def test_symclatron_accepts_bacteria_archaea_and_the_uni56_floor(
    tmp_path: Path,
    classification: str,
) -> None:
    """Both prokaryotic domains pass and UNI56 equal to 50 is interpretable."""
    phenotype = _symclatron_phenotype(
        tmp_path,
        "50",
        gtdbtk={"bin_6": GtdbtkCall(classification, "N/A")},
    )

    assert phenotype["is_applicable"] is True
    assert phenotype["applicability_reason"] is None


def test_symclatron_reports_the_first_failed_check(tmp_path: Path) -> None:
    """Checks run in the study's order: CheckEUK, GVClass, domain, UNI56."""
    failing = {
        "eukaryotic_bins": frozenset({"bin_6"}),
        "viral_bins": frozenset({"bin_6"}),
        "gtdbtk": {"bin_6": GtdbtkCall("Unclassified Bacteria", "domain questionable")},
    }
    expected = [
        ("CheckEUK assigns eukaryotic taxonomy", {}),
        ("GVClass majority is viral", {"eukaryotic_bins": frozenset()}),
        (
            "GTDB-Tk flags the domain as questionable",
            {"eukaryotic_bins": frozenset(), "viral_bins": frozenset()},
        ),
    ]
    for reason, relaxed in expected:
        phenotype = _symclatron_phenotype(tmp_path, "10", **{**failing, **relaxed})
        assert phenotype["applicability_reason"] == reason

    unclassified = {
        **failing,
        "eukaryotic_bins": frozenset(),
        "viral_bins": frozenset(),
    }
    unclassified["gtdbtk"] = {"bin_6": GtdbtkCall("Unclassified Bacteria", "N/A")}
    phenotype = _symclatron_phenotype(tmp_path, "10", **unclassified)
    assert phenotype["applicability_reason"] == (
        "No GTDB-Tk bacterial or archaeal classification"
    )


def test_symclatron_requires_applicability_evidence(tmp_path: Path) -> None:
    """Symclatron cannot be parsed without the evidence that decides applicability."""
    report = _symclatron_report(tmp_path)

    with pytest.raises(ValueError, match="symclatron_evidence"):
        parse_quality("Symclatron", report)


def test_applicability_evidence_is_only_for_symclatron(tmp_path: Path) -> None:
    """Other tools reject Symclatron evidence rather than ignoring it."""
    report = _write_tsv(tmp_path / "checkm1.tsv", CHECKM1_FIELDS, [])

    with pytest.raises(ValueError, match="symclatron_evidence"):
        parse_quality(
            "CheckM1",
            report,
            symclatron_evidence=SymclatronEvidence(**APPLICABLE),
        )


def test_read_gtdbtk_calls_keeps_classification_and_warnings(tmp_path: Path) -> None:
    """GTDB-Tk calls come from every summary with their native warning text."""
    first = _write_tsv(
        tmp_path / "gtdbtk.bac120.summary.tsv",
        GTDBTK_FIELDS,
        [
            {
                "user_genome": "bin_1.fa",
                "classification": BACTERIA,
                "classification_method": "ani_screen",
                "warnings": "N/A",
            }
        ],
    )
    second = _write_tsv(
        tmp_path / "gtdbtk.ar53.summary.tsv",
        GTDBTK_FIELDS,
        [
            {
                "user_genome": "bin_2",
                "classification": "d__Archaea;p__Asgardarchaeota",
                "classification_method": "taxonomic classification defined by topology",
                "warnings": (
                    "Genome has more than 20% of markers with multiple hits; "
                    "Domain questionable"
                ),
            }
        ],
    )

    calls = read_gtdbtk_calls([first, second])

    assert calls == {
        "bin_1": GtdbtkCall(BACTERIA, "N/A"),
        "bin_2": GtdbtkCall(
            "d__Archaea;p__Asgardarchaeota",
            "Genome has more than 20% of markers with multiple hits; "
            "Domain questionable",
        ),
    }
    with pytest.raises(ValueError, match="duplicate GTDB-Tk target"):
        read_gtdbtk_calls([first, first])


def test_quickclade_uses_exact_basename_and_per_rank_scores(tmp_path: Path) -> None:
    """QuickClade IDs and confidence vectors map without substring matching."""
    report = _write_tsv(
        tmp_path / "quickclade.tsv",
        QUICKCLADE_FIELDS,
        [
            {
                "#QueryName": "/bins/bin_0.fa",
                "lineage": "d_NCLDV;k_Megaviricetes;p_Nucleocytoviricota",
                "ConfLevel": "kingdom:98.4",
                "Confidence": "d:100;k:98.4;p:75",
            },
            {
                "#QueryName": "bin_00.fa",
                "lineage": "d_NCLDV",
                "ConfLevel": "domain:99",
                "Confidence": "d:99",
            },
        ],
    )

    parsed = parse_quality("QuickClade", report)

    assert parsed.target_ids == ("bin_0", "bin_00")
    first = [row for row in parsed.taxonomy if row["target_id"] == "bin_0"]
    assert [row["score"] for row in first] == [100.0, 98.4, 75.0]
    assert first[1]["record_id"] == "quickclade:taxonomy:bin_0:kingdom"


def test_valid_header_only_report_is_empty(tmp_path: Path) -> None:
    """A valid native header may represent a successful zero-result stage."""
    report = _write_tsv(tmp_path / "quality_report.tsv", CHECKM2_FIELDS, [])

    parsed = parse_quality("checkm2", report)

    assert parsed.target_ids == ()
    assert parsed.qc == []
    assert parsed.taxonomy == []
    assert parsed.phenotypes == []


def test_missing_native_columns_are_rejected(tmp_path: Path) -> None:
    """A nonempty table missing its native contract fails at the boundary."""
    report = _write_tsv(
        tmp_path / "quality_report.tsv",
        ("Name", "Completeness"),
        [{"Name": "bin_1", "Completeness": "90"}],
    )

    with pytest.raises(ValueError, match="missing required columns"):
        parse_quality("checkm2", report)


def test_duplicate_normalized_target_ids_are_rejected(tmp_path: Path) -> None:
    """A filename and its base ID cannot create two records for one target."""
    report = _write_tsv(
        tmp_path / "quality_report.tsv",
        CHECKM2_FIELDS,
        [
            {
                "Name": "bin_1.fa",
                "Completeness": "90",
                "Contamination": "1",
                "Completeness_Model_Used": "model-a",
            },
            {
                "Name": "bin_1",
                "Completeness": "80",
                "Contamination": "2",
                "Completeness_Model_Used": "model-a",
            },
        ],
    )

    with pytest.raises(ValueError, match="duplicate CheckM2 target"):
        parse_quality("checkm2", report)


@pytest.mark.parametrize(
    ("tool", "fields", "identifier", "basis"),
    [
        ("checkm1", CHECKM1_FIELDS, "Bin Id", "Marker lineage"),
        ("checkm2", CHECKM2_FIELDS, "Name", "Completeness_Model_Used"),
    ],
)
@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("Completeness", "-0.1", "Completeness is below 0"),
        ("Completeness", "inf", "Completeness is not finite"),
        ("Completeness", "100.1", "Completeness exceeds 100"),
        ("Contamination", "-0.1", "Contamination is below 0"),
        ("Contamination", "inf", "Contamination is not finite"),
        ("Contamination", "nan", "Contamination is missing"),
    ],
)
def test_invalid_checkm_quality_value_is_rejected(
    tmp_path: Path,
    *,
    tool: str,
    fields: tuple[str, ...],
    identifier: str,
    basis: str,
    field: str,
    value: str,
    message: str,
) -> None:
    """CheckM requires finite nonnegative values and completeness at most 100."""
    report = _write_tsv(
        tmp_path / "quality_report.tsv",
        fields,
        [
            {
                identifier: "bin_1",
                "Completeness": "90",
                "Contamination": "1",
                basis: "model-a",
                field: value,
            }
        ],
    )

    with pytest.raises(ValueError, match=message):
        parse_quality(tool, report)


def test_checkeuk_detail_requires_exact_target_coverage(tmp_path: Path) -> None:
    """A supplied detail report must cover the same targets as the summary."""
    summary = _write_tsv(
        tmp_path / "checkeuk_report.tsv",
        CHECKEUK_FIELDS,
        [
            {
                "genome": "bin_1.fa",
                "status": "ok",
                "completeness": "80",
                "completeness_basis": "pfam_union1085_ridge",
                "contamination": "0",
                "lineage": "d__Eukaryota",
                "completeness_note": "markers 50/60",
                "pfam_union_families": "1085",
                "pfam_observed_families": "700",
                "pfam_completeness": "80",
            }
        ],
    )
    detail = _write_tsv(
        tmp_path / "checkeuk_report_detail.tsv",
        ("genome", "n_present_core", "n_present_expected", "n_expected"),
        [],
    )

    with pytest.raises(ValueError, match="summary and detail target sets differ"):
        parse_quality("checkeuk", summary, detail)


BUSCO_FIELDS = (
    "Input_file",
    "Dataset",
    "Complete",
    "Single",
    "Duplicated",
    "Fragmented",
    "Missing",
    "n_markers",
    "Scaffold N50",
)


def _busco_row(**changes: str) -> dict[str, str]:
    return {
        "Input_file": "bin_1.fna",
        "Dataset": "eukaryota_odb12",
        "Complete": "91.0",
        "Single": "89.0",
        "Duplicated": "2.0",
        "Fragmented": "3.0",
        "Missing": "6.0",
        "n_markers": "129",
        "Scaffold N50": "1000",
        **changes,
    }


def test_busco_maps_complete_and_duplicated_onto_the_qc_columns(
    tmp_path: Path,
) -> None:
    """Unused columns such as Single and Scaffold N50 are tolerated."""
    report = _write_tsv(tmp_path / "batch_summary.txt", BUSCO_FIELDS, [_busco_row()])

    parsed = parse_quality("busco", report)

    assert parsed.target_ids == ("bin_1",)
    assert parsed.qc == [
        {
            "record_id": "busco:qc:bin_1",
            "target_kind": "bin",
            "target_id": "bin_1",
            "completeness_percent": 91.0,
            "completeness_basis": "BUSCO eukaryota_odb12 complete (n=129)",
            "contamination_percent": 2.0,
            "contamination_basis": "BUSCO duplicated (not a contamination estimate)",
        }
    ]
    assert parsed.taxonomy == []
    assert "bin_1: BUSCO fragmented 3.0%, missing 6.0%" in parsed.evidence_gaps


@pytest.mark.parametrize("column", ["Complete", "Duplicated", "Fragmented", "Missing"])
def test_busco_rejects_a_percentage_above_100(tmp_path: Path, column: str) -> None:
    report = _write_tsv(
        tmp_path / "batch_summary.txt", BUSCO_FIELDS, [_busco_row(**{column: "100.5"})]
    )

    with pytest.raises(ValueError, match=f"{column} exceeds 100"):
        parse_quality("busco", report)


BUSCO_EMPTY_METRICS = {
    "Complete": "",
    "Single": "",
    "Duplicated": "",
    "Fragmented": "",
    "Missing": "",
    "n_markers": "",
    "Scaffold N50": "",
}


def _busco_failure_row(input_file: str, dataset: str) -> dict[str, str]:
    return _busco_row(Input_file=input_file, Dataset=dataset, **BUSCO_EMPTY_METRICS)


@pytest.mark.parametrize("dataset", ["No genes found", "Run failed; check logs"])
def test_busco_failure_row_is_observed_without_a_qc_row(
    tmp_path: Path, dataset: str
) -> None:
    report = _write_tsv(
        tmp_path / "batch_summary.txt",
        BUSCO_FIELDS,
        [_busco_row(), _busco_failure_row("bin_2.fna", dataset)],
    )

    parsed = parse_quality("busco", report)

    assert parsed.target_ids == ("bin_1", "bin_2")
    assert [row["target_id"] for row in parsed.qc] == ["bin_1"]
    assert f"bin_2: BUSCO {dataset}" in parsed.evidence_gaps
    assert "bin_1: BUSCO fragmented 3.0%, missing 6.0%" in parsed.evidence_gaps


def test_busco_table_of_only_failure_rows_has_no_qc_rows(tmp_path: Path) -> None:
    report = _write_tsv(
        tmp_path / "batch_summary.txt",
        BUSCO_FIELDS,
        [
            _busco_failure_row("bin_1.fna", "No genes found"),
            _busco_failure_row("bin_2.fna", "Run failed; check logs"),
        ],
    )

    parsed = parse_quality("busco", report)

    assert parsed.qc == []
    assert parsed.target_ids == ("bin_1", "bin_2")


def test_busco_rejects_a_normal_row_with_an_empty_metric(tmp_path: Path) -> None:
    report = _write_tsv(
        tmp_path / "batch_summary.txt", BUSCO_FIELDS, [_busco_row(Complete="")]
    )

    with pytest.raises(ValueError, match="Complete is missing"):
        parse_quality("busco", report)


def test_busco_rejects_a_failure_row_with_a_metric_value(tmp_path: Path) -> None:
    report = _write_tsv(
        tmp_path / "batch_summary.txt",
        BUSCO_FIELDS,
        [_busco_failure_row("bin_1.fna", "No genes found") | {"Complete": "5.0"}],
    )

    with pytest.raises(ValueError, match="BUSCO failure row for bin_1"):
        parse_quality("busco", report)


@pytest.mark.parametrize("prefix", ["#", "# "])
def test_busco_accepts_a_hash_on_the_header_line(tmp_path: Path, prefix: str) -> None:
    report = tmp_path / "batch_summary.txt"
    report.write_text(
        prefix
        + "\t".join(BUSCO_FIELDS)
        + "\n"
        + "\t".join(_busco_row().values())
        + "\n",
        encoding="utf-8",
    )

    assert parse_quality("busco", report).target_ids == ("bin_1",)


def test_busco_names_the_missing_column(tmp_path: Path) -> None:
    report = _write_tsv(
        tmp_path / "batch_summary.txt",
        tuple(field for field in BUSCO_FIELDS if field != "Missing"),
        [{key: value for key, value in _busco_row().items() if key != "Missing"}],
    )

    with pytest.raises(ValueError, match=r"missing required columns \['Missing'\]"):
        parse_quality("busco", report)
