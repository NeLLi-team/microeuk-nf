import csv
import json
from pathlib import Path

import pytest

from protist_meta.cli import main
from protist_meta.collect import StageManifest
from protist_meta.routing import route_bins

FALLBACK_REASON = (
    "GVClass supports a eukaryotic route; CheckEUK has too few markers to name "
    "a lineage"
)


def _write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _evidence_dirs(root: Path) -> dict[str, Path]:
    return {
        name: root / name
        for name in (
            "quickclade",
            "checkeuk",
            "gvclass",
            "ssuextract",
            "viral",
            "checkm2",
            "tiara",
        )
    }


def _write_reports(
    directories: dict[str, Path],
    quickclade: list[dict[str, str]],
    checkeuk: list[dict[str, str]],
    gvclass: list[dict[str, str]],
    ssu: list[dict[str, str]],
    viruses: list[dict[str, str]],
    voting_markers: str = "74",
    *,
    tiara: tuple[dict[str, str], ...] = (),
    checkm2: tuple[dict[str, str], ...] | None = None,
    ssu_lineage_columns: bool = True,
):
    ssu_fields = (
        "name",
        "model",
        "contig_name",
        "taxonomy",
        "taxonomy_domain",
        "reference_taxonomy",
        "blast_candidate_taxonomy",
    )
    if not ssu_lineage_columns:
        ssu_fields = ssu_fields[:5]
        ssu = [{field: row[field] for field in ssu_fields} for row in ssu]
    _write_tsv(
        directories["quickclade"] / "quickclade.tsv",
        ("#QueryName", "Q_Bases", "lineage", "ConfLevel", "Confidence"),
        quickclade,
    )
    _write_tsv(
        directories["checkeuk"] / "result_checkeuk/checkeuk_report.tsv",
        ("genome", "status", "lineage"),
        checkeuk,
    )
    _write_tsv(
        directories["checkeuk"] / "result_checkeuk/checkeuk_report_detail.tsv",
        ("genome", "marker_tax_host_lineage_selected_markers"),
        [
            {
                "genome": row["genome"],
                "marker_tax_host_lineage_selected_markers": voting_markers,
            }
            for row in checkeuk
        ],
    )
    _write_tsv(
        directories["gvclass"] / "output/gvclass_summary.tsv",
        ("query", "taxonomy_majority", "taxonomy_confidence", "domain"),
        gvclass,
    )
    _write_tsv(
        directories["ssuextract"] / "output/cmsearch_summary.tsv",
        ssu_fields,
        ssu,
    )
    _write_tsv(
        directories["viral"] / "genomad/input_summary/input_virus_summary.tsv",
        ("seq_name", "virus_score", "taxonomy"),
        viruses,
    )
    _write_tsv(
        directories["tiara"] / "tiara.tsv",
        ("sequence_id", "class_fst_stage", "class_snd_stage"),
        list(tiara),
    )
    if checkm2 is not None:
        _write_tsv(
            directories["checkm2"] / "output/quality_report.tsv",
            ("Name", "Completeness", "Contamination", "Completeness_Model_Used"),
            list(checkm2),
        )


def _read_evidence(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["bin_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _sample(path: Path):
    path.write_text(
        json.dumps({"sample_id": "sample-a", "run_id": "run-a"}), encoding="utf-8"
    )


def test_route_bins_uses_concordant_and_physically_linked_evidence(tmp_path):
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    bins = {
        "euk_gv": "contig_euk_gv",
        "euk_ssu": "contig_euk_ssu",
        "shared_taxon": "contig_shared_taxon",
        "prok": "contig_prok",
        "archaeal": "contig_archaeal",
        "archaeal_euk_conflict": "contig_archaeal_euk_conflict",
        "viral": "contig_viral",
        "euk_provirus": "contig_euk_provirus",
        "prok_provirus": "contig_prok_provirus",
        "provirus_only": "contig_provirus_only",
        "prok_viral_conflict": "contig_prok_viral",
        "unclassified": "contig_unclassified",
    }
    for bin_id, contig_id in bins.items():
        (bin_dir / f"{bin_id}.fa").write_text(
            f">{contig_id} full description\nACGTACGT\n", encoding="utf-8"
        )

    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "archaeal.fa",
                "Q_Bases": "8",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            },
            {
                "#QueryName": "archaeal_euk_conflict.fa",
                "Q_Bases": "8",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            },
            {
                "#QueryName": "euk_gv.fa",
                "Q_Bases": "8",
                "lineage": "d__Bacteria;k__Bacillati",
                "ConfLevel": "kingdom:98.4",
                "Confidence": "d:100;k:98.4",
            },
            {
                "#QueryName": "euk_ssu.fa",
                "Q_Bases": "8",
                "lineage": "d__Eukaryota;k__Chromista",
                "ConfLevel": "domain:99.0",
                "Confidence": "d:99.0",
            },
            {
                "#QueryName": "shared_taxon.fa",
                "Q_Bases": "8",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            },
            {
                "#QueryName": "prok.fa",
                "Q_Bases": "8",
                "lineage": "d__Archaea;k__Thermoproteati",
                "ConfLevel": "domain:100",
                "Confidence": "d:100",
            },
            {
                "#QueryName": "viral.fa",
                "Q_Bases": "8",
                "lineage": "k__Heunggongvirae",
                "ConfLevel": "kingdom:94",
                "Confidence": "k:94",
            },
            {
                "#QueryName": "euk_provirus.fa",
                "Q_Bases": "8",
                "lineage": "d__Eukaryota;k__Chromista",
                "ConfLevel": "domain:99",
                "Confidence": "d:99",
            },
            {
                "#QueryName": "prok_provirus.fa",
                "Q_Bases": "8",
                "lineage": "d__Bacteria;k__Pseudomonadati",
                "ConfLevel": "domain:99",
                "Confidence": "d:99",
            },
            {
                "#QueryName": "provirus_only.fa",
                "Q_Bases": "8",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            },
            {
                "#QueryName": "prok_viral_conflict.fa",
                "Q_Bases": "8",
                "lineage": "d__Bacteria;k__Pseudomonadati",
                "ConfLevel": "domain:99",
                "Confidence": "d:99",
            },
            {
                "#QueryName": "unclassified.fa",
                "Q_Bases": "8",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            },
        ],
        checkeuk=[
            {"genome": "archaeal", "status": "filtered", "lineage": ""},
            {"genome": "archaeal_euk_conflict", "status": "ok", "lineage": "sg__Sar"},
            {"genome": "euk_gv", "status": "ok", "lineage": "sg__Sar"},
            {"genome": "euk_ssu", "status": "ok", "lineage": "sg__Sar"},
            {"genome": "shared_taxon", "status": "ok", "lineage": "sg__Sar"},
            {"genome": "prok", "status": "filtered", "lineage": ""},
            {"genome": "viral", "status": "filtered", "lineage": ""},
            {
                "genome": "euk_provirus",
                "status": "ok",
                "lineage": "sg__Sar",
            },
            {
                "genome": "prok_provirus",
                "status": "filtered",
                "lineage": "",
            },
            {
                "genome": "provirus_only",
                "status": "filtered",
                "lineage": "",
            },
            {
                "genome": "prok_viral_conflict",
                "status": "filtered",
                "lineage": "",
            },
            {"genome": "unclassified", "status": "filtered", "lineage": ""},
        ],
        gvclass=[
            {
                "query": "archaeal",
                "taxonomy_majority": "d_ARC;p_Nanobdellota",
                "taxonomy_confidence": "low_support",
                "domain": "",
            },
            {
                "query": "archaeal_euk_conflict",
                "taxonomy_majority": "d_ARC;p_Nanobdellota",
                "taxonomy_confidence": "low_support",
                "domain": "",
            },
            {
                "query": "euk_gv",
                "taxonomy_majority": "d_EUK;p_Ochrophyta",
                "taxonomy_confidence": "high",
                "domain": "EUK:12(100%)",
            },
            {
                "query": "euk_ssu",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            },
            {
                "query": "shared_taxon",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            },
            {
                "query": "prok",
                "taxonomy_majority": "d_BAC;p_Pseudomonadota",
                "taxonomy_confidence": "high",
                "domain": "BAC:9(100%)",
            },
            {
                "query": "viral",
                "taxonomy_majority": "d_NCLDV;p_Nucleocytoviricota",
                "taxonomy_confidence": "high",
                "domain": "NCLDV:8(100%)",
            },
            {
                "query": "euk_provirus",
                "taxonomy_majority": "d_EUK;p_Ochrophyta",
                "taxonomy_confidence": "high",
                "domain": "EUK:12(100%)",
            },
            {
                "query": "prok_provirus",
                "taxonomy_majority": "d_BAC;p_Pseudomonadota",
                "taxonomy_confidence": "high",
                "domain": "BAC:9(100%)",
            },
            {
                "query": "provirus_only",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            },
            {
                "query": "prok_viral_conflict",
                "taxonomy_majority": "d_BAC;p_Pseudomonadota",
                "taxonomy_confidence": "high",
                "domain": "BAC:9(100%)",
            },
            {
                "query": "unclassified",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            },
        ],
        ssu=[
            {
                "name": "archaeal-euk18s-linked",
                "model": "RF01960",
                "contig_name": "contig_archaeal_euk_conflict",
                "taxonomy": "Eukaryota;Sar",
                "taxonomy_domain": "Eukaryota",
            },
            {
                "name": "euk18s-linked",
                "model": "RF01960",
                "contig_name": "contig_euk_ssu",
                "taxonomy": "Eukaryota;Sar;Ochrophyta",
                "taxonomy_domain": "Eukaryota",
            },
            {
                "name": "same-taxon-other-contig",
                "model": "RF01960",
                "contig_name": "contig_euk_ssu",
                "taxonomy": "Eukaryota;Sar",
                "taxonomy_domain": "Eukaryota",
            },
        ],
        viruses=[
            {
                "seq_name": "contig_euk_provirus|provirus_2_7",
                "virus_score": "0.98",
                "taxonomy": "Viruses;Duplodnaviria",
            },
            {
                "seq_name": "contig_prok_provirus|provirus_2_7",
                "virus_score": "0.97",
                "taxonomy": "Viruses;Duplodnaviria",
            },
            {
                "seq_name": "contig_provirus_only|provirus_2_7",
                "virus_score": "0.96",
                "taxonomy": "Viruses;Duplodnaviria",
            },
            {
                "seq_name": "contig_prok_viral",
                "virus_score": "0.99",
                "taxonomy": "Viruses;Varidnaviria",
            },
        ],
    )

    output = tmp_path / "routing"
    assert route_bins(sample_path, bin_dir, directories, output) == output
    evidence = _read_evidence(output / "evidence.tsv")

    assert evidence["euk_gv"]["route"] == "eukaryotic"
    assert evidence["euk_gv"]["candidate_class"] == "eukaryotic_candidate"
    assert json.loads(evidence["euk_gv"]["conflicting_evidence"])
    assert evidence["euk_ssu"]["route"] == "eukaryotic"
    assert len(json.loads(evidence["euk_ssu"]["linked_eukaryotic_ssu"])) == 2
    assert evidence["shared_taxon"]["route"] == "unresolved"
    assert json.loads(evidence["shared_taxon"]["linked_eukaryotic_ssu"]) == []
    assert evidence["prok"]["candidate_class"] == "bacterial_archaeal_candidate"
    assert evidence["archaeal"]["route"] == "prokaryotic"
    assert evidence["archaeal"]["candidate_class"] == "bacterial_archaeal_candidate"
    assert evidence["archaeal"]["gvclass_domain"] == "d_ARC"
    assert evidence["archaeal"]["gvclass_taxonomy_confidence"] == "low_support"
    assert evidence["archaeal"]["reason"] == (
        "GVClass supports a prokaryotic route without a "
        "supported eukaryotic or viral conflict"
    )
    assert json.loads(evidence["archaeal"]["supporting_evidence"]) == [
        "GVClass:d_ARC;p_Nanobdellota"
    ]
    assert json.loads(evidence["archaeal"]["conflicting_evidence"]) == []
    assert evidence["archaeal_euk_conflict"]["route"] == "unresolved"
    assert evidence["archaeal_euk_conflict"]["candidate_class"] == "conflicting"
    assert json.loads(evidence["archaeal_euk_conflict"]["supporting_evidence"]) == [
        "CheckEUK:sg__Sar",
        "SSU:archaeal-euk18s-linked@contig_archaeal_euk_conflict",
        "GVClass:d_ARC;p_Nanobdellota",
    ]
    assert json.loads(evidence["archaeal_euk_conflict"]["conflicting_evidence"]) == [
        "GVClass:d_ARC;p_Nanobdellota"
    ]
    assert evidence["viral"]["candidate_class"] == "viral_candidate"
    assert evidence["euk_provirus"]["candidate_class"] == "eukaryotic_candidate"
    assert len(json.loads(evidence["euk_provirus"]["linked_genomad_viruses"])) == 1
    assert json.loads(evidence["euk_provirus"]["supporting_evidence"]) == [
        "CheckEUK:sg__Sar",
        "GVClass:d_EUK;p_Ochrophyta",
    ]
    assert json.loads(evidence["euk_provirus"]["conflicting_evidence"]) == []
    assert evidence["prok_provirus"]["candidate_class"] == (
        "bacterial_archaeal_candidate"
    )
    assert len(json.loads(evidence["prok_provirus"]["linked_genomad_viruses"])) == 1
    assert json.loads(evidence["prok_provirus"]["supporting_evidence"]) == [
        "GVClass:d_BAC;p_Pseudomonadota",
        "QuickClade:d__Bacteria;k__Pseudomonadati",
    ]
    assert json.loads(evidence["prok_provirus"]["conflicting_evidence"]) == []
    assert evidence["provirus_only"]["candidate_class"] == "unresolved"
    assert len(json.loads(evidence["provirus_only"]["linked_genomad_viruses"])) == 1
    assert json.loads(evidence["provirus_only"]["supporting_evidence"]) == []
    assert json.loads(evidence["provirus_only"]["conflicting_evidence"]) == []
    assert evidence["prok_viral_conflict"]["candidate_class"] == "conflicting"
    assert "geNomad:contig_prok_viral" in json.loads(
        evidence["prok_viral_conflict"]["supporting_evidence"]
    )
    assert json.loads(evidence["prok_viral_conflict"]["conflicting_evidence"]) == [
        "prokaryotic evidence conflicts with viral evidence"
    ]
    assert evidence["unclassified"]["candidate_class"] == "unresolved"
    assert set(evidence) == set(bins)

    copies = list(output.glob("*/*.fna"))
    assert {path.stem for path in copies} == set(bins)
    assert len(copies) == len(bins)
    for path in copies:
        assert path.read_bytes() == (bin_dir / f"{path.stem}.fa").read_bytes()
    manifest = StageManifest.model_validate_json(
        (output / "stage.json").read_text(encoding="utf-8")
    )
    assert manifest.stage == "routing"
    assert manifest.status == "completed"
    assert manifest.outputs["evidence"] == "evidence.tsv"


@pytest.mark.parametrize("domain", ["d_PLASTID", "d_MITO"])
@pytest.mark.parametrize(
    ("quickclade_lineage", "native_support", "expected_support", "expected"),
    [
        ("", ("", [], []), [], ("unresolved", [])),
        (
            "d__Bacteria",
            ("", [], []),
            [],
            ("unresolved", []),
        ),
        (
            "",
            (
                "sg__Sar",
                [
                    {
                        "name": "linked18s",
                        "model": "RF01960",
                        "contig_name": "contig_1",
                        "taxonomy": "Eukaryota;Sar",
                        "taxonomy_domain": "Eukaryota",
                    }
                ],
                [],
            ),
            ["CheckEUK:sg__Sar", "SSU:linked18s@contig_1"],
            (
                "conflicting",
                ["organelle evidence conflicts with supported routing evidence"],
            ),
        ),
        (
            "",
            (
                "",
                [],
                [
                    {
                        "seq_name": "contig_1",
                        "virus_score": "0.99",
                        "taxonomy": "Viruses;Duplodnaviria",
                    }
                ],
            ),
            ["geNomad:contig_1"],
            (
                "conflicting",
                ["organelle evidence conflicts with supported routing evidence"],
            ),
        ),
    ],
    ids=["alone", "quickclade-only", "eukaryotic", "viral"],
)
def test_route_bins_retains_organelle_evidence_without_promoting_a_route(
    tmp_path: Path,
    domain: str,
    quickclade_lineage: str,
    native_support: tuple[str, list[dict[str, str]], list[dict[str, str]]],
    expected_support: list[str],
    expected: tuple[str, list[str]],
) -> None:
    """Keep each organelle unresolved with independently supported conflicts."""
    checkeuk_lineage, ssu, viruses = native_support
    expected_class, expected_conflicts = expected
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "4",
                "lineage": quickclade_lineage,
                "ConfLevel": "domain:100",
                "Confidence": "d:100",
            }
        ],
        checkeuk=[{"genome": "bin_1", "status": "ok", "lineage": checkeuk_lineage}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": f"{domain};p_example",
                "taxonomy_confidence": "low_support",
                "domain": "",
            }
        ],
        ssu=ssu,
        viruses=viruses,
    )

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")
    evidence = _read_evidence(output / "evidence.tsv")["bin_1"]

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == expected_class
    assert evidence["reason"] == (
        "GVClass reports an organelle class without a gene-calling route"
    )
    assert evidence["gvclass_domain"] == domain
    assert evidence["gvclass_taxonomy_majority"] == f"{domain};p_example"
    assert evidence["gvclass_taxonomy_confidence"] == "low_support"
    assert evidence["checkeuk_lineage"] == checkeuk_lineage
    assert evidence["checkeuk_status"] == "ok"
    assert evidence["quickclade_lineage"] == quickclade_lineage
    assert evidence["quickclade_confidence"] == "d:100"
    assert json.loads(evidence["linked_genomad_viruses"]) == viruses
    assert json.loads(evidence["supporting_evidence"]) == [
        f"GVClass:{domain};p_example",
        *expected_support,
    ]
    assert json.loads(evidence["conflicting_evidence"]) == expected_conflicts
    assert list(output.glob("*/*.fna")) == [output / "unresolved/bin_1.fna"]


@pytest.mark.parametrize(
    ("quickclade_lineage", "checkeuk_lineage", "expected"),
    [
        ("", "", ("viral", "viral_candidate", [])),
        (
            "d__Bacteria",
            "",
            ("viral", "viral_candidate", []),
        ),
        (
            "",
            "sg__Sar",
            (
                "unresolved",
                "conflicting",
                ["eukaryotic evidence conflicts with viral evidence"],
            ),
        ),
    ],
    ids=["alone", "quickclade-only", "eukaryotic-conflict"],
)
def test_route_bins_applies_viral_rules_to_phage(
    tmp_path: Path,
    quickclade_lineage: str,
    checkeuk_lineage: str,
    expected: tuple[str, str, list[str]],
) -> None:
    """Retain low-support phage calls under the existing viral policy."""
    expected_route, expected_class, expected_conflicts = expected
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "4",
                "lineage": quickclade_lineage,
                "ConfLevel": "domain:100",
                "Confidence": "d:100",
            }
        ],
        checkeuk=[{"genome": "bin_1", "status": "ok", "lineage": checkeuk_lineage}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": "d_PHAGE;p_example",
                "taxonomy_confidence": "low_support",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=[],
    )

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")
    evidence = _read_evidence(output / "evidence.tsv")["bin_1"]

    assert evidence["route"] == expected_route
    assert evidence["candidate_class"] == expected_class
    assert evidence["gvclass_domain"] == "d_PHAGE"
    assert evidence["gvclass_taxonomy_majority"] == "d_PHAGE;p_example"
    assert evidence["quickclade_lineage"] == quickclade_lineage
    assert evidence["gvclass_taxonomy_confidence"] == "low_support"
    assert "GVClass:d_PHAGE;p_example" in json.loads(evidence["supporting_evidence"])
    assert json.loads(evidence["conflicting_evidence"]) == expected_conflicts
    assert list(output.glob("*/*.fna")) == [output / expected_route / "bin_1.fna"]


def test_route_bins_accepts_an_empty_bin_set_without_native_tables(
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    directories = _evidence_dirs(tmp_path / "native")

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")

    assert _read_evidence(output / "evidence.tsv") == {}
    with (output / "evidence.tsv").open(encoding="utf-8", newline="") as handle:
        assert tuple(csv.DictReader(handle, delimiter="\t").fieldnames or ()) == (
            "sample_id",
            "run_id",
            "bin_id",
            "route",
            "candidate_class",
            "reason",
            "checkeuk_status",
            "checkeuk_lineage",
            "taxonomy",
            "taxonomy_source",
            "quickclade_domain",
            "quickclade_lineage",
            "quickclade_confidence",
            "gvclass_domain",
            "gvclass_taxonomy_majority",
            "gvclass_taxonomy_confidence",
            "linked_eukaryotic_ssu",
            "linked_genomad_viruses",
            "supporting_evidence",
            "conflicting_evidence",
            "bin_bp",
            "tiara_eukaryote_fraction",
            "tiara_prokaryote_fraction",
            "tiara_organelle_fraction",
            "tiara_unclassified_fraction",
            "viral_contig_fraction",
            "checkm2_completeness",
            "checkm2_contamination",
            "duplicated_contig_fraction",
        )
    assert all(
        (output / route).is_dir()
        for route in (
            "prokaryotic",
            "eukaryotic",
            "viral",
            "unresolved",
        )
    )
    assert not list(output.glob("*/*.fna"))
    manifest = StageManifest.model_validate_json(
        (output / "stage.json").read_text(encoding="utf-8")
    )
    assert manifest.status == "completed"
    assert set(manifest.outputs) == {
        "evidence",
        "prokaryotic",
        "eukaryotic",
        "viral",
        "unresolved",
    }


def test_route_bins_records_an_unrouted_quickclade_domain(
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "4",
                "lineage": "sk__Viruses;k__Orthornavirae",
                "ConfLevel": "superkingdom:100",
                "Confidence": "sk:100",
            }
        ],
        checkeuk=[{"genome": "bin_1", "status": "filtered", "lineage": ""}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=[],
    )

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")
    evidence = _read_evidence(output / "evidence.tsv")["bin_1"]

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "unresolved"
    assert evidence["quickclade_domain"] == "sk__Viruses"
    assert json.loads(evidence["conflicting_evidence"]) == [
        "QuickClade-unresolved:sk__Viruses;k__Orthornavirae"
    ]
    assert (output / "unresolved/bin_1.fna").is_file()


def test_route_bins_rejects_a_malformed_quickclade_domain(
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "4",
                "lineage": "sk__",
                "ConfLevel": "superkingdom:100",
                "Confidence": "sk:100",
            }
        ],
        checkeuk=[{"genome": "bin_1", "status": "filtered", "lineage": ""}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=[],
    )

    with pytest.raises(ValueError, match="malformed QuickClade domain"):
        route_bins(sample_path, bin_dir, directories, tmp_path / "routing")


def test_route_bins_rejects_an_unknown_gvclass_domain(tmp_path: Path) -> None:
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "4",
                "lineage": "d__Bacteria",
                "ConfLevel": "domain:100",
                "Confidence": "d:100",
            }
        ],
        checkeuk=[{"genome": "bin_1", "status": "filtered", "lineage": ""}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": "d_UNKNOWN;p_unknown",
                "taxonomy_confidence": "low",
                "domain": "UNKNOWN:1(100%)",
            }
        ],
        ssu=[],
        viruses=[],
    )

    with pytest.raises(ValueError, match="unsupported GVClass domain: d_UNKNOWN"):
        route_bins(sample_path, bin_dir, directories, tmp_path / "routing")


@pytest.mark.parametrize(
    ("gvclass_lineage", "reason", "viruses", "expected"),
    [
        (
            "d_EUK;p_example",
            "eukaryotic evidence does not meet the combined support rule",
            [],
            ("unresolved", "unresolved", [], []),
        ),
        (
            "d_",
            "no supported routing evidence",
            [],
            ("unresolved", "unresolved", [], []),
        ),
        (
            "d_",
            "GVClass or linked whole-contig geNomad evidence supports a viral route",
            [{"seq_name": "contig_3", "virus_score": "0.99", "taxonomy": "Viruses"}],
            ("viral", "viral_candidate", ["geNomad:contig_3"], []),
        ),
        (
            "d_EUK;p_example",
            "supported evidence conflicts across routes",
            [{"seq_name": "contig_3", "virus_score": "0.99", "taxonomy": "Viruses"}],
            (
                "unresolved",
                "conflicting",
                ["geNomad:contig_3"],
                ["eukaryotic evidence conflicts with viral evidence"],
            ),
        ),
    ],
)
def test_route_bins_ignores_quickclade_routing(
    tmp_path: Path,
    gvclass_lineage: str,
    reason: str,
    viruses: list[dict[str, str]],
    expected: tuple[str, str, list[str], list[str]],
) -> None:
    """Route bin_3-like blank lineages; Pfam 69.936/BOM 65.25 are not taxonomy."""
    route, candidate_class, supporting, conflicting = expected
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_3.fa").write_text(">contig_3\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_3.fa",
                "Q_Bases": "4",
                "lineage": "d__Bacteria;k__Bacillati",
                "ConfLevel": "kingdom:97.7",
                "Confidence": "d:100;k:97.7",
            }
        ],
        checkeuk=[{"genome": "bin_3", "status": "ok", "lineage": ""}],
        gvclass=[
            {
                "query": "bin_3",
                "taxonomy_majority": gvclass_lineage,
                "taxonomy_confidence": "low_support",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=viruses,
    )

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")
    evidence = _read_evidence(output / "evidence.tsv")["bin_3"]

    assert evidence["route"] == route
    assert evidence["candidate_class"] == candidate_class
    assert evidence["reason"] == reason
    assert evidence["quickclade_lineage"] == "d__Bacteria;k__Bacillati"
    assert evidence["quickclade_confidence"] == "d:100;k:97.7"
    assert evidence["checkeuk_status"] == "ok"
    assert evidence["checkeuk_lineage"] == ""
    assert evidence["gvclass_taxonomy_majority"] == gvclass_lineage
    assert json.loads(evidence["supporting_evidence"]) == supporting
    assert json.loads(evidence["conflicting_evidence"]) == conflicting
    assert list(output.glob("*/*.fna")) == [output / route / "bin_3.fna"]
    assert (output / route / "bin_3.fna").read_bytes() == (
        bin_dir / "bin_3.fa"
    ).read_bytes()


@pytest.mark.parametrize(
    ("voting_markers", "gvclass_lineage", "expected"),
    [
        (
            "4",
            "d_EUK;p_example",
            (
                "eukaryotic",
                "gvclass_fallback",
                "d_EUK;p_example",
                FALLBACK_REASON,
            ),
        ),
        (
            "",
            "d_EUK;p_example",
            (
                "eukaryotic",
                "gvclass_fallback",
                "d_EUK;p_example",
                FALLBACK_REASON,
            ),
        ),
        (
            "5",
            "d_EUK;p_example",
            (
                "unresolved",
                "",
                "",
                "eukaryotic evidence does not meet the combined support rule",
            ),
        ),
        (
            "4",
            "d_BAC;p_example",
            ("prokaryotic", "", "", None),
        ),
    ],
    ids=["few-markers", "empty-markers", "enough-markers", "prokaryotic"],
)
def test_route_bins_falls_back_to_gvclass_for_few_checkeuk_markers(
    tmp_path: Path,
    voting_markers: str,
    gvclass_lineage: str,
    expected: tuple[str, str, str, str | None],
) -> None:
    route, source, taxonomy, reason = expected
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "4",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            }
        ],
        checkeuk=[{"genome": "bin_1", "status": "ok", "lineage": ""}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": gvclass_lineage,
                "taxonomy_confidence": "low_support",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=[],
        voting_markers=voting_markers,
    )

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")
    evidence = _read_evidence(output / "evidence.tsv")["bin_1"]

    assert evidence["route"] == route
    assert evidence["taxonomy_source"] == source
    assert evidence["taxonomy"] == taxonomy
    if reason is not None:
        assert evidence["reason"] == reason


def test_route_bins_keeps_checkeuk_taxonomy_when_a_lineage_exists(
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "4",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            }
        ],
        checkeuk=[
            {"genome": "bin_1", "status": "ok", "lineage": "sg__SAR;p__Ochrophyta"}
        ],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": "d_EUK;p_other",
                "taxonomy_confidence": "low_support",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=[],
        voting_markers="2",
    )

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")
    evidence = _read_evidence(output / "evidence.tsv")["bin_1"]

    assert evidence["route"] == "eukaryotic"
    assert evidence["reason"] == "CheckEUK and GVClass support a eukaryotic route"
    assert evidence["taxonomy_source"] == "checkeuk"
    assert evidence["taxonomy"] == "sg__SAR;p__Ochrophyta"


EUK_18S = {
    "name": "locus",
    "model": "RF01960",
    "contig_name": "contig_1",
    "taxonomy": "Unclassified",
    "taxonomy_domain": "Unclassified",
    "reference_taxonomy": "Eukaryota;Amoebozoa",
    "blast_candidate_taxonomy": "Unclassified",
}
WHOLE_VIRUS = {"seq_name": "contig_1", "virus_score": "0.99", "taxonomy": "Viruses"}


def _tiara(contig: str, first: str, second: str = "n/a") -> dict[str, str]:
    return {
        "sequence_id": contig,
        "class_fst_stage": first,
        "class_snd_stage": second,
    }


def _checkm2(completeness: str, contamination: str) -> tuple[dict[str, str], ...]:
    return (
        {
            "Name": "bin_1",
            "Completeness": completeness,
            "Contamination": contamination,
            "Completeness_Model_Used": "Neural Network (Specific Model)",
        },
    )


def _route_one(
    tmp_path: Path,
    *,
    contigs: dict[str, int] | None = None,
    checkeuk_lineage: str = "",
    gvclass: str = "d_",
    ssu: tuple[dict[str, str], ...] = (),
    viruses: tuple[dict[str, str], ...] = (),
    tiara: tuple[dict[str, str], ...] = (),
    checkm2: tuple[dict[str, str], ...] | None = None,
    voting_markers: str = "74",
    ssu_lineage_columns: bool = True,
) -> dict[str, str]:
    """Route one bin named bin_1 and return its evidence row."""
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(
        "".join(
            f">{name} description\n{'A' * length}\n"
            for name, length in (contigs or {"contig_1": 100}).items()
        ),
        encoding="utf-8",
    )
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {
                "#QueryName": "bin_1.fa",
                "Q_Bases": "100",
                "lineage": "",
                "ConfLevel": "",
                "Confidence": "",
            }
        ],
        checkeuk=[{"genome": "bin_1", "status": "ok", "lineage": checkeuk_lineage}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": gvclass,
                "taxonomy_confidence": "high",
                "domain": "",
            }
        ],
        ssu=list(ssu),
        viruses=list(viruses),
        voting_markers=voting_markers,
        tiara=tiara,
        checkm2=checkm2,
        ssu_lineage_columns=ssu_lineage_columns,
    )
    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")
    return _read_evidence(output / "evidence.tsv")["bin_1"]


@pytest.mark.parametrize(
    "locus",
    [
        {},
        {"reference_taxonomy": "", "blast_candidate_taxonomy": "Eukaryota;Sar"},
        {
            "taxonomy_domain": "Eukaryota",
            "reference_taxonomy": "",
        },
    ],
    ids=["reference", "blast-candidate", "native-domain"],
)
def test_route_bins_links_an_18s_locus_when_ssuextract_leaves_the_domain_unclassified(
    tmp_path: Path, locus: dict[str, str]
) -> None:
    """SSUextract 1.2.2 writes Unclassified in taxonomy_domain for every locus."""
    evidence = _route_one(
        tmp_path, checkeuk_lineage="sg__Sar", ssu=({**EUK_18S, **locus},)
    )

    assert evidence["route"] == "eukaryotic"
    assert evidence["reason"] == (
        "CheckEUK and a physically linked eukaryotic 18S locus support "
        "a eukaryotic route"
    )
    assert json.loads(evidence["supporting_evidence"]) == [
        "CheckEUK:sg__Sar",
        "SSU:locus@contig_1",
    ]
    assert len(json.loads(evidence["linked_eukaryotic_ssu"])) == 1


def test_route_bins_links_a_locus_from_an_ssu_table_without_lineage_columns(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        checkeuk_lineage="sg__Sar",
        ssu=({**EUK_18S, "taxonomy_domain": "Eukaryota"},),
        ssu_lineage_columns=False,
    )

    assert evidence["route"] == "eukaryotic"
    assert len(json.loads(evidence["linked_eukaryotic_ssu"])) == 1


def test_route_bins_accepts_a_header_only_ssu_table_without_lineage_columns(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path, checkeuk_lineage="sg__Sar", ssu=(), ssu_lineage_columns=False
    )

    assert json.loads(evidence["linked_eukaryotic_ssu"]) == []


@pytest.mark.parametrize(
    "locus",
    [
        {"reference_taxonomy": "Bacteria;Pseudomonadota"},
        {"model": "RF00177"},
        {"contig_name": "contig_elsewhere"},
    ],
    ids=["bacterial", "16s-model", "other-contig"],
)
def test_route_bins_does_not_link_a_locus_without_eukaryotic_18s_evidence(
    tmp_path: Path, locus: dict[str, str]
) -> None:
    evidence = _route_one(
        tmp_path, checkeuk_lineage="sg__Sar", ssu=({**EUK_18S, **locus},)
    )

    assert evidence["route"] == "unresolved"
    assert json.loads(evidence["linked_eukaryotic_ssu"]) == []


def test_route_bins_ignores_a_minority_viral_contig(tmp_path: Path) -> None:
    """One small whole-contig geNomad call does not veto a eukaryotic route."""
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 960, "contig_2": 40},
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        viruses=({**WHOLE_VIRUS, "seq_name": "contig_2"},),
    )

    assert evidence["route"] == "eukaryotic"
    assert evidence["candidate_class"] == "eukaryotic_candidate"
    assert evidence["viral_contig_fraction"] == "0.0400"
    assert json.loads(evidence["supporting_evidence"]) == [
        "CheckEUK:sg__Sar",
        "GVClass:d_EUK;p_example",
    ]
    assert json.loads(evidence["conflicting_evidence"]) == ["geNomad-minority:0.0400"]


def test_route_bins_lets_a_majority_viral_contig_support_a_viral_route(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 600, "contig_2": 400},
        viruses=(WHOLE_VIRUS,),
    )

    assert evidence["route"] == "viral"
    assert evidence["viral_contig_fraction"] == "0.6000"
    assert json.loads(evidence["supporting_evidence"]) == ["geNomad:contig_1"]
    assert json.loads(evidence["conflicting_evidence"]) == []


def test_route_bins_counts_only_whole_contig_calls_toward_the_viral_fraction(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        viruses=({**WHOLE_VIRUS, "seq_name": "contig_1|provirus_5_80"},),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["viral_contig_fraction"] == "0.0000"


def test_route_bins_routes_a_complete_bacterial_genome_despite_eukaryote_labels(
    tmp_path: Path,
) -> None:
    """CheckEUK names a lineage and GVClass says d_EUK for a single-contig genome."""
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 4000},
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("98.5", "1.2"),
    )

    assert evidence["route"] == "prokaryotic"
    assert evidence["candidate_class"] == "bacterial_archaeal_candidate"
    assert evidence["reason"] == "Tiara and CheckM2 support a prokaryotic route"
    assert json.loads(evidence["supporting_evidence"]) == [
        "Tiara:1.0000",
        "CheckM2:98.50/1.20",
    ]
    assert json.loads(evidence["conflicting_evidence"]) == [
        "CheckEUK:sg__Sar",
        "GVClass:d_EUK;p_example",
    ]
    assert evidence["bin_bp"] == "4000"
    assert evidence["tiara_prokaryote_fraction"] == "1.0000"
    assert evidence["checkm2_completeness"] == "98.50"
    assert evidence["checkm2_contamination"] == "1.20"


@pytest.mark.parametrize(
    ("tiara", "checkm2", "reason"),
    [
        ((_tiara("contig_1", "bacteria"),), None, "no supported routing evidence"),
        (
            (_tiara("contig_1", "bacteria"),),
            _checkm2("9.9", "1"),
            "no supported routing evidence",
        ),
        ((), _checkm2("99", "0.5"), "no supported routing evidence"),
        (
            (_tiara("contig_1", "eukarya"),),
            _checkm2("99", "0.5"),
            "eukaryotic evidence does not meet the combined support rule",
        ),
    ],
    ids=[
        "no-checkm2",
        "incomplete",
        "no-tiara-line",
        "tiara-eukaryote",
    ],
)
def test_route_bins_needs_both_tiara_and_checkm2_for_the_pair_route(
    tmp_path: Path,
    tiara: tuple[dict[str, str], ...],
    checkm2: tuple[dict[str, str], ...] | None,
    reason: str,
) -> None:
    evidence = _route_one(tmp_path, tiara=tiara, checkm2=checkm2)

    assert evidence["route"] == "unresolved"
    assert evidence["reason"] == reason


@pytest.mark.parametrize("contamination", ["0", "10", "60"])
def test_route_bins_accepts_ten_percent_checkm2_completeness_at_any_contamination(
    tmp_path: Path, contamination: str
) -> None:
    evidence = _route_one(
        tmp_path,
        tiara=(_tiara("contig_1", "archaea"),),
        checkm2=_checkm2("10", contamination),
    )

    assert evidence["route"] == "prokaryotic"
    assert evidence["reason"] == "Tiara and CheckM2 support a prokaryotic route"


@pytest.mark.parametrize(
    ("prokaryote_bp", "route"),
    [(800, "prokaryotic"), (799, "unresolved")],
    ids=["at-threshold", "below-threshold"],
)
def test_route_bins_needs_a_tiara_prokaryote_share_of_eighty_percent(
    tmp_path: Path, prokaryote_bp: int, route: str
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": prokaryote_bp, "contig_2": 1000 - prokaryote_bp},
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("95", "2"),
    )

    assert evidence["route"] == route


def test_route_bins_prefers_the_gvclass_reason_and_lists_every_prokaryote_item(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        gvclass="d_BAC;p_Pseudomonadota",
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("40", "3"),
    )

    assert evidence["route"] == "prokaryotic"
    assert evidence["reason"] == (
        "GVClass supports a prokaryotic route without a supported eukaryotic "
        "or viral conflict"
    )
    assert json.loads(evidence["supporting_evidence"]) == [
        "GVClass:d_BAC;p_Pseudomonadota",
        "Tiara:1.0000",
        "CheckM2:40.00/3.00",
    ]


def test_route_bins_routes_a_partial_bacterial_genome_with_false_eukaryote_labels(
    tmp_path: Path,
) -> None:
    """CheckM2 at 27 percent still sends the bin to GTDB-Tk."""
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 4000},
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("27", "40"),
    )

    assert evidence["route"] == "prokaryotic"
    assert evidence["reason"] == "Tiara and CheckM2 support a prokaryotic route"
    assert json.loads(evidence["conflicting_evidence"]) == [
        "CheckEUK:sg__Sar",
        "GVClass:d_EUK;p_example",
    ]


def test_route_bins_conflicts_for_bacterial_tiara_calls_with_five_percent_checkm2(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 4000},
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("5", "40"),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "conflicting"
    assert evidence["reason"] == "supported evidence conflicts across routes"
    assert json.loads(evidence["conflicting_evidence"]) == ["Tiara:1.0000"]


def test_route_bins_does_not_veto_eukaryotes_for_a_partial_tiara_prokaryote_share(
    tmp_path: Path,
) -> None:
    """A ciliate that Tiara partly mislabels as bacterial stays eukaryotic."""
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 520, "contig_2": 480},
        checkeuk_lineage="sg__Sar;p__Ciliophora",
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("0.04", "5"),
    )

    assert evidence["route"] == "eukaryotic"
    assert evidence["reason"] == "CheckEUK and GVClass support a eukaryotic route"
    assert evidence["tiara_prokaryote_fraction"] == "0.5200"
    assert json.loads(evidence["conflicting_evidence"]) == []


def test_route_bins_counts_organelle_sequence_toward_the_tiara_eukaryote_line(
    tmp_path: Path,
) -> None:
    """Ciliate nanochromosomes look like organelle DNA to Tiara."""
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 1_350_000, "contig_2": 425_000, "contig_3": 725_000},
        checkeuk_lineage="sg__Sar;p__Ciliophora",
        gvclass="d_EUK;p_example",
        tiara=(
            _tiara("contig_1", "organelle", "plastid"),
            _tiara("contig_2", "eukarya"),
        ),
    )

    assert evidence["route"] == "eukaryotic"
    assert evidence["reason"] == "CheckEUK and GVClass support a eukaryotic route"
    assert evidence["tiara_organelle_fraction"] == "0.5400"
    assert evidence["tiara_eukaryote_fraction"] == "0.1700"
    assert json.loads(evidence["supporting_evidence"]) == [
        "CheckEUK:sg__Sar;p__Ciliophora",
        "GVClass:d_EUK;p_example",
        "Tiara:0.7100",
    ]


def test_route_bins_applies_the_tiara_organelle_rule_up_to_two_megabases(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 2_000_000},
        tiara=(_tiara("contig_1", "organelle", "mitochondrion"),),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["reason"] == (
        "Tiara classifies most of the bin as organelle sequence "
        "without a gene-calling route"
    )


def test_route_bins_ignores_tiara_organelle_calls_above_two_megabases(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 2_000_001},
        tiara=(_tiara("contig_1", "organelle", "mitochondrion"),),
    )

    assert evidence["reason"] == (
        "eukaryotic evidence does not meet the combined support rule"
    )


def test_route_bins_ignores_tiara_organelle_calls_when_gvclass_is_prokaryotic(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 60_000},
        gvclass="d_BAC;p_Pseudomonadota",
        tiara=(_tiara("contig_1", "organelle", "plastid"),),
    )

    assert evidence["route"] == "prokaryotic"
    assert evidence["tiara_organelle_fraction"] == "1.0000"


def test_route_bins_withholds_a_route_from_a_small_all_organelle_bin(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 60_000},
        tiara=(_tiara("contig_1", "organelle", "unknown"),),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "unresolved"
    assert evidence["reason"] == (
        "Tiara classifies most of the bin as organelle sequence "
        "without a gene-calling route"
    )


def test_route_bins_withholds_the_pair_route_when_an_18s_locus_is_linked(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("95", "2"),
        ssu=(EUK_18S,),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["reason"] == (
        "eukaryotic evidence does not meet the combined support rule"
    )


def test_route_bins_conflicts_when_the_pair_route_meets_viral_support(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("95", "2"),
        viruses=(WHOLE_VIRUS,),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "conflicting"
    assert evidence["reason"] == "supported evidence conflicts across routes"
    assert json.loads(evidence["conflicting_evidence"]) == [
        "prokaryotic evidence conflicts with viral evidence"
    ]


def test_route_bins_leaves_a_viral_bin_holding_an_18s_locus_unresolved(
    tmp_path: Path,
) -> None:
    """Before the 18S fix the locus was invisible and this bin routed viral."""
    evidence = _route_one(
        tmp_path,
        gvclass="d_NCLDV;p_Nucleocytoviricota",
        viruses=(WHOLE_VIRUS,),
        ssu=(EUK_18S,),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "conflicting"
    assert json.loads(evidence["conflicting_evidence"]) == [
        "eukaryotic evidence conflicts with viral evidence"
    ]


def test_route_bins_keeps_a_giant_virus_bin_viral_when_tiara_calls_it_eukaryotic(
    tmp_path: Path,
) -> None:
    """Tiara has no virus class and labels giant-virus contigs eukarya."""
    evidence = _route_one(
        tmp_path,
        gvclass="d_NCLDV;p_Nucleocytoviricota",
        tiara=(_tiara("contig_1", "eukarya"),),
    )

    assert evidence["route"] == "viral"
    assert evidence["candidate_class"] == "viral_candidate"
    assert evidence["tiara_eukaryote_fraction"] == "1.0000"
    assert json.loads(evidence["supporting_evidence"]) == [
        "GVClass:d_NCLDV;p_Nucleocytoviricota"
    ]
    assert json.loads(evidence["conflicting_evidence"]) == []


def test_route_bins_conflicts_when_a_giant_virus_bin_also_has_a_checkeuk_lineage(
    tmp_path: Path,
) -> None:
    """CheckEUK and Tiara form two eukaryote lines, which meet viral support."""
    evidence = _route_one(
        tmp_path,
        checkeuk_lineage="sg__Sar",
        gvclass="d_NCLDV;p_Nucleocytoviricota",
        tiara=(_tiara("contig_1", "eukarya"),),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "conflicting"
    assert evidence["reason"] == "supported evidence conflicts across routes"
    assert json.loads(evidence["conflicting_evidence"]) == [
        "eukaryotic evidence conflicts with viral evidence"
    ]


def test_route_bins_keeps_a_giant_virus_bin_viral_below_ten_percent_checkm2(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        gvclass="d_NCLDV;p_Nucleocytoviricota",
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("9.9", "1"),
    )

    assert evidence["route"] == "viral"
    assert json.loads(evidence["conflicting_evidence"]) == []


def test_route_bins_conflicts_when_a_giant_virus_bin_meets_the_pair_route(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        gvclass="d_NCLDV;p_Nucleocytoviricota",
        tiara=(_tiara("contig_1", "bacteria"),),
        checkm2=_checkm2("49.9", "1"),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "conflicting"


def test_route_bins_matches_a_tiara_row_that_carries_the_full_myloasm_header(
    tmp_path: Path,
) -> None:
    contig = "u1ctg_len-100_circular-no_depth-2-2-2_duplicated-no"
    evidence = _route_one(
        tmp_path,
        contigs={contig: 100},
        tiara=(_tiara(f"{contig} mult=1.00", "eukarya"),),
    )

    assert evidence["tiara_eukaryote_fraction"] == "1.0000"
    assert evidence["tiara_unclassified_fraction"] == "0.0000"


def test_route_bins_leaves_the_checkm2_columns_blank_without_a_table(
    tmp_path: Path,
) -> None:
    evidence = _route_one(tmp_path, gvclass="d_BAC;p_example")

    assert evidence["route"] == "prokaryotic"
    assert evidence["checkm2_completeness"] == ""
    assert evidence["checkm2_contamination"] == ""


def test_route_bins_leaves_the_checkm2_columns_blank_for_a_bin_without_a_row(
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {"#QueryName": "bin_1.fa", "Q_Bases": "4", "lineage": "", "ConfLevel": ""}
        ],
        checkeuk=[{"genome": "bin_1", "status": "filtered", "lineage": ""}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=[],
        checkm2=(),
    )

    output = route_bins(sample_path, bin_dir, directories, tmp_path / "routing")

    evidence = _read_evidence(output / "evidence.tsv")["bin_1"]
    assert evidence["checkm2_completeness"] == ""
    assert evidence["checkm2_contamination"] == ""


def test_route_bins_counts_contigs_absent_from_tiara_as_unclassified(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 600, "contig_2": 400},
        tiara=(_tiara("contig_1", "eukarya"),),
    )

    assert evidence["bin_bp"] == "1000"
    assert evidence["tiara_eukaryote_fraction"] == "0.6000"
    assert evidence["tiara_prokaryote_fraction"] == "0.0000"
    assert evidence["tiara_organelle_fraction"] == "0.0000"
    assert evidence["tiara_unclassified_fraction"] == "0.4000"


def test_route_bins_counts_tiara_unknown_as_unclassified(tmp_path: Path) -> None:
    evidence = _route_one(tmp_path, tiara=(_tiara("contig_1", "unknown", "unknown"),))

    assert evidence["tiara_unclassified_fraction"] == "1.0000"


def test_route_bins_matches_tiara_rows_by_the_first_header_token(
    tmp_path: Path,
) -> None:
    evidence = _route_one(tmp_path, tiara=(_tiara("contig_1 description", "bacteria"),))

    assert evidence["tiara_prokaryote_fraction"] == "1.0000"


@pytest.mark.parametrize(
    ("first", "second"),
    [
        ("organelle", "unknown"),
        ("organelle", "plastid"),
        ("organelle", "mitochondrion"),
    ],
)
def test_route_bins_withholds_a_route_for_tiara_organelle_bins(
    tmp_path: Path, first: str, second: str
) -> None:
    evidence = _route_one(tmp_path, tiara=(_tiara("contig_1", first, second),))

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "unresolved"
    assert evidence["reason"] == (
        "Tiara classifies most of the bin as organelle sequence "
        "without a gene-calling route"
    )
    assert evidence["tiara_organelle_fraction"] == "1.0000"
    assert json.loads(evidence["supporting_evidence"]) == ["Tiara:1.0000"]


def test_route_bins_marks_tiara_organelle_bins_with_eukaryote_support_conflicting(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "organelle", "plastid"),),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "conflicting"
    assert json.loads(evidence["conflicting_evidence"]) == [
        "organelle evidence conflicts with supported routing evidence"
    ]


def test_route_bins_keeps_the_gvclass_organelle_reason_when_both_agree(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        gvclass="d_PLASTID;p_example",
        tiara=(_tiara("contig_1", "organelle", "plastid"),),
    )

    assert evidence["reason"] == (
        "GVClass reports an organelle class without a gene-calling route"
    )


def test_route_bins_ignores_a_minority_of_organelle_sequence(tmp_path: Path) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 600, "contig_2": 400},
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_2", "organelle", "mitochondrion"),),
    )

    assert evidence["route"] == "eukaryotic"
    assert evidence["tiara_organelle_fraction"] == "0.4000"


@pytest.mark.parametrize(
    ("checkeuk_lineage", "gvclass", "ssu", "tiara", "reason"),
    [
        (
            "",
            "d_EUK;p_example",
            (),
            (_tiara("contig_1", "eukarya"),),
            "GVClass and Tiara support a eukaryotic route",
        ),
        (
            "sg__Sar",
            "d_",
            (),
            (_tiara("contig_1", "eukarya"),),
            "CheckEUK and Tiara support a eukaryotic route",
        ),
        (
            "",
            "d_EUK;p_example",
            (EUK_18S,),
            (),
            (
                "GVClass and a physically linked eukaryotic 18S locus "
                "support a eukaryotic route"
            ),
        ),
        (
            "",
            "d_",
            (EUK_18S,),
            (_tiara("contig_1", "eukarya"),),
            (
                "a physically linked eukaryotic 18S locus and Tiara "
                "support a eukaryotic route"
            ),
        ),
    ],
    ids=["gvclass-tiara", "checkeuk-tiara", "gvclass-18s", "18s-tiara"],
)
def test_route_bins_names_the_eukaryote_lines_that_support_a_route(
    tmp_path: Path,
    checkeuk_lineage: str,
    gvclass: str,
    ssu: tuple[dict[str, str], ...],
    tiara: tuple[dict[str, str], ...],
    reason: str,
) -> None:
    evidence = _route_one(
        tmp_path,
        checkeuk_lineage=checkeuk_lineage,
        gvclass=gvclass,
        ssu=ssu,
        tiara=tiara,
    )

    assert evidence["route"] == "eukaryotic"
    assert evidence["candidate_class"] == "eukaryotic_candidate"
    assert evidence["reason"] == reason


def test_route_bins_lists_the_tiara_eukaryote_fraction_as_support(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "eukarya"),),
    )

    assert json.loads(evidence["supporting_evidence"]) == [
        "GVClass:d_EUK;p_example",
        "Tiara:1.0000",
    ]


def test_route_bins_does_not_count_a_minority_tiara_fraction_as_a_line(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={"contig_1": 499, "contig_2": 501},
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "eukarya"),),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["tiara_eukaryote_fraction"] == "0.4990"


def test_route_bins_conflicts_when_tiara_prokaryote_meets_eukaryote_support(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        tiara=(_tiara("contig_1", "bacteria"),),
    )

    assert evidence["route"] == "unresolved"
    assert evidence["candidate_class"] == "conflicting"
    assert evidence["reason"] == "supported evidence conflicts across routes"
    assert json.loads(evidence["conflicting_evidence"]) == ["Tiara:1.0000"]


def test_route_bins_lets_checkm2_alone_leave_eukaryote_support_intact(
    tmp_path: Path,
) -> None:
    evidence = _route_one(
        tmp_path,
        checkeuk_lineage="sg__Sar",
        gvclass="d_EUK;p_example",
        checkm2=_checkm2("95", "2"),
    )

    assert evidence["route"] == "eukaryotic"


def test_route_bins_records_the_duplicated_contig_fraction(tmp_path: Path) -> None:
    evidence = _route_one(
        tmp_path,
        contigs={
            "u1ctg_len-250_circular-no_depth-2-2-2_duplicated-yes": 250,
            "u2ctg_len-250_circular-no_depth-2-2-2_duplicated-possibly": 250,
            "u3ctg_len-500_circular-no_depth-9-9-9_duplicated-no": 500,
        },
    )

    assert evidence["duplicated_contig_fraction"] == "0.2500"


def test_route_bins_reports_the_tiara_and_checkm2_evidence_in_the_manifest(
    tmp_path: Path,
) -> None:
    _route_one(tmp_path)

    manifest = json.loads((tmp_path / "routing/stage.json").read_text("utf-8"))

    assert manifest["command"][-6:-2] == [
        "--checkm2-dir",
        str(tmp_path / "native/checkm2"),
        "--tiara-dir",
        str(tmp_path / "native/tiara"),
    ]


def test_route_bins_rejects_a_checkm2_row_for_an_unknown_bin(tmp_path: Path) -> None:
    with pytest.raises(
        ValueError, match="CheckM2 reports other, which is not a physical bin"
    ):
        _route_one(
            tmp_path,
            checkm2=(
                *_checkm2("90", "1"),
                {**_checkm2("90", "1")[0], "Name": "other"},
            ),
        )


def test_route_bins_rejects_a_repeated_checkm2_row(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="repeats identifier bin_1"):
        _route_one(tmp_path, checkm2=(*_checkm2("90", "1"), *_checkm2("91", "2")))


def test_route_bins_rejects_a_nonnumeric_checkm2_completeness(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="CheckM2 Completeness"):
        _route_one(tmp_path, checkm2=_checkm2("high", "1"))


@pytest.mark.parametrize("value", ["inf", "-inf", "nan"])
def test_route_bins_rejects_a_nonfinite_checkm2_completeness(
    tmp_path: Path, value: str
) -> None:
    with pytest.raises(ValueError, match="CheckM2 Completeness"):
        _route_one(tmp_path, checkm2=_checkm2(value, "1"))


def test_route_bins_rejects_a_tiara_table_without_a_used_column(
    tmp_path: Path,
) -> None:
    bad = {"sequence_id": "contig_1", "class_fst_stage": "bacteria"}
    with pytest.raises(ValueError, match=r"lacks columns: \['class_snd_stage'\]"):
        _route_one_with_bad_tiara(tmp_path, bad)


def test_route_bins_rejects_an_unsupported_tiara_class(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unsupported Tiara class: viruses"):
        _route_one(tmp_path, tiara=(_tiara("contig_1", "viruses"),))


def test_route_bins_rejects_a_repeated_tiara_sequence(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="repeats identifier contig_1"):
        _route_one(
            tmp_path,
            tiara=(_tiara("contig_1", "bacteria"), _tiara("contig_1", "eukarya")),
        )


def _route_one_with_bad_tiara(tmp_path: Path, row: dict[str, str]) -> None:
    sample_path = tmp_path / "sample.json"
    _sample(sample_path)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    (bin_dir / "bin_1.fa").write_text(">contig_1\nACGT\n", encoding="utf-8")
    directories = _evidence_dirs(tmp_path / "native")
    _write_reports(
        directories,
        quickclade=[
            {"#QueryName": "bin_1.fa", "Q_Bases": "4", "lineage": "", "ConfLevel": ""}
        ],
        checkeuk=[{"genome": "bin_1", "status": "filtered", "lineage": ""}],
        gvclass=[
            {
                "query": "bin_1",
                "taxonomy_majority": "d_",
                "taxonomy_confidence": "low",
                "domain": "",
            }
        ],
        ssu=[],
        viruses=[],
    )
    _write_tsv(
        directories["tiara"] / "tiara.tsv", ("sequence_id", "class_fst_stage"), [row]
    )
    route_bins(sample_path, bin_dir, directories, tmp_path / "routing")


def test_route_command_requires_the_checkm2_and_tiara_directories(
    tmp_path: Path,
) -> None:
    arguments = [
        "route",
        "--sample-json",
        "sample.json",
        "--bins-dir",
        "bins",
        "--quickclade-dir",
        "quickclade",
        "--checkeuk-dir",
        "checkeuk",
        "--gvclass-dir",
        "gvclass",
        "--ssuextract-dir",
        "ssuextract",
        "--viral-dir",
        "viral",
        "--output-dir",
        str(tmp_path / "out"),
    ]

    with pytest.raises(SystemExit) as error:
        main(arguments)

    assert error.value.code == 2
