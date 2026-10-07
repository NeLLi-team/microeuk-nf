import csv
import json
from pathlib import Path

import pytest

from protist_meta.collect import StageManifest
from protist_meta.routing import route_bins


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
        )
    }


def _write_reports(
    directories: dict[str, Path],
    quickclade: list[dict[str, str]],
    checkeuk: list[dict[str, str]],
    gvclass: list[dict[str, str]],
    ssu: list[dict[str, str]],
    viruses: list[dict[str, str]],
):
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
        directories["gvclass"] / "output/gvclass_summary.tsv",
        ("query", "taxonomy_majority", "taxonomy_confidence", "domain"),
        gvclass,
    )
    _write_tsv(
        directories["ssuextract"] / "output/cmsearch_summary.tsv",
        ("name", "model", "contig_name", "taxonomy", "taxonomy_domain"),
        ssu,
    )
    _write_tsv(
        directories["viral"] / "genomad/input_summary/input_virus_summary.tsv",
        ("seq_name", "virus_score", "taxonomy"),
        viruses,
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
