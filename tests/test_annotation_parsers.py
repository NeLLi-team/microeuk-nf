import json
import sqlite3
from pathlib import Path

import pytest

from protist_meta.annotation_parsers import parse_eggnog, parse_interpro
from protist_meta.catalog import build_catalog

EGGNOG_HEADER = (
    "#query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
)


def _write_eggnog(path: Path, rows: list[dict[str, str]]):
    lines = [
        "## emapper-2.1.15 / Expected eggNOG DB version: 5.0.2",
        "\t".join(EGGNOG_HEADER),
    ]
    lines.extend(
        "\t".join(row.get(field, "-") for field in EGGNOG_HEADER) for row in rows
    )
    lines.append("## 1 queries scanned")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_parse_eggnog_keeps_exact_supported_terms_without_coordinates(tmp_path):
    path = tmp_path / "annotations.emapper.annotations"
    gene_id = "braker:bin_1:transcript_1.p1"
    _write_eggnog(
        path,
        [
            {
                "#query": gene_id,
                "seed_ortholog": "9606.ENSP000001",
                "evalue": "1.2e-40",
                "score": "212.5",
                "eggNOG_OGs": "KOG0002@2759|Eukaryota,COG0001@1|root",
                "max_annot_lvl": "2759|Eukaryota",
                "COG_category": "O",
                "Description": "ATP-dependent chaperone",
                "Preferred_name": "hsp70",
                "GOs": "GO:0005524,GO:0000166",
                "EC": "3.6.4.12",
                "KEGG_ko": "ko:K04043",
            }
        ],
    )

    rows = parse_eggnog(path, {gene_id, "prodigal:bin_2:no_hit"}, "eggNOG 5.0.2")

    assert rows == parse_eggnog(
        path, {gene_id, "prodigal:bin_2:no_hit"}, "eggNOG 5.0.2"
    )
    assert len(rows) == 6
    assert len({row["annotation_id"] for row in rows}) == len(rows)
    assert {row["database_name"] for row in rows} == {
        "eggNOG",
        "Gene Ontology via eggNOG",
        "Enzyme Commission via eggNOG",
        "KEGG Orthology via eggNOG",
    }
    assert {row["accession"] for row in rows} == {
        "KOG0002@2759",
        "COG0001@1",
        "GO:0005524",
        "GO:0000166",
        "3.6.4.12",
        "ko:K04043",
    }
    assert all(row["gene_id"] == gene_id for row in rows)
    assert all(row["database_version"] == "eggNOG 5.0.2" for row in rows)
    assert all("start" not in row and "end" not in row for row in rows)
    assert all(row["score"] == 212.5 and row["evalue"] == 1.2e-40 for row in rows)
    orthologous_groups = [row for row in rows if row["database_name"] == "eggNOG"]
    assert {row["function_name"] for row in orthologous_groups} == {
        "ATP-dependent chaperone"
    }


def test_parse_eggnog_accepts_a_header_without_hits(tmp_path):
    path = tmp_path / "empty.emapper.annotations"
    _write_eggnog(path, [])

    assert parse_eggnog(path, {"prodigal:bin_1:protein_1"}, "eggNOG 5.0.2") == []


def test_parse_eggnog_normalizes_taxon_labels_without_changing_description(
    tmp_path: Path,
) -> None:
    path = tmp_path / "annotations.emapper.annotations"
    _write_eggnog(
        path,
        [
            {
                "#query": "protein_1",
                "eggNOG_OGs": ("4QV4U@35237|dsDNA viruses  no RNA stage,4QV4U@35237"),
                "Description": "Native protein function",
            }
        ],
    )

    rows = parse_eggnog(path, {"protein_1"}, "5.0.2")

    assert len(rows) == 1
    assert rows[0]["accession"] == "4QV4U@35237"
    assert rows[0]["function_name"] == "Native protein function"


@pytest.mark.parametrize(
    "term",
    [
        "COG0001",
        "@1|root",
        "COG0001@root|root",
        "COG0001@1|",
        "COG0001@1|   ",
        "COG0001@1|root|extra",
        "COG0001 @1|root",
        "COG0001@1|root,",
    ],
)
def test_parse_eggnog_rejects_malformed_orthologous_groups(
    tmp_path: Path,
    term: str,
) -> None:
    path = tmp_path / "annotations.emapper.annotations"
    _write_eggnog(path, [{"#query": "protein_1", "eggNOG_OGs": term}])

    with pytest.raises(ValueError, match="invalid eggNOG orthologous group"):
        parse_eggnog(path, {"protein_1"}, "5.0.2")


@pytest.mark.parametrize(
    ("field", "term"),
    [
        ("GOs", "GO:0005524 extra"),
        ("EC", "3.6.4.12 extra"),
        ("KEGG_ko", "ko:K04043 extra"),
    ],
)
def test_parse_eggnog_keeps_other_accessions_whitespace_free(
    tmp_path: Path,
    field: str,
    term: str,
) -> None:
    path = tmp_path / "annotations.emapper.annotations"
    _write_eggnog(path, [{"#query": "protein_1", field: term}])

    with pytest.raises(ValueError, match="malformed functional accession list"):
        parse_eggnog(path, {"protein_1"}, "5.0.2")


def test_parse_interpro_keeps_member_and_integrated_coordinates(tmp_path):
    path = tmp_path / "proteins.tsv"
    gene_id = "braker:bin_1:transcript_1.p1"
    path.write_text(
        f"{gene_id}\t0123456789abcdef0123456789abcdef\t100\tPfam\tPF00001\t"
        "Example ATPase domain\t3\t45\t2.4E-7\tT\t19-09-2026\tIPR000001\t"
        "Example integrated domain\n"
        f"{gene_id}\t0123456789abcdef0123456789abcdef\t100\tMobiDBLite\t"
        "mobidb-lite\t-\t60\t90\t-\tT\t19-09-2026\n",
        encoding="utf-8",
    )

    rows = parse_interpro(path, {gene_id, "prodigal:bin_2:no_hit"}, "5.76-107.0")

    assert len(rows) == 3
    assert len({row["annotation_id"] for row in rows}) == 3
    assert {row["database_name"] for row in rows} == {
        "InterProScan:Pfam",
        "InterProScan:MobiDBLite",
        "InterPro",
    }
    assert {row["accession"] for row in rows} == {
        "PF00001",
        "mobidb-lite",
        "IPR000001",
    }
    pfam = next(row for row in rows if row["accession"] == "PF00001")
    assert pfam["start"] == 3
    assert pfam["end"] == 45
    assert pfam["score"] == 2.4e-7
    assert "evalue" not in pfam
    mobidb = next(row for row in rows if row["accession"] == "mobidb-lite")
    assert mobidb["start"] == 60
    assert mobidb["end"] == 90
    assert "score" not in mobidb
    assert "function_name" not in mobidb


def test_parse_interpro_collapses_integrated_keys_without_member_scores(
    tmp_path: Path,
) -> None:
    path = tmp_path / "proteins.tsv"
    path.write_text(
        "protein_1\t0123456789abcdef0123456789abcdef\t100\tPfam\tPF00001\t"
        "First member\t12\t66\t6.4e-13\tT\t20-09-2026\tIPR001387\t"
        "Cro/C1-type, helix-turn-helix domain\n"
        "protein_1\t0123456789abcdef0123456789abcdef\t100\tProSiteProfiles\t"
        "PS00001\tSecond member\t12\t66\t14.363916\tT\t20-09-2026\tIPR001387\t"
        "Cro/C1-type, helix-turn-helix domain\n",
        encoding="utf-8",
    )

    rows = parse_interpro(path, {"protein_1"}, "5.76-107.0")

    assert len(rows) == 3
    assert len({row["annotation_id"] for row in rows}) == 3
    assert {row["database_name"]: row.get("score") for row in rows} == {
        "InterProScan:Pfam": 6.4e-13,
        "InterProScan:ProSiteProfiles": 14.363916,
        "InterPro": None,
    }
    integrated = next(row for row in rows if row["database_name"] == "InterPro")
    assert "score" not in integrated
    assert (integrated["start"], integrated["end"]) == (12, 66)
    assert integrated["function_name"] == "Cro/C1-type, helix-turn-helix domain"


def test_parse_interpro_accepts_no_hit_output(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text("# InterProScan TSV has no native header\n", encoding="utf-8")

    assert parse_interpro(path, {"prodigal:bin_1:protein_1"}, "5.76-107.0") == []


def test_corrected_annotations_publish_in_existing_catalog_fixture(
    tmp_path: Path,
) -> None:
    """Publish parser rows in the existing fixture, without claiming native science."""
    eggnog = tmp_path / "annotations.emapper.annotations"
    _write_eggnog(
        eggnog,
        [
            {
                "#query": "gene1",
                "eggNOG_OGs": "4QV4U@35237|dsDNA viruses  no RNA stage",
                "Description": "Native protein function",
            }
        ],
    )
    interpro = tmp_path / "proteins.tsv"
    interpro.write_text(
        "gene1\t0123456789abcdef0123456789abcdef\t100\tPfam\tPF00001\t"
        "First member\t12\t66\t6.4e-13\tT\t20-09-2026\tIPR001387\tDomain\n"
        "gene1\t0123456789abcdef0123456789abcdef\t100\tProSiteProfiles\tPS00001\t"
        "Second member\t12\t66\t14.363916\tT\t20-09-2026\tIPR001387\tDomain\n",
        encoding="utf-8",
    )
    rows = parse_eggnog(eggnog, {"gene1"}, "5.0.2") + parse_interpro(
        interpro, {"gene1"}, "5.76-107.0"
    )
    fixture = Path(__file__).parent / "fixtures/catalog.json"
    bundle = json.loads(fixture.read_text(encoding="utf-8"))
    scope = {
        key: bundle["annotations"][0][key]
        for key in (
            "sample_id",
            "run_id",
            "stage_id",
            "artifact_id",
        )
    }
    bundle["annotations"] = [{**scope, **row} for row in rows]
    annotation_stage = next(
        stage for stage in bundle["stages"] if stage["stage_id"] == "annotations"
    )
    annotation_stage["expected_result_keys"] = [
        f"annotations:{row['annotation_id']}" for row in rows
    ]
    source = tmp_path / "bundle.json"
    source.write_text(json.dumps(bundle), encoding="utf-8")
    database = build_catalog(source, tmp_path / "catalog.sqlite")

    with sqlite3.connect(database) as connection:
        assert connection.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute(
            "SELECT database_name, accession, score "
            "FROM annotations ORDER BY database_name"
        ).fetchall() == [
            ("InterPro", "IPR001387", None),
            ("InterProScan:Pfam", "PF00001", 6.4e-13),
            ("InterProScan:ProSiteProfiles", "PS00001", 14.363916),
            ("eggNOG", "4QV4U@35237", None),
        ]


def test_annotation_parsers_reject_unknown_protein_ids(tmp_path):
    eggnog = tmp_path / "annotations.emapper.annotations"
    _write_eggnog(
        eggnog,
        [
            {
                "#query": "unknown",
                "evalue": "1e-5",
                "score": "20",
                "eggNOG_OGs": "COG0001@1|root",
            }
        ],
    )
    interpro = tmp_path / "proteins.tsv"
    interpro.write_text(
        "unknown\t0123456789abcdef0123456789abcdef\t50\tPfam\tPF00001\tDomain\t1\t20\t1e-5\tT\t19-09-2026\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="unknown protein unknown"):
        parse_eggnog(eggnog, {"known"}, "eggNOG 5.0.2")
    with pytest.raises(ValueError, match="unknown protein unknown"):
        parse_interpro(interpro, {"known"}, "5.76-107.0")
