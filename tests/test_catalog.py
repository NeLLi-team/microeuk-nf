import copy
import json
import sqlite3
from pathlib import Path

import pytest

from protist_meta.catalog import (
    CatalogValidationError,
    build_catalog,
    validate_bundle,
)

FIXTURE = Path(__file__).parent / "fixtures" / "catalog.json"


def load_fixture():
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def write_bundle(tmp_path, bundle):
    bundle_path = tmp_path / "catalog.json"
    bundle_path.write_text(json.dumps(bundle), encoding="utf-8")
    return bundle_path


def test_build_catalog_normalizes_and_checks_valid_bundle(tmp_path):
    output = tmp_path / "catalog.sqlite"

    published = build_catalog(FIXTURE, output)

    assert published == output.resolve()
    with sqlite3.connect(output) as connection:
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_schema WHERE type = 'table'"
            )
        }
        assert len(tables) == 19
        assert connection.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute(
            "SELECT schema_version, workflow_name FROM catalog_metadata"
        ).fetchone() == ("1.0.0", "protist-meta-nf")
        assert connection.execute(
            "SELECT checkeuk_pfam_ridge_completeness_percent, "
            "checkeuk_bom_completeness_percent, checkeuk_inventory74_present, "
            "checkeuk_expected_markers_present, checkeuk_expected_markers_total, "
            "checkeuk_pfam1085_distinct FROM qc"
        ).fetchone() == (81.5, 78.4, 61, 58, 70, 742)
        assert connection.execute(
            "SELECT count(*) FROM stage_expected_keys WHERE stage_id = 'empty-checkv'"
        ).fetchone() == (0,)


def test_catalog_preserves_model_when_caller_is_unknown(tmp_path):
    bundle = load_fixture()
    run = bundle["runs"][0]
    run["basecaller"] = None
    run["basecaller_model"] = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
    output = tmp_path / "catalog.sqlite"

    build_catalog(write_bundle(tmp_path, bundle), output)

    with sqlite3.connect(output) as connection:
        assert connection.execute(
            "SELECT basecaller, basecaller_model FROM runs"
        ).fetchone() == (None, run["basecaller_model"])


def test_validate_bundle_rejects_unknown_field_and_strict_type(tmp_path):
    unknown = load_fixture()
    unknown["samples"][0]["unexpected"] = "value"
    wrong_type = load_fixture()
    wrong_type["samples"][0]["genetic_code"] = "1"

    with pytest.raises(CatalogValidationError, match="extra_forbidden"):
        validate_bundle(write_bundle(tmp_path, unknown))
    with pytest.raises(CatalogValidationError, match="int_type"):
        validate_bundle(write_bundle(tmp_path, wrong_type))


def test_catalog_accepts_one_native_report_with_multiple_roles(tmp_path):
    bundle = load_fixture()
    source = bundle["artifacts"][0]
    second_role = {**source, "artifact_id": "reads-provenance", "kind": "other"}
    bundle["artifacts"].append(second_role)
    output = tmp_path / "catalog.sqlite"

    build_catalog(write_bundle(tmp_path, bundle), output)

    with sqlite3.connect(output) as connection:
        assert connection.execute(
            "SELECT count(*) FROM artifacts WHERE path = ?", (source["path"],)
        ).fetchone() == (2,)


def test_catalog_rejects_conflicting_content_for_one_path(tmp_path):
    bundle = load_fixture()
    source = bundle["artifacts"][0]
    bundle["artifacts"].append(
        {**source, "artifact_id": "changed-reads", "sha256": "f" * 64}
    )

    with pytest.raises(CatalogValidationError, match="conflicting artifact content"):
        validate_bundle(write_bundle(tmp_path, bundle))


def test_build_catalog_keeps_absent_qc_measurement_null(tmp_path):
    bundle = load_fixture()
    bundle["qc"][0].pop("checkeuk_bom_completeness_percent")
    output = tmp_path / "catalog.sqlite"

    build_catalog(write_bundle(tmp_path, bundle), output)

    with sqlite3.connect(output) as connection:
        assert connection.execute(
            "SELECT q.checkeuk_bom_completeness_percent, "
            "q.completeness_query_support, q.completeness_model_reliability, "
            "q.contamination_query_support, q.contamination_model_reliability, "
            "t.lineage_confidence FROM qc q CROSS JOIN taxonomy t"
        ).fetchone() == (None, None, None, None, None, None)


def test_build_catalog_preserves_gvclass_confidence_labels(tmp_path):
    bundle = load_fixture()
    bundle["qc"][0].update(
        {
            "completeness_query_support": "unassigned",
            "completeness_model_reliability": "unavailable",
            "contamination_query_support": "supported",
            "contamination_model_reliability": "unvalidated_legacy",
        }
    )
    bundle["taxonomy"][0]["lineage_confidence"] = "low_support,reduced_fastmode"
    output = tmp_path / "catalog.sqlite"

    build_catalog(write_bundle(tmp_path, bundle), output)

    with sqlite3.connect(output) as connection:
        assert connection.execute(
            "SELECT q.completeness_query_support, "
            "q.completeness_model_reliability, q.contamination_query_support, "
            "q.contamination_model_reliability, t.lineage_confidence "
            "FROM qc q CROSS JOIN taxonomy t"
        ).fetchone() == (
            "unassigned",
            "unavailable",
            "supported",
            "unvalidated_legacy",
            "low_support,reduced_fastmode",
        )


def test_validate_bundle_accepts_annotation_without_coordinates(
    tmp_path: Path,
) -> None:
    bundle = load_fixture()
    annotation = bundle["annotations"][0]
    annotation.pop("start")
    annotation.pop("end")

    validated = validate_bundle(write_bundle(tmp_path, bundle))

    assert validated.annotations[0].start is None
    assert validated.annotations[0].end is None


def test_validate_bundle_rejects_half_populated_annotation_coordinates(
    tmp_path: Path,
) -> None:
    bundle = load_fixture()
    bundle["annotations"][0].pop("end")

    with pytest.raises(
        CatalogValidationError,
        match="annotation coordinates must appear together",
    ):
        validate_bundle(write_bundle(tmp_path, bundle))


def test_genetic_code_override_only_constrains_prokaryotic_bins(
    tmp_path: Path,
) -> None:
    """A mixed sample retains the separate eukaryotic translation table."""
    bundle = load_fixture()
    bundle["samples"][0]["genetic_code"] = 11
    bundle["genes"][0]["genetic_code"] = 1
    bundle["bins"][0]["candidate_class"] = "eukaryotic_candidate"

    validate_bundle(write_bundle(tmp_path, bundle))

    bundle["bins"][0]["candidate_class"] = "bacterial_archaeal_candidate"
    with pytest.raises(CatalogValidationError, match="prokaryotic override"):
        validate_bundle(write_bundle(tmp_path, bundle))


def test_validate_bundle_rejects_duplicate_and_foreign_key(tmp_path):
    duplicate = load_fixture()
    duplicate["bins"].append(copy.deepcopy(duplicate["bins"][0]))
    broken_reference = load_fixture()
    broken_reference["memberships"][0]["bin_id"] = "missing-bin"

    with pytest.raises(CatalogValidationError, match="duplicate bins key"):
        validate_bundle(write_bundle(tmp_path, duplicate))
    with pytest.raises(CatalogValidationError, match="missing membership bin"):
        validate_bundle(write_bundle(tmp_path, broken_reference))


def test_validate_bundle_requires_exact_completed_stage_coverage(tmp_path):
    bundle = load_fixture()
    assemble = next(
        stage for stage in bundle["stages"] if stage["stage_id"] == "assemble"
    )
    assemble["expected_result_keys"].remove("contigs:ctg1")

    with pytest.raises(CatalogValidationError, match="key coverage differs"):
        validate_bundle(write_bundle(tmp_path, bundle))


def test_invalid_rebuild_preserves_existing_catalog(tmp_path):
    output = tmp_path / "catalog.sqlite"
    build_catalog(FIXTURE, output)
    original = output.read_bytes()
    bundle = load_fixture()
    bundle["contigs"][0]["assembly_id"] = "missing-assembly"

    with pytest.raises(CatalogValidationError, match="missing contig assembly"):
        build_catalog(write_bundle(tmp_path, bundle), output)

    assert output.read_bytes() == original
