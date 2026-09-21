"""Execute the shipped checker against native-shaped GVClass receipts."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/check-gvclass-reference.py"


def _write_receipt(path: Path, database: Path, version: str = "v2.0.0") -> Path:
    path.write_text(
        json.dumps({"database": {"path": str(database), "version": version}}),
        encoding="utf-8",
    )
    return path


def _run(receipt: Path, database: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(receipt), str(database), "v2.0.0"],
        check=False,
        capture_output=True,
        text=True,
    )


def test_accepts_matching_receipt_without_modifying_it(tmp_path: Path) -> None:
    database = tmp_path / "resources"
    database.mkdir()
    receipt = _write_receipt(tmp_path / "run_status.json", database)
    before = receipt.read_bytes(), receipt.stat().st_mtime_ns

    completed = _run(receipt, database)

    assert completed.returncode == 0
    assert completed.stderr == ""
    assert (receipt.read_bytes(), receipt.stat().st_mtime_ns) == before


def test_accepts_resolved_path_alias(tmp_path: Path) -> None:
    database = tmp_path / "resources"
    database.mkdir()
    alias = tmp_path / "reference-alias"
    alias.symlink_to(database, target_is_directory=True)
    receipt = _write_receipt(tmp_path / "run_status.json", alias)

    assert _run(receipt, database).returncode == 0


def test_rejects_different_reference_path(tmp_path: Path) -> None:
    receipt = _write_receipt(tmp_path / "run_status.json", tmp_path / "other")

    completed = _run(receipt, tmp_path / "resources")

    assert completed.returncode != 0
    assert "database path differs" in completed.stderr


@pytest.mark.parametrize("version", ["v1.0.0", "v2.1.0", "2.0.0"])
def test_requires_exact_reference_version(tmp_path: Path, version: str) -> None:
    database = tmp_path / "resources"
    receipt = _write_receipt(tmp_path / "run_status.json", database, version)

    completed = _run(receipt, database)

    assert completed.returncode != 0
    assert "database version differs" in completed.stderr


@pytest.mark.parametrize(
    "contents",
    [
        "{}",
        '{"database": {}}',
        '{"database": {"path": "/reference"}}',
        '{"database": {"version": "v2.0.0"}}',
        '{"database": {"path": null, "version": "v2.0.0"}}',
        '{"database": []}',
        "[]",
        "{broken",
    ],
)
def test_rejects_missing_or_malformed_metadata(tmp_path: Path, contents: str) -> None:
    receipt = tmp_path / "run_status.json"
    receipt.write_text(contents, encoding="utf-8")

    completed = _run(receipt, tmp_path / "resources")

    assert completed.returncode != 0
    assert "GVClass reference check failed:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_rejects_missing_receipt(tmp_path: Path) -> None:
    completed = _run(tmp_path / "missing.json", tmp_path / "resources")

    assert completed.returncode != 0
    assert "GVClass reference check failed:" in completed.stderr
    assert "Traceback" not in completed.stderr
