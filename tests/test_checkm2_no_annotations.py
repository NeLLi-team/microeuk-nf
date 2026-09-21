"""Exercise the shipped CheckM2 no-annotation checker on native-shaped outputs."""

import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/checkm2-no-annotations.sh"
NO_HITS = (
    "[09/19/2026 08:10:22 PM] ERROR: No DIAMOND annotation was generated. Exiting\n"
)


@pytest.fixture
def native_case(tmp_path: Path) -> tuple[Path, Path]:
    output = tmp_path / "output"
    bins = tmp_path / "bins"
    bins.mkdir()
    (output / "protein_files").mkdir(parents=True)
    (output / "diamond_output").mkdir()
    (output / "checkm2.log").write_text(
        "[09/19/2026 08:10:22 PM] INFO: Processing DIAMOND output\n" + NO_HITS,
        encoding="utf-8",
    )
    (output / "diamond_output/DIAMOND_RESULTS.tsv").touch()
    (bins / "alpha.fa").write_text(">alpha\nATGAAA\n", encoding="utf-8")
    (bins / "beta.fa").write_text(">beta\nATGCCC\n", encoding="utf-8")
    (output / "protein_files/alpha.faa").write_text(">alpha_1\nMK\n", encoding="utf-8")
    (output / "protein_files/beta.faa").write_text(">beta_1\nMP\n", encoding="utf-8")
    return output, bins


def _run(native_case: tuple[Path, Path], exit_code: str = "1") -> int:
    output, bins = native_case
    return subprocess.run(
        ["bash", str(SCRIPT), exit_code, str(output), str(bins)],
        check=False,
        capture_output=True,
        text=True,
    ).returncode


def _snapshot(root: Path) -> dict[Path, tuple[bytes, int]]:
    return {
        path: (path.read_bytes(), path.stat().st_mtime_ns)
        for path in root.rglob("*")
        if path.is_file()
    }


def test_accepts_known_no_hit_without_modifying_files(
    native_case: tuple[Path, Path], tmp_path: Path
) -> None:
    output, _bins = native_case
    (output / "diamond_output/DIAMOND_RESULTS_1.tsv").touch()
    before = _snapshot(tmp_path)

    assert _run(native_case) == 0
    assert _snapshot(tmp_path) == before


@pytest.mark.parametrize("exit_code", ["0", "2", "137", "01", ""])
def test_rejects_unexpected_exit(
    native_case: tuple[Path, Path], exit_code: str
) -> None:
    assert _run(native_case, exit_code) != 0


@pytest.mark.parametrize(
    "log",
    [
        "[09/19/2026 08:10:22 PM] INFO: No DIAMOND annotation was generated. Exiting\n",
        "[09/19/2026 08:10:22 PM] ERROR: DIAMOND database was not found\n",
        NO_HITS + NO_HITS,
        NO_HITS + "[09/19/2026 08:10:23 PM] ERROR: Failed to read proteins\n",
        NO_HITS.rstrip() + " unexpectedly\n",
    ],
)
def test_rejects_unexpected_log(native_case: tuple[Path, Path], log: str) -> None:
    output, _bins = native_case
    (output / "checkm2.log").write_text(log, encoding="utf-8")

    assert _run(native_case) != 0


@pytest.mark.parametrize(
    "relative_path",
    ["checkm2.log", "diamond_output/DIAMOND_RESULTS.tsv", "protein_files/alpha.faa"],
)
def test_rejects_missing_evidence(
    native_case: tuple[Path, Path], relative_path: str
) -> None:
    output, _bins = native_case
    (output / relative_path).unlink()

    assert _run(native_case) != 0


def test_rejects_nonempty_diamond_shard(native_case: tuple[Path, Path]) -> None:
    output, _bins = native_case
    (output / "diamond_output/DIAMOND_RESULTS_1.tsv").write_text(
        "protein\treference\n", encoding="utf-8"
    )

    assert _run(native_case) != 0


@pytest.mark.parametrize("report", ["", "Name\tCompleteness\nalpha\t0\n"])
def test_rejects_partial_report(native_case: tuple[Path, Path], report: str) -> None:
    output, _bins = native_case
    (output / "quality_report.tsv").write_text(report, encoding="utf-8")

    assert _run(native_case) != 0


def test_rejects_extra_protein(native_case: tuple[Path, Path]) -> None:
    output, _bins = native_case
    (output / "protein_files/extra.faa").write_text(">extra_1\nMK\n", encoding="utf-8")

    assert _run(native_case) != 0


def test_rejects_mismatched_protein_names(native_case: tuple[Path, Path]) -> None:
    output, _bins = native_case
    (output / "protein_files/alpha.faa").rename(output / "protein_files/other.faa")

    assert _run(native_case) != 0


def test_rejects_empty_protein(native_case: tuple[Path, Path]) -> None:
    output, _bins = native_case
    (output / "protein_files/alpha.faa").write_text("", encoding="utf-8")

    assert _run(native_case) != 0


def test_rejects_nonregular_diamond_result(native_case: tuple[Path, Path]) -> None:
    output, _bins = native_case
    (output / "diamond_output/DIAMOND_RESULTS_1.tsv").mkdir()

    assert _run(native_case) != 0
