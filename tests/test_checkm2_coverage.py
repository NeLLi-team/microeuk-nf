"""Require native evidence for CheckM2 bins omitted after gene calling."""

from pathlib import Path

import pytest

from protist_meta.checkm2 import validate_checkm2_coverage


@pytest.fixture
def output(tmp_path: Path) -> Path:
    """Provide a scored bin and a native empty-protein omission."""
    output = tmp_path / "output"
    proteins = output / "protein_files"
    proteins.mkdir(parents=True)
    (proteins / "scored.faa").write_text(">scored_1\nMPEPTIDE\n", encoding="utf-8")
    (proteins / "empty.faa").touch()
    (output / "quality_report.tsv").write_text(
        "Name\tCompleteness\tContamination\tCompleteness_Model_Used\n"
        "scored\t90\t103.62\tNeural Network (Specific Model)\n",
        encoding="utf-8",
    )
    (output / "checkm2.log").write_text(
        "[10/14/2026 03:21:22 AM] WARNING: Skipping protein file "
        "output/protein_files/empty.faa as it was empty.\n",
        encoding="utf-8",
    )
    return output


def test_validates_native_empty_protein_omission(output: Path) -> None:
    assert validate_checkm2_coverage(output, {"scored", "empty"}, 0) == {"empty"}


@pytest.mark.parametrize(
    ("relative", "contents", "exit_code"),
    [
        ("protein_files/empty.faa", ">empty_1\nMPEPTIDE\n", 0),
        ("protein_files/empty.faa", None, 0),
        ("checkm2.log", "", 0),
        ("checkm2.log", "WARNING: empty protein file empty.faa\n", 0),
        (
            "checkm2.log",
            "WARNING: Skipping protein file foreign.faa as it was empty.\n",
            0,
        ),
        (
            "checkm2.log",
            "WARNING: Skipping protein file empty.faa as it was empty.\n" * 2,
            0,
        ),
        ("checkm2.log", "ERROR: native failure\n", 0),
        ("protein_files/foreign.faa", "", 0),
        (
            "quality_report.tsv",
            (
                "Name\tCompleteness\tContamination\tCompleteness_Model_Used\n"
                "foreign\t90\t1\tNeural Network (Specific Model)\n"
            ),
            0,
        ),
        ("unused", "", 1),
    ],
)
def test_rejects_unproven_omissions(
    output: Path, relative: str, contents: str | None, exit_code: int
) -> None:
    path = output / relative
    if contents is None:
        path.unlink()
    else:
        path.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match="CheckM2"):
        validate_checkm2_coverage(output, {"scored", "empty"}, exit_code)
