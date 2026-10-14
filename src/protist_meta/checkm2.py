"""Validate native CheckM2 coverage, including proven empty-protein omissions."""

import argparse
import re
from collections.abc import Sequence, Set as AbstractSet
from pathlib import Path

from protist_meta.quality_parsers import parse_quality

__all__ = ["validate_checkm2_coverage"]

EMPTY_WARNING = re.compile(
    r"^(?:\[[^\]\n]+\] )?WARNING: Skipping protein file "
    r"(?P<path>[^\n]+\.faa) as it was empty\.$",
    re.MULTILINE,
)
NATIVE_ERROR = re.compile(r"\b(?:ERROR|CRITICAL):|Traceback \(most recent call last\):")


def validate_checkm2_coverage(
    checkm2_output_dir: Path, expected_bin_ids: AbstractSet[str], exit_code: int
) -> set[str]:
    """Return omitted IDs only when native success and empty-protein evidence agree."""
    if exit_code != 0:
        raise ValueError(f"CheckM2 prediction exited with {exit_code}")
    reported = set(
        parse_quality("checkm2", checkm2_output_dir / "quality_report.tsv").target_ids
    )
    proteins = list((checkm2_output_dir / "protein_files").glob("*.faa"))
    if {path.stem for path in proteins} != expected_bin_ids or any(
        not path.is_file() or path.is_symlink() for path in proteins
    ):
        raise ValueError("CheckM2 protein files do not cover exactly its input bins")
    log = (checkm2_output_dir / "checkm2.log").read_text(encoding="utf-8")
    if NATIVE_ERROR.search(log):
        raise ValueError("CheckM2 native log contains an error")
    warned = [Path(match["path"]).stem for match in EMPTY_WARNING.finditer(log)]
    if len(warned) != len(set(warned)):
        raise ValueError("CheckM2 native log contains duplicate empty-protein warnings")
    omitted = expected_bin_ids - reported
    empty = {path.stem for path in proteins if path.stat().st_size == 0}
    if not reported <= expected_bin_ids or omitted != empty or omitted != set(warned):
        raise ValueError(
            "CheckM2 output does not cover exactly its input bins with evidence"
        )
    return set(omitted)


def main(argv: Sequence[str] | None = None) -> int:
    """Validate a native output directory against the workflow's input bin files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Native CheckM2 output directory")
    parser.add_argument("bins", type=Path, help="Input directory of .fa bins")
    args = parser.parse_args(argv)
    expected = {path.stem for path in args.bins.glob("*.fa") if path.is_file()}
    exit_path = args.output.parent / "exit_code.txt"
    exit_code = int(exit_path.read_text(encoding="utf-8").strip())
    validate_checkm2_coverage(args.output, expected, exit_code)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
