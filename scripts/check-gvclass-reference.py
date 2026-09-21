"""Require a native GVClass receipt to match the configured reference."""

import argparse
import json
from collections.abc import Sequence
from pathlib import Path


def main(argv: Sequence[str] | None = None) -> int:
    """Check resolved database paths and exact versions without changing files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run_status_json", type=Path, help="Native GVClass run_status.json receipt."
    )
    parser.add_argument(
        "expected_database_path",
        type=Path,
        help="Configured database directory; symlinks are resolved.",
    )
    parser.add_argument(
        "expected_version", help="Exact configured database version label."
    )
    args = parser.parse_args(argv)
    prefix = "GVClass reference check failed: "
    try:
        receipt = json.loads(args.run_status_json.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        parser.exit(1, f"{prefix}{error}\n")
    if not isinstance(receipt, dict) or not isinstance(receipt.get("database"), dict):
        parser.exit(1, f"{prefix}native database metadata is missing\n")
    database = receipt["database"]
    actual_path = database.get("path")
    actual_version = database.get("version")
    if not isinstance(actual_path, str) or actual_path == "":
        parser.exit(1, f"{prefix}native database path is missing or invalid\n")
    if not isinstance(actual_version, str) or actual_version == "":
        parser.exit(1, f"{prefix}native database version is missing or invalid\n")
    if Path(actual_path).resolve() != args.expected_database_path.resolve():
        parser.exit(
            1, f"{prefix}native database path differs from configured reference\n"
        )
    if actual_version != args.expected_version:
        parser.exit(
            1, f"{prefix}native database version differs from configured reference\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
