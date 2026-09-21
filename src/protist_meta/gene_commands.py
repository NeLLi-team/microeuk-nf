"""Prepare normalized protein inputs for functional annotation."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from protist_meta.gene_parsers import ProteinSource, merge_annotation_proteins

__all__ = ["merge_proteins"]

type TsvRow = dict[str, str]

IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
SOURCE_FIELDS = ("source_id", "bin_id", "status")


def merge_proteins(
    prodigal_dir: Path,
    braker_dir: Path,
    output_dir: Path,
) -> Path:
    """Merge completed Prodigal-GV and BRAKER3 proteins with provenance.

    Args:
        prodigal_dir: Prodigal-GV stage directory with ``inputs/sources.tsv``.
        braker_dir: BRAKER3 stage directory with per-bin ``proteins/*.faa``.
        output_dir: New directory for the merged FASTA and identifier map.

    Returns:
        Path to ``proteins.faa`` inside ``output_dir``.

    Raises:
        ValueError: If sources are missing, duplicated, or incomplete.
    """
    if output_dir.exists():
        raise ValueError(f"protein merge output already exists: {output_dir}")
    sources = _prodigal_sources(prodigal_dir)
    sources.extend(_braker_sources(braker_dir))
    if not sources:
        raise ValueError("no completed protein sources to merge")
    output_fasta = output_dir / "proteins.faa"
    merge_annotation_proteins(
        sources,
        output_fasta,
        output_dir / "protein-map.tsv",
    )
    return output_fasta


def _prodigal_sources(stage_dir: Path) -> list[ProteinSource]:
    rows = _read_tsv(stage_dir / "inputs/sources.tsv", SOURCE_FIELDS)
    indexed: dict[str, TsvRow] = {}
    for row in rows:
        source_id = _identifier(row["source_id"], "Prodigal-GV source_id")
        if source_id in indexed:
            raise ValueError(f"duplicate Prodigal-GV source_id: {source_id}")
        if row["status"] not in {"included", "excluded"}:
            raise ValueError(f"invalid Prodigal-GV source status: {row['status']!r}")
        indexed[source_id] = row

    expected = {
        source_id for source_id, row in indexed.items() if row["status"] == "included"
    }
    protein_dir = stage_dir / "proteins"
    observed = {path.stem for path in protein_dir.glob("*.faa") if path.is_file()}
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        message = (
            f"Prodigal-GV protein/source mismatch; missing={missing}; extra={extra}"
        )
        raise ValueError(message)

    sources: list[ProteinSource] = []
    for source_id in sorted(expected):
        bin_id = indexed[source_id]["bin_id"].strip() or None
        if bin_id is not None:
            _identifier(bin_id, "Prodigal-GV bin_id")
        sources.append(
            ProteinSource(
                protein_dir / f"{source_id}.faa",
                "prodigal-gv",
                source_id,
                bin_id,
            )
        )
    return sources


def _braker_sources(stage_dir: Path) -> list[ProteinSource]:
    sources: list[ProteinSource] = []
    for path in sorted((stage_dir / "proteins").glob("*.faa")):
        if not path.is_file() or path.stat().st_size == 0:
            continue
        bin_id = _identifier(path.stem, "BRAKER3 bin_id")
        sources.append(ProteinSource(path, "braker3", "braker3", bin_id))
    return sources


def _read_tsv(path: Path, required: tuple[str, ...]) -> list[TsvRow]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        missing = set(required).difference(fields)
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"{path} contains a malformed TSV row")
    return rows


def _identifier(value: str, field: str) -> str:
    stripped = value.strip()
    if not IDENTIFIER.fullmatch(stripped):
        raise ValueError(f"invalid {field}: {value!r}")
    return stripped
