"""Prepare exclusive gene-calling FASTAs and coordinate origin ledgers."""

from __future__ import annotations

import csv
import hashlib
import re
import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from protist_meta.sequences import read_fasta

__all__ = ["lift_prodigal_gene", "load_prodigal_origins", "prepare_gene_inputs"]

type GeneRow = dict[str, object]
type TsvRow = dict[str, str]

IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
PROVIRUS_ID = re.compile(
    r"^(?P<parent>.+)\|provirus_(?P<start>[1-9][0-9]*)_(?P<end>[1-9][0-9]*)$"
)
CALLED_ROUTES = ("prokaryotic", "viral", "eukaryotic")
ALL_ROUTES = (*CALLED_ROUTES, "unresolved")
SOURCE_FIELDS = (
    "source_id",
    "category",
    "bin_id",
    "fasta",
    "origin_map",
    "status",
    "reason",
)
ORIGIN_FIELDS = (
    "contig_id",
    "parent_contig",
    "start",
    "end",
    "source_id",
    "bin_id",
    "category",
    "status",
    "reason",
)
MISSING = {"", "na", "n/a", "none", "null", "."}
GENOMAD_SOURCE = "genomad-unbinned"


@dataclass(frozen=True, slots=True)
class _RoutedSource:
    route: str
    bin_id: str
    path: Path
    contigs: tuple[tuple[str, int], ...]


@dataclass(frozen=True, slots=True)
class _Candidate:
    native_id: str
    contig_id: str
    parent_contig: str
    start: int
    end: int
    sequence: str
    bin_id: str
    category: str
    status: str
    reason: str


def prepare_gene_inputs(
    routing_dir: Path,
    viral_dir: Path,
    output_dir: Path,
) -> Path:
    """Stage each selected bin once and add nonduplicating geNomad candidates.

    Routed prokaryotic and viral bins are staged for Prodigal-GV. A geNomad
    whole-contig or proviral sequence is added only when its parent is not in a
    prokaryotic, viral, or eukaryotic bin. Exclusions remain in ``origins.tsv``.

    Args:
        routing_dir: Routing stage output with four route subdirectories.
        viral_dir: Viral stage output containing the geNomad input summary.
        output_dir: New directory for staged FASTAs and provenance tables.

    Returns:
        Path to ``sources.tsv`` inside ``output_dir``.

    Raises:
        ValueError: If routing or geNomad records are empty, duplicated, or
            inconsistent.
    """
    if output_dir.exists():
        raise ValueError(f"gene input output already exists: {output_dir}")
    routed, contig_routes = _routed_sources(routing_dir)
    candidates = _genomad_candidates(viral_dir, contig_routes)
    _validate_candidate_intervals(candidates)

    sources_dir = output_dir / "sources"
    origin_maps_dir = output_dir / "origins"
    sources_dir.mkdir(parents=True)
    origin_maps_dir.mkdir()
    source_rows: list[TsvRow] = []
    origin_rows: list[TsvRow] = []

    for source in routed:
        if source.route not in {"prokaryotic", "viral"}:
            continue
        source_row, rows = _stage_routed_source(
            source, output_dir, sources_dir, origin_maps_dir
        )
        source_rows.append(source_row)
        origin_rows.extend(rows)

    included_candidates = [
        candidate for candidate in candidates if candidate.status == "included"
    ]
    if included_candidates:
        source_row, rows = _stage_genomad_source(
            included_candidates, output_dir, sources_dir, origin_maps_dir
        )
        source_rows.append(source_row)
        origin_rows.extend(rows)
    origin_rows.extend(
        _origin_row(candidate)
        for candidate in candidates
        if candidate.status == "excluded"
    )

    _write_tsv(output_dir / "sources.tsv", SOURCE_FIELDS, source_rows)
    _write_tsv(output_dir / "origins.tsv", ORIGIN_FIELDS, origin_rows)
    return output_dir / "sources.tsv"


def load_prodigal_origins(path: Path) -> dict[str, TsvRow]:
    """Read one source's origin TSV and reject duplicate contig identifiers."""
    return _index_rows(
        _read_tsv(path, ORIGIN_FIELDS), "contig_id", "Prodigal-GV origins"
    )


def lift_prodigal_gene(
    row: Mapping[str, object], origins: Mapping[str, TsvRow]
) -> GeneRow:
    """Lift one Prodigal-GV row from a staged sequence to its parent contig.

    Args:
        row: Parsed gene row using the staged FASTA contig identifier.
        origins: Per-source origins indexed once by ``load_prodigal_origins``.

    Returns:
        A copy with parent-contig coordinates and the originating bin ID.

    Raises:
        ValueError: If the identifier or coordinates lack one exact origin.
    """
    contig_id = row.get("contig_id")
    start = row.get("start")
    end = row.get("end")
    if not isinstance(contig_id, str) or not contig_id:
        raise ValueError("gene row requires a nonempty contig_id")
    if (
        not isinstance(start, int)
        or isinstance(start, bool)
        or not isinstance(end, int)
        or isinstance(end, bool)
        or start < 1
        or end < start
    ):
        raise ValueError("gene row has invalid one-based coordinates")

    origin = origins.get(contig_id)
    if origin is None:
        raise ValueError(f"gene contig requires one origin mapping: {contig_id}")
    if origin["status"] != "included":
        raise ValueError(f"gene contig has an excluded origin: {contig_id}")
    origin_start = _positive_integer(origin["start"], "origin start")
    origin_end = _positive_integer(origin["end"], "origin end")
    if end > origin_end - origin_start + 1:
        raise ValueError(f"gene coordinates exceed staged sequence: {contig_id}")

    lifted = dict(row)
    lifted["contig_id"] = origin["parent_contig"]
    lifted["start"] = origin_start + start - 1
    lifted["end"] = origin_start + end - 1
    lifted["bin_id"] = origin["bin_id"] or None
    return lifted


def _routed_sources(
    routing_dir: Path,
) -> tuple[list[_RoutedSource], dict[str, tuple[str, str]]]:
    sources: list[_RoutedSource] = []
    contig_routes: dict[str, tuple[str, str]] = {}
    for route in ALL_ROUTES:
        route_dir = routing_dir / route
        if not route_dir.is_dir():
            raise ValueError(f"routing directory is missing: {route_dir}")
        for path in sorted(route_dir.glob("*.fna")):
            records = tuple(
                (record.identifier, len(record.sequence)) for record in read_fasta(path)
            )
            if not records:
                raise ValueError(f"routed bin FASTA contains no records: {path}")
            bin_id = path.stem
            for contig_id, _length in records:
                _require_identifier(contig_id, "routed contig ID")
                if contig_id in contig_routes:
                    raise ValueError(
                        f"contig occurs in multiple routed bins: {contig_id}"
                    )
                contig_routes[contig_id] = (route, bin_id)
            sources.append(_RoutedSource(route, bin_id, path, records))
    return sources, contig_routes


def _genomad_candidates(
    viral_dir: Path,
    contig_routes: Mapping[str, tuple[str, str]],
) -> list[_Candidate]:
    summary_path = viral_dir / "genomad/input_summary/input_virus_summary.tsv"
    fasta_path = viral_dir / "genomad/input_summary/input_virus.fna"
    summary_rows = _read_tsv(
        summary_path,
        ("seq_name", "length", "topology", "coordinates"),
    )
    summaries = _index_rows(summary_rows, "seq_name", "geNomad summary")
    records = list(read_fasta(fasta_path))
    record_ids = [record.identifier for record in records]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("geNomad viral FASTA repeats a sequence identifier")
    if set(record_ids) != set(summaries):
        raise ValueError("geNomad summary and viral FASTA identifiers differ")

    candidates: list[_Candidate] = []
    observed_regions: set[tuple[str, int, int]] = set()
    for record in records:
        row = summaries[record.identifier]
        parent, start, end, category = _candidate_coordinates(
            record.identifier, record.sequence, row
        )
        region = (parent, start, end)
        if region in observed_regions:
            raise ValueError(f"geNomad repeats a parent interval: {region}")
        observed_regions.add(region)
        route, bin_id = contig_routes.get(parent, ("unbinned", ""))
        if route in CALLED_ROUTES:
            status = "excluded"
            reason = "parent_assigned_to_called_bin"
        elif route == "unresolved":
            status = "included"
            reason = "unresolved_parent_genomad_candidate"
        else:
            status = "included"
            reason = "unbinned_parent_genomad_candidate"
        candidates.append(
            _Candidate(
                native_id=record.identifier,
                contig_id=_candidate_id(record.identifier, parent, start, end),
                parent_contig=parent,
                start=start,
                end=end,
                sequence=record.sequence,
                bin_id=bin_id,
                category=category,
                status=status,
                reason=reason,
            )
        )
    return candidates


def _candidate_coordinates(
    native_id: str,
    sequence: str,
    summary: TsvRow,
) -> tuple[str, int, int, str]:
    reported_length = _positive_integer(summary["length"], "geNomad length")
    if reported_length != len(sequence):
        raise ValueError(f"geNomad sequence length disagrees for {native_id}")
    match = PROVIRUS_ID.fullmatch(native_id)
    if match is None:
        if summary["coordinates"].strip().casefold() not in MISSING:
            raise ValueError(f"whole-contig geNomad call has coordinates: {native_id}")
        if summary["topology"].strip() == "Provirus":
            raise ValueError(f"geNomad provirus lacks coordinate suffix: {native_id}")
        return native_id, 1, len(sequence), "genomad_full"

    start = int(match.group("start"))
    end = int(match.group("end"))
    if end < start or end - start + 1 != len(sequence):
        raise ValueError(f"geNomad provirus coordinates disagree for {native_id}")
    if summary["coordinates"].strip() != f"{start}-{end}":
        raise ValueError(f"geNomad summary coordinates disagree for {native_id}")
    if summary["topology"].strip() != "Provirus":
        raise ValueError(f"geNomad provirus topology disagrees for {native_id}")
    return match.group("parent"), start, end, "genomad_provirus"


def _validate_candidate_intervals(candidates: list[_Candidate]) -> None:
    by_parent: dict[str, list[_Candidate]] = {}
    for candidate in candidates:
        if candidate.status == "included":
            by_parent.setdefault(candidate.parent_contig, []).append(candidate)
    for parent, records in by_parent.items():
        ordered = sorted(records, key=lambda record: (record.start, record.end))
        for previous, current in pairwise(ordered):
            if current.start <= previous.end:
                raise ValueError(f"geNomad candidate intervals overlap on {parent}")


def _stage_routed_source(
    source: _RoutedSource,
    output_dir: Path,
    sources_dir: Path,
    origin_maps_dir: Path,
) -> tuple[TsvRow, list[TsvRow]]:
    category = f"{source.route}_bin"
    source_id = f"{source.route}-bin-{_safe_component(source.bin_id)}"
    fasta = sources_dir / f"{source_id}.fna"
    shutil.copyfile(source.path, fasta)
    rows = [
        {
            "contig_id": contig_id,
            "parent_contig": contig_id,
            "start": "1",
            "end": str(length),
            "source_id": source_id,
            "bin_id": source.bin_id,
            "category": category,
            "status": "included",
            "reason": f"routed_{source.route}_bin",
        }
        for contig_id, length in source.contigs
    ]
    origin_map = origin_maps_dir / f"{source_id}.tsv"
    _write_tsv(origin_map, ORIGIN_FIELDS, rows)
    return (
        _source_row(
            source_id,
            category,
            source.bin_id,
            fasta.relative_to(output_dir),
            origin_map.relative_to(output_dir),
            f"routed_{source.route}_bin",
        ),
        rows,
    )


def _stage_genomad_source(
    candidates: list[_Candidate],
    output_dir: Path,
    sources_dir: Path,
    origin_maps_dir: Path,
) -> tuple[TsvRow, list[TsvRow]]:
    fasta = sources_dir / f"{GENOMAD_SOURCE}.fna"
    with fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for candidate in candidates:
            handle.write(f">{candidate.contig_id}\n")
            for offset in range(0, len(candidate.sequence), 60):
                handle.write(f"{candidate.sequence[offset : offset + 60]}\n")
    rows = [_origin_row(candidate) for candidate in candidates]
    origin_map = origin_maps_dir / f"{GENOMAD_SOURCE}.tsv"
    _write_tsv(origin_map, ORIGIN_FIELDS, rows)
    return (
        _source_row(
            GENOMAD_SOURCE,
            "genomad_unbinned",
            "",
            fasta.relative_to(output_dir),
            origin_map.relative_to(output_dir),
            "unassigned_parent_genomad_candidates",
        ),
        rows,
    )


def _source_row(
    source_id: str,
    category: str,
    bin_id: str,
    fasta: Path,
    origin_map: Path,
    reason: str,
) -> TsvRow:
    return {
        "source_id": source_id,
        "category": category,
        "bin_id": bin_id,
        "fasta": fasta.as_posix(),
        "origin_map": origin_map.as_posix(),
        "status": "included",
        "reason": reason,
    }


def _origin_row(candidate: _Candidate) -> TsvRow:
    return {
        "contig_id": candidate.contig_id,
        "parent_contig": candidate.parent_contig,
        "start": str(candidate.start),
        "end": str(candidate.end),
        "source_id": GENOMAD_SOURCE,
        "bin_id": candidate.bin_id,
        "category": candidate.category,
        "status": candidate.status,
        "reason": candidate.reason,
    }


def _candidate_id(native_id: str, parent: str, start: int, end: int) -> str:
    digest = hashlib.sha256(native_id.encode()).hexdigest()[:10]
    return f"genomad_{_safe_component(parent)}_{start}_{end}_{digest}"


def _safe_component(value: str) -> str:
    if IDENTIFIER.fullmatch(value):
        return value
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.-")
    if not cleaned or not cleaned[0].isalnum():
        cleaned = f"id_{cleaned}"
    digest = hashlib.sha256(value.encode()).hexdigest()[:10]
    return f"{cleaned}-{digest}"


def _index_rows(rows: list[TsvRow], key: str, label: str) -> dict[str, TsvRow]:
    indexed: dict[str, TsvRow] = {}
    for row in rows:
        identifier = row[key].strip()
        if not identifier:
            raise ValueError(f"{label} contains an empty identifier")
        if identifier in indexed:
            raise ValueError(f"{label} repeats identifier: {identifier}")
        indexed[identifier] = row
    return indexed


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


def _write_tsv(path: Path, fields: tuple[str, ...], rows: list[TsvRow]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _positive_integer(value: str, field: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValueError(f"{field} must be an integer: {value!r}") from error
    if parsed < 1:
        raise ValueError(f"{field} must be at least 1: {value!r}")
    return parsed


def _require_identifier(value: str, field: str) -> None:
    if not IDENTIFIER.fullmatch(value):
        raise ValueError(f"invalid {field}: {value!r}")
