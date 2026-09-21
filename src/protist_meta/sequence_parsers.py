"""Normalize SSUextract, geNomad, and CheckV native result tables."""

from __future__ import annotations

import csv
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "ParsedSequenceRows",
    "parse_checkv",
    "parse_genomad",
    "parse_ssu",
]

IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
COORDINATES = re.compile(r"^([1-9][0-9]*)-([1-9][0-9]*)$")
PROVIRUS_ID = re.compile(
    r"^(?P<parent>[A-Za-z0-9][A-Za-z0-9_.:-]*)"
    r"\|provirus_(?P<start>[1-9][0-9]*)_(?P<end>[1-9][0-9]*)$"
)
MISSING_VALUES = {"", "na", "n/a", "nan", "none", "nd", "unavailable"}
SSU_MODELS = {"RF00177": "rrna_16s", "RF01960": "rrna_18s"}
VIRUS_RANKS = (
    "root",
    "realm",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)

type NormalizedRow = dict[str, object]


@dataclass(frozen=True, slots=True)
class ParsedSequenceRows:
    """Catalog-ready record bodies and native evidence outside the schema."""

    ssu: list[NormalizedRow]
    viruses: list[NormalizedRow]
    qc: list[NormalizedRow]
    taxonomy: list[NormalizedRow]
    target_ids: tuple[str, ...]
    evidence_gaps: tuple[str, ...] = ()


def parse_ssu(path: Path) -> ParsedSequenceRows:
    """Parse current SSUextract loci with one-based inclusive coordinates."""
    required = {
        "name",
        "model",
        "length",
        "coordinates",
        "strand",
        "contig_name",
        "blast_sseqid",
        "blast_pident",
        "taxonomy",
        "reference_identifiers",
        "reference_source",
        "reference_versions",
        "reference_taxonomy",
        "taxonomy_assignment_method",
        "blast_query_coverage",
    }
    rows = _read_tsv(path, required)
    loci: list[NormalizedRow] = []
    target_ids: list[str] = []
    observed: set[str] = set()
    locus_ids: set[str] = set()
    for row in rows:
        native_id = _unique_native_id(row["name"], observed, "SSUextract")
        target_ids.append(native_id)
        locus = _ssu_locus(row)
        locus_id = str(locus["record_id"])
        if locus_id in locus_ids:
            raise ValueError(f"duplicate SSUextract locus: {locus_id}")
        locus_ids.add(locus_id)
        loci.append(locus)
    return ParsedSequenceRows(
        ssu=loci,
        viruses=[],
        qc=[],
        taxonomy=[],
        target_ids=tuple(target_ids),
        evidence_gaps=(
            (
                "SSUextract model coverage, fragment assembly, alternative calls, tied "
                "hits, and raw sequence remain in the source report because the "
                "catalog has no matching fields."
            ),
        ),
    )


def parse_genomad(
    path: Path,
    contig_lengths: Mapping[str, int],
) -> ParsedSequenceRows:
    """Parse geNomad calls against exact parent-contig lengths."""
    required = {
        "seq_name",
        "length",
        "topology",
        "coordinates",
        "n_genes",
        "genetic_code",
        "virus_score",
        "fdr",
        "n_hallmarks",
        "marker_enrichment",
        "taxonomy",
    }
    rows = _read_tsv(path, required)
    viruses: list[NormalizedRow] = []
    taxonomy: list[NormalizedRow] = []
    target_ids: list[str] = []
    observed: set[str] = set()
    for row in rows:
        native_id = _unique_native_id(row["seq_name"], observed, "geNomad")
        target_ids.append(native_id)
        virus = _genomad_virus(row, contig_lengths)
        viruses.append(virus)
        taxonomy.extend(
            _virus_taxonomy(
                native_id,
                str(virus["record_id"]),
                row["taxonomy"],
            )
        )
    return ParsedSequenceRows(
        ssu=[],
        viruses=viruses,
        qc=[],
        taxonomy=taxonomy,
        target_ids=tuple(target_ids),
        evidence_gaps=(
            (
                "geNomad topology, sequence length, gene count, genetic code, FDR, "
                "hallmark count, and marker enrichment remain in the source report "
                "because the catalog has no matching fields."
            ),
        ),
    )


def parse_checkv(
    path: Path,
    virus_ids: Mapping[str, str],
) -> ParsedSequenceRows:
    """Parse CheckV QC and bind every native target to a geNomad virus ID."""
    required = {
        "contig_id",
        "contig_length",
        "provirus",
        "proviral_length",
        "gene_count",
        "viral_genes",
        "host_genes",
        "checkv_quality",
        "miuvig_quality",
        "completeness",
        "completeness_method",
        "contamination",
        "kmer_freq",
        "warnings",
    }
    _validate_virus_ids(virus_ids)
    rows = _read_tsv(path, required)
    qc: list[NormalizedRow] = []
    target_ids: list[str] = []
    observed: set[str] = set()
    for row in rows:
        native_id = _unique_native_id(row["contig_id"], observed, "CheckV")
        target_ids.append(native_id)
        if native_id not in virus_ids:
            raise ValueError(f"CheckV target has no geNomad virus mapping: {native_id}")
        assessment = _checkv_qc(row, virus_ids[native_id])
        if assessment is not None:
            qc.append(assessment)
    if observed != set(virus_ids):
        raise ValueError("CheckV and geNomad virus target sets differ")
    return ParsedSequenceRows(
        ssu=[],
        viruses=[],
        qc=qc,
        taxonomy=[],
        target_ids=tuple(target_ids),
        evidence_gaps=(
            (
                "CheckV quality tiers, provirus flag and length, gene counts, "
                "k-mer frequency, and warnings remain in the "
                "source report because the catalog has no matching fields."
            ),
            (
                "CheckV quality_summary.tsv does not provide parent-contig "
                "coordinates for newly trimmed proviruses."
            ),
        ),
    )


def _ssu_locus(row: dict[str, str]) -> NormalizedRow:
    contig_id = _identifier(row["contig_name"], "SSUextract contig")
    start, end = _coordinates(row["coordinates"], "SSUextract")
    length = _required_integer(row["length"], "SSUextract length", minimum=1)
    if end - start + 1 != length:
        raise ValueError(
            f"SSUextract length disagrees with coordinates for {contig_id}"
        )
    model = row["model"].strip()
    if model not in SSU_MODELS:
        raise ValueError(f"unsupported SSUextract model: {model!r}")
    strand = {"+": "plus", "-": "minus"}.get(row["strand"].strip())
    if strand is None:
        raise ValueError(f"invalid SSUextract strand for {contig_id}")
    locus: NormalizedRow = {
        "record_id": f"ssuextract:ssu:{contig_id}:{start}-{end}:{model}",
        "contig_id": contig_id,
        "ssu_type": SSU_MODELS[model],
        "start": start,
        "end": end,
        "strand": strand,
    }
    _set_nonspace_text(locus, "hit_accession", row["blast_sseqid"])
    _set_text(locus, "reference_identifiers", row["reference_identifiers"])
    _set_text(locus, "reference_source", row["reference_source"])
    _set_text(locus, "reference_versions", row["reference_versions"])
    _set_text(locus, "reference_taxonomy", row["reference_taxonomy"])
    _set_text(locus, "taxon_name", row["taxonomy"])
    _set_text(
        locus,
        "taxonomy_assignment_method",
        row["taxonomy_assignment_method"],
    )
    identity = _float(row["blast_pident"], "blast_pident", 0, 100)
    coverage = _float(row["blast_query_coverage"], "blast_query_coverage", 0, 1)
    if (identity is not None or coverage is not None) and "hit_accession" not in locus:
        raise ValueError(
            f"SSUextract BLAST evidence lacks an accession for {contig_id}"
        )
    _set_number(
        locus,
        "identity_percent",
        identity,
    )
    _set_number(
        locus,
        "query_coverage_percent",
        coverage * 100 if coverage is not None else None,
    )
    return locus


def _genomad_virus(
    row: dict[str, str],
    contig_lengths: Mapping[str, int],
) -> NormalizedRow:
    native_id = row["seq_name"].strip()
    _validate_genomad_evidence(row)
    length = _required_integer(row["length"], "geNomad length", minimum=1)
    parent, start, end, is_provirus = _genomad_coordinates(row, contig_lengths)
    if end - start + 1 != length:
        raise ValueError(f"geNomad length disagrees with coordinates for {native_id}")
    score = _required_float(row["virus_score"], "virus_score", 0, 1)
    record_id = _virus_record_id(parent, start, end, is_provirus=is_provirus)
    virus: NormalizedRow = {
        "record_id": record_id,
        "contig_id": parent,
        "start": start,
        "end": end,
        "strand": "unknown",
        "score": score,
    }
    _set_text(virus, "viral_taxon", row["taxonomy"])
    return virus


def _validate_genomad_evidence(row: dict[str, str]) -> None:
    _required_integer(row["n_genes"], "n_genes", minimum=0)
    _required_integer(row["genetic_code"], "genetic_code", minimum=1)
    _float(row["fdr"], "fdr", 0, 1)
    _required_integer(row["n_hallmarks"], "n_hallmarks", minimum=0)
    _float(row["marker_enrichment"], "marker_enrichment", None, None)
    if _optional_text(row["topology"]) is None:
        raise ValueError("geNomad topology is missing")
    if _optional_text(row["taxonomy"]) is None:
        raise ValueError("geNomad taxonomy is missing")


def _genomad_coordinates(
    row: dict[str, str],
    contig_lengths: Mapping[str, int],
) -> tuple[str, int, int, bool]:
    native_id = row["seq_name"].strip()
    provirus = PROVIRUS_ID.fullmatch(native_id)
    if provirus is not None:
        parent = _known_contig(provirus.group("parent"), contig_lengths)
        name_coordinates = (int(provirus.group("start")), int(provirus.group("end")))
        table_coordinates = _coordinates(row["coordinates"], "geNomad provirus")
        if name_coordinates != table_coordinates:
            raise ValueError(f"geNomad provirus coordinates disagree for {native_id}")
        if row["topology"].strip() != "Provirus":
            raise ValueError(f"geNomad provirus topology disagrees for {native_id}")
        _check_coordinate_bounds(parent, *table_coordinates, contig_lengths[parent])
        return parent, *table_coordinates, True

    parent = _known_contig(native_id, contig_lengths)
    if not _is_missing(row["coordinates"]):
        raise ValueError(f"geNomad whole-contig call has coordinates: {native_id}")
    if row["topology"].strip() == "Provirus":
        raise ValueError(
            f"geNomad provirus lacks an exact coordinate suffix: {native_id}"
        )
    return parent, 1, contig_lengths[parent], False


def _virus_record_id(
    parent: str,
    start: int,
    end: int,
    *,
    is_provirus: bool,
) -> str:
    if is_provirus:
        return f"genomad:virus:{parent}:{start}-{end}"
    return f"genomad:virus:{parent}"


def _virus_taxonomy(
    native_id: str,
    virus_id: str,
    lineage: str,
) -> list[NormalizedRow]:
    text = lineage.strip()
    if _is_missing(text) or text.casefold() == "unclassified":
        return []
    taxa = text.split(";")
    if len(taxa) > len(VIRUS_RANKS):
        raise ValueError(f"geNomad taxonomy has too many ranks for {native_id}")
    if taxa[0].strip() != "Viruses":
        raise ValueError(f"geNomad taxonomy lacks the Viruses root for {native_id}")
    records: list[NormalizedRow] = []
    for rank, taxon in zip(VIRUS_RANKS, taxa, strict=False):
        name = taxon.strip()
        if name == "":
            continue
        records.append(
            {
                "record_id": f"genomad:taxonomy:{virus_id}:{rank}",
                "target_kind": "virus",
                "target_id": virus_id,
                "rank": rank,
                "taxon_name": name,
            }
        )
    return records


def _checkv_qc(row: dict[str, str], virus_id: str) -> NormalizedRow | None:
    _validate_checkv_evidence(row, virus_id)
    completeness = _float(row["completeness"], "completeness", 0, 100)
    contamination = _float(row["contamination"], "contamination", 0, 100)
    method = _optional_text(row["completeness_method"])
    if completeness is not None and method is None:
        raise ValueError(f"CheckV completeness lacks a method for {virus_id}")
    if completeness is None and method is not None:
        raise ValueError(f"CheckV method lacks completeness for {virus_id}")
    if completeness is None and contamination is None:
        return None
    qc: NormalizedRow = {
        "record_id": f"checkv:qc:{virus_id}",
        "target_kind": "virus",
        "target_id": virus_id,
    }
    _set_number(qc, "completeness_percent", completeness)
    if method is not None:
        qc["completeness_basis"] = method
    _set_number(qc, "contamination_percent", contamination)
    if contamination is not None:
        qc["contamination_basis"] = "CheckV predicted host-region fraction"
    return qc


def _validate_checkv_evidence(row: dict[str, str], virus_id: str) -> None:
    contig_length = _required_integer(row["contig_length"], "contig_length", minimum=1)
    is_provirus = _boolean_word(row["provirus"], "provirus")
    proviral_length = _optional_integer(
        row["proviral_length"], "proviral_length", minimum=1
    )
    if is_provirus and proviral_length is None:
        raise ValueError(f"CheckV provirus lacks a length for {virus_id}")
    if not is_provirus and proviral_length is not None:
        raise ValueError(f"CheckV non-provirus has a proviral length for {virus_id}")
    if proviral_length is not None and proviral_length > contig_length:
        raise ValueError(f"CheckV proviral length exceeds contig for {virus_id}")
    gene_count = _required_integer(row["gene_count"], "gene_count", minimum=0)
    viral_genes = _required_integer(row["viral_genes"], "viral_genes", minimum=0)
    host_genes = _required_integer(row["host_genes"], "host_genes", minimum=0)
    if viral_genes + host_genes > gene_count:
        raise ValueError(f"CheckV gene categories exceed gene count for {virus_id}")
    if _optional_text(row["checkv_quality"]) is None:
        raise ValueError(f"CheckV quality tier is missing for {virus_id}")
    if _optional_text(row["miuvig_quality"]) is None:
        raise ValueError(f"CheckV MIUViG tier is missing for {virus_id}")
    _float(row["kmer_freq"], "kmer_freq", 0, None)


def _read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if fields is None:
            raise ValueError(f"TSV has no header: {path}")
        if len(fields) != len(set(fields)):
            raise ValueError(f"TSV has duplicate columns: {path}")
        missing = sorted(required - set(fields))
        if missing:
            raise ValueError(f"TSV is missing required columns {missing}: {path}")
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, start=2):
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError(f"malformed TSV row at line {line_number}: {path}")
            rows.append({key: value for key, value in raw.items() if key is not None})
    return rows


def _unique_native_id(value: str, observed: set[str], tool: str) -> str:
    native_id = value.strip()
    if native_id == "" or any(character.isspace() for character in native_id):
        raise ValueError(f"invalid {tool} target identifier: {value!r}")
    if native_id in observed:
        raise ValueError(f"duplicate {tool} target: {native_id}")
    observed.add(native_id)
    return native_id


def _identifier(value: str, field: str) -> str:
    identifier = value.strip()
    if IDENTIFIER.fullmatch(identifier) is None:
        raise ValueError(f"invalid {field} identifier: {value!r}")
    return identifier


def _known_contig(value: str, contig_lengths: Mapping[str, int]) -> str:
    contig_id = _identifier(value, "geNomad parent contig")
    if contig_id not in contig_lengths:
        raise ValueError(f"geNomad parent contig is unknown: {contig_id}")
    if contig_lengths[contig_id] < 1:
        raise ValueError(f"geNomad parent contig length is invalid: {contig_id}")
    return contig_id


def _coordinates(value: str, field: str) -> tuple[int, int]:
    match = COORDINATES.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"invalid {field} coordinates: {value!r}")
    start, end = int(match.group(1)), int(match.group(2))
    if start > end:
        raise ValueError(f"{field} coordinate start exceeds end")
    return start, end


def _check_coordinate_bounds(
    contig_id: str,
    start: int,
    end: int,
    contig_length: int,
) -> None:
    if start < 1 or end > contig_length:
        raise ValueError(f"geNomad coordinates exceed parent contig {contig_id}")


def _float(
    value: str,
    field: str,
    minimum: float | None,
    maximum: float | None,
) -> float | None:
    text = value.strip()
    if _is_missing(text):
        return None
    try:
        number = float(text)
    except ValueError as error:
        raise ValueError(f"{field} is not numeric") from error
    if not math.isfinite(number):
        raise ValueError(f"{field} is not finite")
    if minimum is not None and number < minimum:
        raise ValueError(f"{field} is below {minimum}")
    if maximum is not None and number > maximum:
        raise ValueError(f"{field} exceeds {maximum}")
    return number


def _required_float(
    value: str,
    field: str,
    minimum: float,
    maximum: float | None,
) -> float:
    number = _float(value, field, minimum, maximum)
    if number is None:
        raise ValueError(f"{field} is missing")
    return number


def _optional_integer(
    value: str,
    field: str,
    *,
    minimum: int,
) -> int | None:
    text = value.strip()
    if _is_missing(text):
        return None
    if not text.isdigit():
        raise ValueError(f"{field} is not an unsigned integer")
    number = int(text)
    if number < minimum:
        raise ValueError(f"{field} is below {minimum}")
    return number


def _required_integer(
    value: str,
    field: str,
    *,
    minimum: int,
) -> int:
    number = _optional_integer(value, field, minimum=minimum)
    if number is None:
        raise ValueError(f"{field} is missing")
    return number


def _boolean_word(value: str, field: str) -> bool:
    normalized = value.strip().casefold()
    if normalized == "yes":
        return True
    if normalized == "no":
        return False
    raise ValueError(f"{field} must be Yes or No")


def _validate_virus_ids(virus_ids: Mapping[str, str]) -> None:
    normalized = [_identifier(value, "geNomad virus") for value in virus_ids.values()]
    if len(normalized) != len(set(normalized)):
        raise ValueError("geNomad virus mapping contains duplicate record IDs")
    for native_id in virus_ids:
        if (
            native_id == ""
            or native_id.strip() != native_id
            or any(character.isspace() for character in native_id)
        ):
            raise ValueError(f"invalid geNomad target identifier: {native_id!r}")


def _is_missing(value: str) -> bool:
    return value.strip().casefold() in MISSING_VALUES


def _optional_text(value: str) -> str | None:
    text = value.strip()
    return None if _is_missing(text) else text


def _set_text(row: NormalizedRow, field: str, value: str) -> None:
    text = _optional_text(value)
    if text is not None:
        row[field] = text


def _set_nonspace_text(row: NormalizedRow, field: str, value: str) -> None:
    text = _optional_text(value)
    if text is None:
        return
    if any(character.isspace() for character in text):
        raise ValueError(f"{field} contains whitespace")
    row[field] = text


def _set_number(row: NormalizedRow, field: str, value: float | None) -> None:
    if value is not None:
        row[field] = value
