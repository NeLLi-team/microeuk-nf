"""Parse native gene calls and prepare uniquely named annotation proteins."""

from __future__ import annotations

import csv
import hashlib
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast
from urllib.parse import unquote

__all__ = [
    "GeneRow",
    "ProteinMapping",
    "ProteinSource",
    "merge_annotation_proteins",
    "parse_braker3",
    "parse_prodigal_gv",
]

type Caller = Literal["braker3", "prodigal-gv"]
type GeneRow = dict[str, object]

IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
BRAKER_TRANSCRIPT = re.compile(r"^(.+)\.t[0-9]+$")
TRANSLATION_TABLE = re.compile(r"(?:^|;)transl_table=([0-9]+)(?:;|$)")


@dataclass(frozen=True, slots=True)
class ProteinSource:
    """One caller- and source-scoped protein FASTA to merge."""

    fasta: Path
    caller: Caller
    source_id: str
    bin_id: str | None = None


@dataclass(frozen=True, slots=True)
class ProteinMapping:
    """Relationship between an annotation identifier and its native call."""

    normalized_id: str
    native_id: str
    native_gene_id: str
    caller: Caller
    source_id: str
    bin_id: str | None


@dataclass(frozen=True, slots=True)
class _FastaRecord:
    identifier: str
    header: str
    sequence: str


@dataclass(frozen=True, slots=True)
class _Feature:
    contig_id: str
    feature_type: str
    start: int
    end: int
    strand: str
    phase: int | None
    attributes: dict[str, str]
    model_code: int | None


@dataclass(frozen=True, slots=True)
class _ProdigalProtein:
    record: _FastaRecord
    gene_id: str
    start: int
    end: int
    strand: str
    genetic_code: int


@dataclass(frozen=True, slots=True)
class _BrakerFeatures:
    genes: dict[str, _Feature]
    transcripts: dict[str, _Feature]
    transcript_cds: dict[str, list[_Feature]]


def parse_prodigal_gv(
    gff: Path,
    proteins: Path,
    *,
    source_id: str,
    bin_id: str | None,
    genetic_code: int | None = None,
) -> list[GeneRow]:
    """Parse Prodigal-GV CDS calls, joined to proteins by the native gene ID.

    Prodigal-GV can select a different translation table for each sequence.
    With ``genetic_code=None``, the per-CDS ``genetic_code`` attribute is
    retained. An explicit value is a validation constraint, not an override.

    Args:
        gff: Native Prodigal-GV GFF3 output.
        proteins: Native Prodigal-GV amino-acid FASTA output.
        source_id: Stable source namespace, such as an assembly or bin ID.
        bin_id: Catalog bin ID, or ``None`` for assembly-derived viral calls.
        genetic_code: Optional required NCBI translation table.

    Returns:
        One catalog-ready gene row per protein.

    Raises:
        ValueError: If native records are inconsistent or duplicated.
    """
    required_code = _optional_genetic_code(genetic_code)
    features = _read_gff(gff, allow_empty=True)
    protein_records = _read_fasta(proteins, allow_empty=True)
    if not features:
        _validate_empty_prodigal_gff(gff, required_code)
        if protein_records:
            raise ValueError("Prodigal-GV proteins exist without GFF CDS records")
        return []

    cds_by_id = _index_prodigal_cds(features, gff, required_code)
    parsed_proteins = [_parse_prodigal_header(record) for record in protein_records]
    protein_gene_ids = {protein.gene_id for protein in parsed_proteins}
    if len(protein_gene_ids) != len(parsed_proteins):
        raise ValueError("duplicate Prodigal-GV gene ID in protein headers")
    if protein_gene_ids != set(cds_by_id):
        raise ValueError(
            _set_mismatch("Prodigal-GV GFF CDS", cds_by_id, protein_gene_ids)
        )

    return [
        _prodigal_row(
            protein,
            cds_by_id[protein.gene_id],
            gff,
            source_id,
            bin_id,
        )
        for protein in parsed_proteins
    ]


def parse_braker3(
    gff3: Path,
    proteins: Path,
    *,
    source_id: str,
    bin_id: str | None,
    genetic_code: int = 1,
) -> list[GeneRow]:
    """Parse BRAKER3 transcripts and their translated proteins.

    Each protein-producing transcript becomes one gene row. This preserves
    alternative isoforms instead of collapsing them to their parent locus.
    Coordinates span the complete transcript, including introns, and phase is
    taken from the first CDS in transcript order.

    Args:
        gff3: Native BRAKER3 GFF3 output.
        proteins: Native BRAKER3 amino-acid FASTA output.
        source_id: Stable source namespace, usually the input bin ID.
        bin_id: Catalog bin ID, or ``None`` for an unbinned source.
        genetic_code: Explicit eukaryotic NCBI translation table.

    Returns:
        One catalog-ready gene row per translated transcript.

    Raises:
        ValueError: If relationships, coordinates, or protein IDs disagree.
    """
    code = _genetic_code(genetic_code)
    indexed = _index_braker_features(_read_gff(gff3), gff3)
    protein_records = _read_fasta(proteins)
    protein_ids = {record.identifier for record in protein_records}
    if protein_ids != set(indexed.transcripts):
        raise ValueError(
            _set_mismatch("BRAKER3 transcripts", indexed.transcripts, protein_ids)
        )

    return [
        _braker_row(record, indexed, gff3, source_id, bin_id, code)
        for record in protein_records
    ]


def _index_braker_features(
    features: list[_Feature],
    gff3: Path,
) -> _BrakerFeatures:
    genes: dict[str, _Feature] = {}
    transcripts: dict[str, _Feature] = {}
    transcript_cds: dict[str, list[_Feature]] = {}
    for feature in features:
        if feature.feature_type == "gene":
            gene_id = _required_attribute(feature, "ID", gff3)
            _add_unique(genes, gene_id, feature, "BRAKER3 gene")
        elif feature.feature_type in {"mRNA", "transcript"}:
            transcript_id = _required_attribute(feature, "ID", gff3)
            _add_unique(transcripts, transcript_id, feature, "BRAKER3 transcript")
        elif feature.feature_type == "CDS":
            parents = _parents(feature, gff3)
            for transcript_id in parents:
                transcript_cds.setdefault(transcript_id, []).append(feature)

    if not transcripts:
        raise ValueError(f"BRAKER3 GFF contains no transcript records: {gff3}")
    return _BrakerFeatures(genes, transcripts, transcript_cds)


def _braker_row(
    protein: _FastaRecord,
    indexed: _BrakerFeatures,
    gff3: Path,
    source_id: str,
    bin_id: str | None,
    genetic_code: int,
) -> GeneRow:
    transcript = indexed.transcripts[protein.identifier]
    parents = _parents(transcript, gff3)
    if len(parents) != 1:
        message = f"BRAKER3 transcript {protein.identifier} must have one parent gene"
        raise ValueError(message)
    parent_id = parents[0]
    if parent_id not in indexed.genes:
        message = (
            f"BRAKER3 transcript {protein.identifier} has unknown parent {parent_id}"
        )
        raise ValueError(message)
    _validate_contained(
        transcript, indexed.genes[parent_id], protein.identifier, parent_id
    )
    cds_records = indexed.transcript_cds.get(protein.identifier, [])
    if not cds_records:
        raise ValueError(f"BRAKER3 transcript has no CDS: {protein.identifier}")
    for cds in cds_records:
        _validate_contained(cds, transcript, "CDS", protein.identifier)
        if cds.phase is None:
            raise ValueError(f"BRAKER3 CDS phase is missing for {protein.identifier}")
    return _gene_row(
        gene_id=_normalized_id("braker3", source_id, bin_id, protein.identifier),
        feature=transcript,
        bin_id=bin_id,
        phase=_first_cds(cds_records, transcript.strand).phase,
        genetic_code=genetic_code,
    )


def _index_prodigal_cds(
    features: list[_Feature],
    gff: Path,
    required_code: int | None,
) -> dict[str, _Feature]:
    cds_by_id: dict[str, _Feature] = {}
    for feature in features:
        if feature.feature_type != "CDS":
            continue
        native_gene_id = _required_attribute(feature, "ID", gff)
        if native_gene_id in cds_by_id:
            raise ValueError(f"duplicate Prodigal-GV CDS ID: {native_gene_id}")
        code = _integer_attribute(feature, "genetic_code", gff)
        if feature.model_code is not None and code != feature.model_code:
            message = (
                f"Prodigal-GV CDS {native_gene_id} uses genetic code {code}, "
                f"but its model declares {feature.model_code}"
            )
            raise ValueError(message)
        if required_code is not None and code != required_code:
            message = (
                f"Prodigal-GV CDS {native_gene_id} uses genetic code {code}, "
                f"expected {required_code}"
            )
            raise ValueError(message)
        cds_by_id[native_gene_id] = feature
    if not cds_by_id:
        raise ValueError(f"Prodigal-GV GFF contains no CDS records: {gff}")
    return cds_by_id


def _prodigal_row(
    protein: _ProdigalProtein,
    feature: _Feature,
    gff: Path,
    source_id: str,
    bin_id: str | None,
) -> GeneRow:
    code = _integer_attribute(feature, "genetic_code", gff)
    expected_id = _prodigal_protein_id(feature.contig_id, protein.gene_id)
    if protein.record.identifier != expected_id:
        message = (
            f"Prodigal-GV protein {protein.record.identifier} does not match "
            f"GFF contig/CDS {feature.contig_id}/{protein.gene_id}"
        )
        raise ValueError(message)
    if (
        protein.start != feature.start
        or protein.end != feature.end
        or protein.strand != feature.strand
        or protein.genetic_code != code
    ):
        message = f"Prodigal-GV protein header disagrees with GFF CDS {protein.gene_id}"
        raise ValueError(message)
    return _gene_row(
        gene_id=_normalized_id(
            "prodigal-gv", source_id, bin_id, protein.record.identifier
        ),
        feature=feature,
        bin_id=bin_id,
        phase=feature.phase,
        genetic_code=code,
    )


def merge_annotation_proteins(
    inputs: Sequence[ProteinSource],
    output_fasta: Path,
    mapping_tsv: Path,
) -> list[ProteinMapping]:
    """Merge protein FASTAs and normalize one terminal stop for annotation.

    Args:
        inputs: Native protein FASTAs and their provenance.
        output_fasta: Destination FASTA for functional annotation.
        mapping_tsv: Destination identifier/provenance mapping.

    Returns:
        Mapping rows in output FASTA order.

    Raises:
        ValueError: If an input is missing, malformed, contains an invalid stop
            marker, or has a normalized ID collision.
    """
    if not inputs:
        raise ValueError("at least one protein source is required")
    if output_fasta.resolve() == mapping_tsv.resolve():
        raise ValueError("protein FASTA and mapping TSV paths must differ")

    merged: list[tuple[ProteinMapping, str]] = []
    observed: set[str] = set()
    for source in inputs:
        _caller(source.caller)
        records = _read_fasta(source.fasta, allow_empty=True)
        for record in records:
            sequence = record.sequence.removesuffix("*")
            if not sequence or "*" in sequence:
                message = (
                    f"invalid annotation protein {record.identifier!r} from "
                    f"{source.fasta}: expected at most one terminal stop marker"
                )
                raise ValueError(message)
            normalized_id = _normalized_id(
                source.caller, source.source_id, source.bin_id, record.identifier
            )
            if normalized_id in observed:
                raise ValueError(f"duplicate normalized protein ID: {normalized_id}")
            observed.add(normalized_id)
            mapping = ProteinMapping(
                normalized_id=normalized_id,
                native_id=record.identifier,
                native_gene_id=_native_gene_id(source.caller, record),
                caller=source.caller,
                source_id=source.source_id,
                bin_id=source.bin_id,
            )
            merged.append((mapping, sequence))
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    mapping_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for mapping, sequence in merged:
            handle.write(f">{mapping.normalized_id}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(f"{sequence[offset : offset + 60]}\n")
    with mapping_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "normalized_id",
                "native_id",
                "native_gene_id",
                "caller",
                "source_id",
                "bin_id",
            ]
        )
        for mapping, _sequence in merged:
            writer.writerow(
                [
                    mapping.normalized_id,
                    mapping.native_id,
                    mapping.native_gene_id,
                    mapping.caller,
                    mapping.source_id,
                    mapping.bin_id or "",
                ]
            )
    return [mapping for mapping, _sequence in merged]


def _read_fasta(path: Path, *, allow_empty: bool = False) -> list[_FastaRecord]:
    records: list[_FastaRecord] = []
    header: str | None = None
    sequence: list[str] = []
    observed: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(_fasta_record(header, sequence, path, observed))
                header = line[1:].strip()
                sequence = []
                if not header:
                    raise ValueError(f"empty FASTA header at {path}:{line_number}")
            elif header is None:
                raise ValueError(
                    f"sequence precedes FASTA header at {path}:{line_number}"
                )
            else:
                sequence.append("".join(line.split()))
    if header is not None:
        records.append(_fasta_record(header, sequence, path, observed))
    if not records and not allow_empty:
        raise ValueError(f"protein FASTA contains no records: {path}")
    return records


def _fasta_record(
    header: str,
    sequence: list[str],
    path: Path,
    observed: set[str],
) -> _FastaRecord:
    identifier = header.split(maxsplit=1)[0]
    residues = "".join(sequence)
    if not identifier or not residues:
        raise ValueError(f"empty FASTA identifier or sequence in {path}")
    if identifier in observed:
        raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
    observed.add(identifier)
    return _FastaRecord(identifier, header, residues)


def _read_gff(path: Path, *, allow_empty: bool = False) -> list[_Feature]:
    features: list[_Feature] = []
    model_code: int | None = None
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if line.startswith("##FASTA"):
                break
            if line.startswith("# Model Data:"):
                match = TRANSLATION_TABLE.search(
                    line.removeprefix("# Model Data:").strip()
                )
                model_code = int(match.group(1)) if match else None
                continue
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) != 9:
                raise ValueError(f"expected 9 GFF columns at {path}:{line_number}")
            (
                contig_id,
                _source,
                feature_type,
                start,
                end,
                _score,
                strand,
                phase,
                attrs,
            ) = columns
            _require_identifier(contig_id, "GFF contig ID")
            parsed_start = _positive_integer(start, "GFF start")
            parsed_end = _positive_integer(end, "GFF end")
            if parsed_end < parsed_start:
                raise ValueError(f"GFF end precedes start at {path}:{line_number}")
            features.append(
                _Feature(
                    contig_id=contig_id,
                    feature_type=feature_type,
                    start=parsed_start,
                    end=parsed_end,
                    strand=_strand(strand, path, line_number),
                    phase=_phase(phase, path, line_number),
                    attributes=_attributes(attrs, path, line_number),
                    model_code=model_code,
                )
            )
    if not features and not allow_empty:
        raise ValueError(f"GFF contains no feature records: {path}")
    return features


def _validate_empty_prodigal_gff(path: Path, required_code: int | None) -> None:
    has_gff3_header = False
    has_sequence = False
    has_model = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("##gff-version"):
                has_gff3_header = line.split()[-1] == "3"
            elif line.startswith("# Sequence Data:"):
                has_sequence = True
            elif line.startswith("# Model Data:"):
                model = line.removeprefix("# Model Data:").strip()
                match = TRANSLATION_TABLE.search(model)
                if match is None:
                    raise ValueError(f"Prodigal-GV model lacks transl_table: {path}")
                model_code = int(match.group(1))
                if required_code is not None and model_code != required_code:
                    message = (
                        f"Prodigal-GV model uses genetic code {model_code}, "
                        f"expected {required_code}"
                    )
                    raise ValueError(message)
                has_model = True
    if not (has_gff3_header and has_sequence and has_model):
        raise ValueError(f"invalid zero-call Prodigal-GV GFF: {path}")


def _attributes(value: str, path: Path, line_number: int) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for field in value.rstrip(";").split(";"):
        if not field:
            continue
        key, separator, raw_value = field.partition("=")
        if not separator or not key or not raw_value:
            raise ValueError(f"invalid GFF attribute at {path}:{line_number}")
        if key in attributes:
            raise ValueError(f"duplicate GFF attribute {key!r} at {path}:{line_number}")
        attributes[key] = unquote(raw_value)
    return attributes


def _parse_prodigal_header(record: _FastaRecord) -> _ProdigalProtein:
    fields = [field.strip() for field in record.header.split(" # ")]
    if len(fields) != 5 or fields[0] != record.identifier:
        raise ValueError(f"invalid Prodigal-GV protein header: {record.header!r}")
    attributes = _attributes(fields[4], Path("<protein FASTA>"), 1)
    gene_id = attributes.get("ID")
    if not gene_id:
        raise ValueError(
            f"Prodigal-GV protein lacks native gene ID: {record.identifier}"
        )
    return _ProdigalProtein(
        record=record,
        gene_id=gene_id,
        start=_positive_integer(fields[1], "Prodigal-GV protein start"),
        end=_positive_integer(fields[2], "Prodigal-GV protein end"),
        strand={"1": "plus", "-1": "minus"}.get(fields[3])
        or _invalid_prodigal_strand(fields[3]),
        genetic_code=_positive_integer(
            attributes.get("genetic_code", ""), "Prodigal-GV protein genetic_code"
        ),
    )


def _invalid_prodigal_strand(value: str) -> str:
    raise ValueError(f"invalid Prodigal-GV protein strand: {value!r}")


def _parents(feature: _Feature, path: Path) -> tuple[str, ...]:
    value = _required_attribute(feature, "Parent", path)
    parents = tuple(parent for parent in value.split(",") if parent)
    if not parents:
        raise ValueError("empty GFF Parent attribute")
    return parents


def _required_attribute(feature: _Feature, key: str, path: Path) -> str:
    value = feature.attributes.get(key)
    if not value:
        raise ValueError(f"{feature.feature_type} lacks {key} attribute in {path}")
    return value


def _integer_attribute(feature: _Feature, key: str, path: Path) -> int:
    return _positive_integer(_required_attribute(feature, key, path), key)


def _add_unique(
    records: dict[str, _Feature],
    identifier: str,
    feature: _Feature,
    label: str,
) -> None:
    if identifier in records:
        raise ValueError(f"duplicate {label} ID: {identifier}")
    records[identifier] = feature


def _validate_contained(
    child: _Feature,
    parent: _Feature,
    child_id: str,
    parent_id: str,
) -> None:
    if (
        child.contig_id != parent.contig_id
        or child.strand != parent.strand
        or child.start < parent.start
        or child.end > parent.end
    ):
        raise ValueError(f"{child_id} is inconsistent with parent {parent_id}")


def _first_cds(cds_records: list[_Feature], strand: str) -> _Feature:
    if strand == "plus":
        return min(cds_records, key=lambda feature: (feature.start, feature.end))
    if strand == "minus":
        return max(cds_records, key=lambda feature: (feature.end, feature.start))
    raise ValueError("BRAKER3 transcript strand must be plus or minus")


def _gene_row(
    *,
    gene_id: str,
    feature: _Feature,
    bin_id: str | None,
    phase: int | None,
    genetic_code: int,
) -> GeneRow:
    return {
        "gene_id": gene_id,
        "contig_id": feature.contig_id,
        "bin_id": bin_id,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand,
        "phase": phase,
        "genetic_code": genetic_code,
    }


def _native_gene_id(caller: Caller, record: _FastaRecord) -> str:
    if caller == "braker3":
        match = BRAKER_TRANSCRIPT.fullmatch(record.identifier)
        if match is None:
            raise ValueError(
                f"BRAKER3 protein ID lacks .t<number> suffix: {record.identifier}"
            )
        return match.group(1)
    return _parse_prodigal_header(record).gene_id


def _normalized_id(
    caller: Caller,
    source_id: str,
    bin_id: str | None,
    native_id: str,
) -> str:
    _caller(caller)
    components = [caller, source_id]
    if bin_id is not None:
        components.append(bin_id)
    components.append(native_id)
    normalized = ":".join(_identifier_component(value) for value in components)
    _require_identifier(normalized, "normalized protein ID")
    return normalized


def _identifier_component(value: str) -> str:
    stripped = value.strip()
    if not stripped:
        raise ValueError("identifier component cannot be empty")
    if IDENTIFIER.fullmatch(stripped):
        return stripped
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", stripped).strip("_.-")
    if not cleaned or not cleaned[0].isalnum():
        cleaned = f"id_{cleaned}"
    digest = hashlib.sha256(stripped.encode()).hexdigest()[:10]
    return f"{cleaned}-{digest}"


def _caller(caller: str) -> Caller:
    if caller not in {"braker3", "prodigal-gv"}:
        raise ValueError(f"unsupported gene caller: {caller!r}")
    return cast("Caller", caller)


def _prodigal_protein_id(contig_id: str, gene_id: str) -> str:
    _sequence_number, separator, ordinal = gene_id.rpartition("_")
    if not separator or not ordinal.isdigit():
        raise ValueError(f"invalid Prodigal-GV CDS ID: {gene_id!r}")
    return f"{contig_id}_{ordinal}"


def _set_mismatch(
    left_label: str,
    left: dict[str, object],
    right: set[str],
) -> str:
    left_ids = set(left)
    missing = sorted(left_ids - right)
    extra = sorted(right - left_ids)
    return f"{left_label}/protein ID mismatch; missing={missing}; extra={extra}"


def _strand(value: str, path: Path, line_number: int) -> str:
    strands = {"+": "plus", "-": "minus", ".": "unknown", "?": "unknown"}
    if value not in strands:
        raise ValueError(f"invalid GFF strand at {path}:{line_number}: {value!r}")
    return strands[value]


def _phase(value: str, path: Path, line_number: int) -> int | None:
    if value == ".":
        return None
    try:
        phase = int(value)
    except ValueError as error:
        raise ValueError(
            f"invalid GFF phase at {path}:{line_number}: {value!r}"
        ) from error
    if phase not in {0, 1, 2}:
        raise ValueError(f"invalid GFF phase at {path}:{line_number}: {value!r}")
    return phase


def _positive_integer(value: str, field: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValueError(f"{field} must be an integer: {value!r}") from error
    if parsed < 1:
        raise ValueError(f"{field} must be at least 1: {value!r}")
    return parsed


def _optional_genetic_code(value: int | None) -> int | None:
    return None if value is None else _genetic_code(value)


def _genetic_code(value: int) -> int:
    if isinstance(value, bool) or value < 1:
        raise ValueError(f"genetic_code must be a positive integer: {value!r}")
    return value


def _require_identifier(value: str, field: str) -> None:
    if not IDENTIFIER.fullmatch(value):
        raise ValueError(f"invalid {field}: {value!r}")
