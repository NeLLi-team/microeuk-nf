"""Normalize native protein annotation tables into catalog row bodies."""

import hashlib
import math
import re
from pathlib import Path

__all__ = ["parse_eggnog", "parse_interpro"]

type Annotation = dict[str, object]
type Reference = tuple[str, str, str]

EGGNOG_HEADER = (
    "#query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
)
EGGNOG_TERM_FIELDS = (
    ("GOs", "Gene Ontology via eggNOG", re.compile(r"GO:\d{7}")),
    (
        "EC",
        "Enzyme Commission via eggNOG",
        re.compile(r"\d+(?:\.(?:\d+|-)){3}"),
    ),
    ("KEGG_ko", "KEGG Orthology via eggNOG", re.compile(r"ko:K\d{5}")),
)
OG_PATTERN = re.compile(r"[^@|\s]+@[0-9]+")
INTERPRO_PATTERN = re.compile(r"IPR\d{6}")
MISSING = {"", "-"}


def parse_eggnog(
    path: Path,
    protein_ids: set[str],
    database_version: str,
) -> list[Annotation]:
    """Parse an eggNOG-mapper 2.1.15 annotations table.

    Protein-level eggNOG records omit ``start`` and ``end`` because the native
    annotations table does not report coordinates. Orthologous-group accessions
    retain ``OG@taxid``; the native artifact retains the taxonomic display label.

    Args:
        path: Native ``*.emapper.annotations`` file.
        protein_ids: Protein FASTA identifiers accepted for this stage.
        database_version: Combined eggNOG database release label.

    Returns:
        Deterministically ordered functional annotation row bodies.

    Raises:
        ValueError: If the native table is malformed or names an unknown protein.
    """
    _validate_version(database_version)
    rows = _eggnog_rows(path)
    annotations: list[Annotation] = []
    seen_queries: set[str] = set()
    for row in rows:
        gene_id = row["#query"]
        _validate_gene(gene_id, protein_ids, path)
        if gene_id in seen_queries:
            raise ValueError(f"eggNOG repeats query {gene_id}: {path}")
        seen_queries.add(gene_id)
        score = _float(row["score"], "eggNOG score", gene_id)
        evalue = _float(row["evalue"], "eggNOG evalue", gene_id)
        if evalue is not None and evalue < 0:
            raise ValueError(f"eggNOG evalue is negative for {gene_id}")
        evidence = _numeric_evidence(score, evalue)
        description = _optional(row["Description"])
        for accession in _orthologous_groups(row["eggNOG_OGs"]):
            annotation = _annotation(
                ("eggnog", "eggNOG", database_version),
                gene_id,
                accession,
                function_name=description,
            )
            annotation.update(evidence)
            annotations.append(annotation)
        for field, database_name, pattern in EGGNOG_TERM_FIELDS:
            for accession in _terms(row[field]):
                if pattern.fullmatch(accession) is None:
                    raise ValueError(f"invalid {field} accession: {accession}")
                annotation = _annotation(
                    ("eggnog", database_name, database_version),
                    gene_id,
                    accession,
                )
                annotation.update(evidence)
                annotations.append(annotation)
    return _deduplicate(annotations)


def parse_interpro(
    path: Path,
    protein_ids: set[str],
    database_version: str,
) -> list[Annotation]:
    """Parse InterProScan 5.76 TSV member and integrated annotations.

    Native scores belong to member analyses. Integrated InterPro records do
    not inherit those scores.

    Args:
        path: Native InterProScan TSV output.
        protein_ids: Protein FASTA identifiers accepted for this stage.
        database_version: Combined InterProScan and data release label.

    Returns:
        Deterministically ordered functional annotation row bodies.

    Raises:
        ValueError: If a row is malformed or names an unknown protein.
    """
    _validate_version(database_version)
    annotations: list[Annotation] = []
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) not in {11, 13, 14, 15}:
                raise ValueError(
                    f"InterProScan row {line_number} has {len(fields)} columns: {path}"
                )
            annotations.extend(
                _interpro_annotations(
                    fields, protein_ids, database_version, path, lengths
                )
            )
    return _deduplicate(annotations)


def _interpro_annotations(
    fields: list[str],
    protein_ids: set[str],
    database_version: str,
    path: Path,
    lengths: dict[str, int],
) -> list[Annotation]:
    gene_id = fields[0]
    _validate_gene(gene_id, protein_ids, path)
    length = _positive_integer(fields[2], "protein length", gene_id)
    previous_length = lengths.setdefault(gene_id, length)
    if previous_length != length:
        raise ValueError(f"InterProScan protein length changes for {gene_id}")
    analysis = _required(fields[3], "analysis", gene_id)
    accession = _required(fields[4], "signature accession", gene_id)
    start = _positive_integer(fields[6], "start", gene_id)
    end = _positive_integer(fields[7], "end", gene_id)
    if end < start or end > length:
        raise ValueError(f"invalid InterProScan coordinates for {gene_id}")
    if fields[9] != "T":
        raise ValueError(f"unsupported InterProScan status for {gene_id}: {fields[9]}")
    score = _float(fields[8], "InterProScan score", gene_id)
    member = _annotation(
        ("interproscan", f"InterProScan:{analysis}", database_version),
        gene_id,
        accession,
        function_name=_optional(fields[5]),
        coordinates=(start, end),
    )
    if score is not None:
        member["score"] = score
    annotations = [member]
    if len(fields) < 13 or fields[11] in MISSING:
        return annotations
    integrated_accession = fields[11]
    if INTERPRO_PATTERN.fullmatch(integrated_accession) is None:
        raise ValueError(
            f"invalid InterPro accession for {gene_id}: {integrated_accession}"
        )
    integrated = _annotation(
        ("interproscan", "InterPro", database_version),
        gene_id,
        integrated_accession,
        function_name=_optional(fields[12]),
        coordinates=(start, end),
    )
    annotations.append(integrated)
    return annotations


def _eggnog_rows(path: Path) -> list[dict[str, str]]:
    header: tuple[str, ...] | None = None
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("##") or not raw.strip():
                continue
            values = tuple(raw.rstrip("\r\n").split("\t"))
            if header is None:
                if values not in {EGGNOG_HEADER, (*EGGNOG_HEADER, "md5")}:
                    raise ValueError(f"unexpected eggNOG header: {path}")
                header = values
                continue
            if raw.startswith("#"):
                raise ValueError(
                    f"unexpected eggNOG comment at row {line_number}: {path}"
                )
            if len(values) != len(header):
                raise ValueError(f"malformed eggNOG row {line_number}: {path}")
            rows.append(dict(zip(header, values, strict=True)))
    if header is None:
        raise ValueError(f"eggNOG annotations file has no header: {path}")
    return rows


def _annotation(
    reference: Reference,
    gene_id: str,
    accession: str,
    *,
    function_name: str | None = None,
    coordinates: tuple[int, int] | None = None,
) -> Annotation:
    source, database_name, database_version = reference
    start, end = coordinates or (None, None)
    identity = "\x1f".join(
        (source, gene_id, database_name, accession, str(start), str(end))
    )
    digest = hashlib.sha256(identity.encode()).hexdigest()
    row: Annotation = {
        "annotation_id": f"ann:{source}:{digest}",
        "gene_id": gene_id,
        "database_name": database_name,
        "database_version": database_version,
        "accession": accession,
    }
    if function_name is not None:
        row["function_name"] = function_name
    if start is not None and end is not None:
        row.update({"start": start, "end": end})
    return row


def _deduplicate(rows: list[Annotation]) -> list[Annotation]:
    unique: dict[str, Annotation] = {}
    for row in rows:
        identifier = str(row["annotation_id"])
        previous = unique.setdefault(identifier, row)
        if previous != row:
            raise ValueError(f"conflicting annotation rows share {identifier}")
    return sorted(
        unique.values(),
        key=lambda row: (
            str(row["gene_id"]),
            str(row["database_name"]),
            str(row["accession"]),
            int(row.get("start", 0)),
            int(row.get("end", 0)),
        ),
    )


def _orthologous_groups(value: str) -> tuple[str, ...]:
    if value in MISSING:
        return ()
    accessions: set[str] = set()
    for term in value.split(","):
        accession, separator, label = term.strip().partition("|")
        if OG_PATTERN.fullmatch(accession) is None or (
            separator and (not label.strip() or "|" in label)
        ):
            raise ValueError(f"invalid eggNOG orthologous group: {term}")
        accessions.add(accession)
    return tuple(sorted(accessions))


def _terms(value: str) -> tuple[str, ...]:
    if value in MISSING:
        return ()
    terms = tuple(part.strip() for part in value.split(","))
    if any(
        not term or any(character.isspace() for character in term) for term in terms
    ):
        raise ValueError(f"malformed functional accession list: {value}")
    return tuple(sorted(set(terms)))


def _float(value: str, field: str, gene_id: str) -> float | None:
    if value in MISSING:
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"invalid {field} for {gene_id}: {value}") from error
    if not math.isfinite(number):
        raise ValueError(f"non-finite {field} for {gene_id}: {value}")
    return number


def _numeric_evidence(score: float | None, evalue: float | None) -> Annotation:
    evidence: Annotation = {}
    if score is not None:
        evidence["score"] = score
    if evalue is not None:
        evidence["evalue"] = evalue
    return evidence


def _positive_integer(value: str, field: str, gene_id: str) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise ValueError(f"invalid {field} for {gene_id}: {value}") from error
    if number < 1:
        raise ValueError(f"invalid {field} for {gene_id}: {value}")
    return number


def _optional(value: str) -> str | None:
    text = value.strip()
    return None if text in MISSING else text


def _required(value: str, field: str, gene_id: str) -> str:
    if value in MISSING or any(character.isspace() for character in value):
        raise ValueError(f"invalid {field} for {gene_id}: {value}")
    return value


def _validate_gene(gene_id: str, protein_ids: set[str], path: Path) -> None:
    if gene_id not in protein_ids:
        raise ValueError(f"annotation names unknown protein {gene_id}: {path}")


def _validate_version(database_version: str) -> None:
    if not database_version.strip() or database_version != database_version.strip():
        raise ValueError("database_version must be nonempty without outer whitespace")
