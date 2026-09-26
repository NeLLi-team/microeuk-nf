"""Parse native quality and taxonomy tables into catalog-ready record bodies."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

__all__ = ["ParsedRows", "parse_quality"]

MISSING_VALUES = {"", "na", "n/a", "nan", "none", "nd", "unavailable"}
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
SCORED_MARKERS = re.compile(r"markers\s+([0-9]+)/([0-9]+)")
CONFIDENCE_LEVEL = re.compile(r"^([^:]+):([0-9]+(?:\.[0-9]+)?)$")
FASTA_SUFFIXES = (".fasta", ".fna", ".faa", ".fa")
SYMCLATRON_THRESHOLD = 0.725
RANKS = {
    "d": "domain",
    "k": "kingdom",
    "sg": "supergroup",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

type NormalizedRow = dict[str, object]


@dataclass(frozen=True, slots=True)
class ParsedRows:
    """Normalized record bodies plus native evidence that lacks catalog fields."""

    qc: list[NormalizedRow]
    taxonomy: list[NormalizedRow]
    phenotypes: list[NormalizedRow]
    target_ids: tuple[str, ...]
    evidence_gaps: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class _ParsedToolRows:
    qc: list[NormalizedRow]
    taxonomy: list[NormalizedRow]
    phenotypes: list[NormalizedRow]
    target_ids: tuple[str, ...]
    evidence_gaps: tuple[str, ...] = ()


def parse_quality(
    tool: str,
    path: Path,
    detail_path: Path | None = None,
) -> ParsedRows:
    """Parse one native tool table without adding workflow envelope fields.

    Args:
        tool: Canonical tool name from the stage manifest.
        path: Native TSV output from that tool.
        detail_path: CheckEUK detail TSV with structured marker inventory fields.

    Raises:
        ValueError: If the tool, header, value, or target identifiers are invalid.
    """
    parser_name = _tool_name(tool)
    if detail_path is not None and parser_name != "checkeuk":
        raise ValueError("detail_path is supported only for CheckEUK")

    if parser_name == "checkeuk":
        parsed = _parse_checkeuk(path, detail_path)
    else:
        parser = {
            "checkm1": _parse_checkm1,
            "checkm2": _parse_checkm2,
            "gvclass": _parse_gvclass,
            "gtdbtk": _parse_gtdbtk,
            "symclatron": _parse_symclatron,
            "quickclade": _parse_quickclade,
        }[parser_name]
        parsed = parser(path)
    return ParsedRows(
        qc=parsed.qc,
        taxonomy=parsed.taxonomy,
        phenotypes=parsed.phenotypes,
        target_ids=parsed.target_ids,
        evidence_gaps=parsed.evidence_gaps + _evidence_gaps(parser_name, detail_path),
    )


def _tool_name(tool: str) -> str:
    names = {
        "checkeuk": "checkeuk",
        "checkm1": "checkm1",
        "checkm2": "checkm2",
        "gvclass": "gvclass",
        "gtdb-tk": "gtdbtk",
        "gtdbtk": "gtdbtk",
        "symclatron": "symclatron",
        "quickclade": "quickclade",
    }
    name = tool.strip().casefold()
    if name not in names:
        supported = ", ".join(sorted(names))
        message = f"unsupported quality tool {tool!r}; expected one of {supported}"
        raise ValueError(message)
    return names[name]


def _parse_checkeuk(path: Path, detail_path: Path | None) -> _ParsedToolRows:
    required = {
        "genome",
        "status",
        "completeness",
        "completeness_basis",
        "completeness_bom",
        "contamination",
        "lineage",
        "completeness_note",
        "pfam_status",
        "pfam_union_families",
        "pfam_observed_families",
        "pfam_completeness",
    }
    rows, _ = _read_tsv(path, required)
    detail = _checkeuk_detail(detail_path) if detail_path is not None else {}
    qc_rows: list[NormalizedRow] = []
    taxonomy_rows: list[NormalizedRow] = []
    evidence_gaps: list[str] = []
    observed: set[str] = set()
    for row in rows:
        target = _unique_target(row["genome"], observed, "CheckEUK")
        status = row["status"].strip()
        if status == "filtered":
            reason = row["completeness_note"].strip() or "reason not reported"
            evidence_gaps.append(
                f"CheckEUK filtered {target} before scoring: {reason}."
            )
            continue
        if status != "ok":
            raise ValueError(
                f"CheckEUK reported unsupported status for {target}: {status}"
            )
        qc_rows.append(_checkeuk_qc(row, detail.get(target), target))
        taxonomy_rows.extend(_taxonomy("checkeuk", target, row["lineage"]))
    if detail_path is not None and set(detail) != observed:
        raise ValueError("CheckEUK summary and detail target sets differ")
    return _ParsedToolRows(
        qc_rows,
        taxonomy_rows,
        [],
        tuple(sorted(observed)),
        tuple(evidence_gaps),
    )


def _checkeuk_qc(
    row: dict[str, str],
    detail: dict[str, str] | None,
    target: str,
) -> NormalizedRow:
    completeness = _float(row["completeness"], "completeness", 0, 100)
    contamination = _float(row["contamination"], "contamination", 0, 100)
    pfam_completeness, pfam_observed = _checkeuk_pfam(row, target)
    markers, inventory = _checkeuk_markers(row, detail, target)
    qc = _qc_row("checkeuk", target)
    _set(qc, "completeness_percent", completeness)
    if completeness is not None:
        qc["completeness_basis"] = _required_text(
            row["completeness_basis"], "completeness_basis"
        )
    _set(qc, "contamination_percent", contamination)
    if contamination is not None:
        qc["contamination_basis"] = "CheckEUK report contamination"
    _set(qc, "checkeuk_pfam_ridge_completeness_percent", pfam_completeness)
    _set(
        qc,
        "checkeuk_bom_completeness_percent",
        _float(row["completeness_bom"], "completeness_bom", 0, 100),
    )
    _set(qc, "checkeuk_inventory74_present", inventory)
    _set(qc, "checkeuk_expected_markers_present", markers[0])
    _set(qc, "checkeuk_expected_markers_total", markers[1])
    _set(qc, "checkeuk_pfam1085_distinct", pfam_observed)
    if len(qc) == 3:
        raise ValueError(f"CheckEUK row contains no quality measurements for {target}")
    return qc


def _checkeuk_pfam(
    row: dict[str, str],
    target: str,
) -> tuple[float | None, int | None]:
    completeness = _float(row["pfam_completeness"], "pfam_completeness", 0, 100)
    union = _integer(row["pfam_union_families"], "pfam_union_families")
    observed = _integer(row["pfam_observed_families"], "pfam_observed_families")
    if observed is not None and union != 1085:
        raise ValueError(f"CheckEUK Pfam panel is not 1,085 families for {target}")
    if observed is not None and observed > 1085:
        raise ValueError(f"CheckEUK Pfam recovery exceeds 1,085 for {target}")
    return completeness, observed


def _checkeuk_markers(
    row: dict[str, str],
    detail: dict[str, str] | None,
    target: str,
) -> tuple[tuple[int | None, int | None], int | None]:
    markers = _scored_markers(row["completeness_note"])
    if detail is None:
        return markers, None
    structured = (
        _required_integer(detail["n_present_expected"], "n_present_expected"),
        _required_integer(detail["n_expected"], "n_expected"),
    )
    _validate_marker_fraction(*structured, target)
    if markers not in ((None, None), structured):
        raise ValueError(f"CheckEUK summary/detail marker counts differ for {target}")
    inventory = _required_integer(detail["n_present_core"], "n_present_core")
    if not 0 <= inventory <= 74:
        raise ValueError(f"CheckEUK 74-marker inventory is invalid for {target}")
    return structured, inventory


def _checkeuk_detail(path: Path) -> dict[str, dict[str, str]]:
    rows, _ = _read_tsv(
        path,
        {"genome", "n_present_core", "n_present_expected", "n_expected"},
    )
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        target = _target_id(row["genome"])
        if target in indexed:
            raise ValueError(f"duplicate CheckEUK detail target: {target}")
        indexed[target] = row
    return indexed


def _parse_checkm1(path: Path) -> _ParsedToolRows:
    required = {"Bin Id", "Marker lineage", "Completeness", "Contamination"}
    rows, _ = _read_tsv(path, required)
    qc_rows: list[NormalizedRow] = []
    observed: set[str] = set()
    for row in rows:
        target = _unique_target(row["Bin Id"], observed, "CheckM1")
        completeness = _required_float(row["Completeness"], "Completeness", 0, 100)
        contamination = _required_float(row["Contamination"], "Contamination", 0, 100)
        basis = row["Marker lineage"].strip()
        if basis == "":
            raise ValueError(f"CheckM1 marker lineage is empty for {target}")
        qc_rows.append(
            {
                **_qc_row("checkm1", target),
                "completeness_percent": completeness,
                "contamination_percent": contamination,
                "completeness_basis": basis,
                "contamination_basis": "CheckM1 marker-set multi-copy estimate",
            }
        )
    return _ParsedToolRows(qc_rows, [], [], tuple(sorted(observed)))


def _parse_checkm2(path: Path) -> _ParsedToolRows:
    required = {
        "Name",
        "Completeness",
        "Contamination",
        "Completeness_Model_Used",
    }
    rows, _ = _read_tsv(path, required)
    qc_rows: list[NormalizedRow] = []
    observed: set[str] = set()
    for row in rows:
        target = _unique_target(row["Name"], observed, "CheckM2")
        basis = row["Completeness_Model_Used"].strip()
        if basis == "":
            raise ValueError(f"CheckM2 completeness model is empty for {target}")
        qc_rows.append(
            {
                **_qc_row("checkm2", target),
                "completeness_percent": _required_float(
                    row["Completeness"], "Completeness", 0, 100
                ),
                "contamination_percent": _required_float(
                    row["Contamination"], "Contamination", 0, None
                ),
                "completeness_basis": basis,
                "contamination_basis": "CheckM2 reported contamination",
            }
        )
    return _ParsedToolRows(qc_rows, [], [], tuple(sorted(observed)))


def _parse_gvclass(path: Path) -> _ParsedToolRows:
    required = {
        "query",
        "taxonomy_majority",
        "taxonomy_confidence",
        "estimated_completeness",
        "estimated_completeness_strategy",
        "completeness_query_support",
        "completeness_model_group",
        "completeness_model_reliability",
        "estimated_contamination",
        "estimated_contamination_strategy",
        "contamination_query_support",
        "contamination_model_reliability",
        "contamination_type",
    }
    rows, _ = _read_tsv(path, required)
    qc_rows: list[NormalizedRow] = []
    taxonomy_rows: list[NormalizedRow] = []
    observed: set[str] = set()
    for row in rows:
        target = _unique_target(row["query"], observed, "GVClass")
        completeness = _float(
            row["estimated_completeness"], "estimated_completeness", 0, 100
        )
        contamination = _float(
            row["estimated_contamination"], "estimated_contamination", 0, 100
        )
        if completeness is not None or contamination is not None:
            qc = _qc_row("gvclass", target)
            _set(qc, "completeness_percent", completeness)
            _set(qc, "contamination_percent", contamination)
            if completeness is not None:
                qc["completeness_basis"] = _required_text(
                    row["estimated_completeness_strategy"],
                    "estimated_completeness_strategy",
                )
            if contamination is not None:
                qc["contamination_basis"] = _required_text(
                    row["estimated_contamination_strategy"],
                    "estimated_contamination_strategy",
                )
            for field in (
                "completeness_query_support",
                "completeness_model_reliability",
                "contamination_query_support",
                "contamination_model_reliability",
            ):
                qc[field] = _required_text(row[field], field, allow_unavailable=True)
            qc_rows.append(qc)
        taxonomy = _taxonomy("gvclass", target, row["taxonomy_majority"])
        if taxonomy:
            lineage_confidence = _required_text(
                row["taxonomy_confidence"], "taxonomy_confidence"
            )
            for taxonomy_row in taxonomy:
                taxonomy_row["lineage_confidence"] = lineage_confidence
        taxonomy_rows.extend(taxonomy)
    return _ParsedToolRows(qc_rows, taxonomy_rows, [], tuple(sorted(observed)))


def _parse_gtdbtk(path: Path) -> _ParsedToolRows:
    required = {"user_genome", "classification", "classification_method"}
    rows, _ = _read_tsv(path, required)
    taxonomy_rows: list[NormalizedRow] = []
    observed: set[str] = set()
    for row in rows:
        target = _unique_target(row["user_genome"], observed, "GTDB-Tk")
        classification = _required_text(row["classification"], "classification")
        parsed = _taxonomy("gtdbtk", target, classification)
        if not parsed:
            parsed = [_taxonomy_row("gtdbtk", target, "domain", classification, None)]
        taxonomy_rows.extend(parsed)
    return _ParsedToolRows([], taxonomy_rows, [], tuple(sorted(observed)))


def _parse_symclatron(path: Path) -> _ParsedToolRows:
    required = {
        "taxon_oid",
        "completeness_UNI56",
        "classification",
        "confidence",
        "passes_confidence_threshold",
        "classification_thresholded",
    }
    rows, _ = _read_tsv(path, required)
    qc_rows: list[NormalizedRow] = []
    phenotype_rows: list[NormalizedRow] = []
    observed: set[str] = set()
    for row in rows:
        target = _unique_target(row["taxon_oid"], observed, "Symclatron")
        confidence = _required_float(row["confidence"], "confidence", 0, 1)
        passes = _boolean(row["passes_confidence_threshold"])
        if passes != (confidence >= SYMCLATRON_THRESHOLD):
            raise ValueError(f"Symclatron threshold flag disagrees for {target}")
        qc_rows.append(
            {
                **_qc_row("symclatron", target),
                "completeness_percent": _required_float(
                    row["completeness_UNI56"], "completeness_UNI56", 0, 100
                ),
                "completeness_basis": "UNI56 marker inventory (56 markers)",
            }
        )
        phenotype_rows.append(
            {
                "record_id": f"symclatron:phenotype:{target}",
                "bin_id": target,
                "raw_class": _required_text(row["classification"], "classification"),
                "raw_score": confidence,
                "threshold": SYMCLATRON_THRESHOLD,
                "passes_threshold": passes,
                "thresholded_class": _required_text(
                    row["classification_thresholded"],
                    "classification_thresholded",
                ),
                "is_applicable": False,
                "applicability_reason": (
                    "Requires supported prokaryotic routing and completeness policy"
                ),
            }
        )
    return _ParsedToolRows(qc_rows, [], phenotype_rows, tuple(sorted(observed)))


def _parse_quickclade(path: Path) -> _ParsedToolRows:
    required = {"#QueryName", "lineage", "ConfLevel", "Confidence"}
    rows, _ = _read_tsv(path, required)
    taxonomy_rows: list[NormalizedRow] = []
    observed: set[str] = set()
    for row in rows:
        target = _unique_target(row["#QueryName"], observed, "QuickClade")
        score_by_rank = _quickclade_scores(row["Confidence"], target)
        confidence = row["ConfLevel"].strip()
        if confidence != "":
            match = CONFIDENCE_LEVEL.fullmatch(confidence)
            if match is None:
                raise ValueError(f"invalid QuickClade ConfLevel for {target}")
            rank = _rank(match.group(1))
            score = _bounded_score(match.group(2), "ConfLevel")
            if rank in score_by_rank and score_by_rank[rank] != score:
                raise ValueError(f"QuickClade confidence fields disagree for {target}")
            score_by_rank[rank] = score
        taxonomy_rows.extend(
            _taxonomy("quickclade", target, row["lineage"], score_by_rank)
        )
    return _ParsedToolRows([], taxonomy_rows, [], tuple(sorted(observed)))


def _read_tsv(
    path: Path,
    required: set[str],
) -> tuple[list[dict[str, str]], tuple[str, ...]]:
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
    return rows, tuple(fields)


def _target_id(value: str) -> str:
    name = Path(value.strip()).name
    lowered = name.casefold()
    for suffix in FASTA_SUFFIXES:
        if lowered.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if not IDENTIFIER.fullmatch(name):
        raise ValueError(f"invalid native target identifier: {value!r}")
    return name


def _unique_target(value: str, observed: set[str], tool: str) -> str:
    target = _target_id(value)
    if target in observed:
        raise ValueError(f"duplicate {tool} target: {target}")
    observed.add(target)
    return target


def _float(
    value: str,
    field: str,
    minimum: float | None,
    maximum: float | None,
) -> float | None:
    text = value.strip()
    if text.casefold() in MISSING_VALUES:
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
    minimum: float | None,
    maximum: float | None,
) -> float:
    number = _float(value, field, minimum, maximum)
    if number is None:
        raise ValueError(f"{field} is missing")
    return number


def _integer(value: str, field: str) -> int | None:
    text = value.strip()
    if text.casefold() in MISSING_VALUES:
        return None
    if not text.isdigit():
        raise ValueError(f"{field} is not an unsigned integer")
    return int(text)


def _required_integer(value: str, field: str) -> int:
    number = _integer(value, field)
    if number is None:
        raise ValueError(f"{field} is missing")
    return number


def _required_text(
    value: str,
    field: str,
    *,
    allow_unavailable: bool = False,
) -> str:
    text = value.strip()
    normalized = text.casefold()
    if normalized in MISSING_VALUES and not (
        allow_unavailable and normalized == "unavailable"
    ):
        raise ValueError(f"{field} is missing")
    return text


def _boolean(value: str) -> bool:
    normalized = value.strip().casefold()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ValueError("boolean field must be True or False")


def _scored_markers(note: str) -> tuple[int | None, int | None]:
    match = SCORED_MARKERS.search(note)
    if match is None:
        return None, None
    present, total = int(match.group(1)), int(match.group(2))
    _validate_marker_fraction(present, total)
    return present, total


def _validate_marker_fraction(present: int, total: int, target: str = "") -> None:
    if 0 <= present <= total <= 74 and total > 0:
        return
    suffix = f" for {target}" if target else ""
    raise ValueError(f"CheckEUK expected-marker counts are invalid{suffix}")


def _qc_row(tool: str, target: str) -> NormalizedRow:
    return {
        "record_id": f"{tool}:qc:{target}",
        "target_kind": "bin",
        "target_id": target,
    }


def _set(row: NormalizedRow, field: str, value: object | None) -> None:
    if value is not None:
        row[field] = value


def _taxonomy(
    tool: str,
    target: str,
    lineage: str,
    score_by_rank: dict[str, float] | None = None,
) -> list[NormalizedRow]:
    text = lineage.strip()
    if text.casefold() in MISSING_VALUES:
        return []
    scores = score_by_rank or {}
    records: list[NormalizedRow] = []
    observed_ranks: set[str] = set()
    for token in text.split(";"):
        match = re.fullmatch(r"([a-z]+)_{1,2}(.*)", token.strip())
        if match is None or match.group(2).strip() == "":
            continue
        rank = RANKS.get(match.group(1).casefold(), match.group(1).casefold())
        if rank in observed_ranks:
            raise ValueError(f"duplicate taxonomy rank {rank} for {target}")
        observed_ranks.add(rank)
        records.append(
            _taxonomy_row(
                tool,
                target,
                rank,
                match.group(2).strip(),
                scores.get(rank),
            )
        )
    return records


def _quickclade_scores(value: str, target: str) -> dict[str, float]:
    text = value.strip()
    if text.casefold() in MISSING_VALUES:
        return {}
    scores: dict[str, float] = {}
    for field in text.split(";"):
        key, separator, raw_score = field.partition(":")
        if separator == "" or key.strip() == "":
            raise ValueError(f"invalid QuickClade Confidence for {target}")
        rank = _rank(key)
        if rank in scores:
            message = f"duplicate QuickClade confidence rank {rank} for {target}"
            raise ValueError(message)
        scores[rank] = _bounded_score(raw_score, "Confidence")
    return scores


def _rank(value: str) -> str:
    normalized = value.strip().casefold()
    return RANKS.get(normalized, normalized)


def _bounded_score(value: str, field: str) -> float:
    return _required_float(value, field, 0, 100)


def _taxonomy_row(
    tool: str,
    target: str,
    rank: str,
    taxon: str,
    score: float | None,
) -> NormalizedRow:
    row: NormalizedRow = {
        "record_id": f"{tool}:taxonomy:{target}:{rank}",
        "target_kind": "bin",
        "target_id": target,
        "rank": rank,
        "taxon_name": taxon,
    }
    _set(row, "score", score)
    return row


def _evidence_gaps(tool: str, detail_path: Path | None) -> tuple[str, ...]:
    gaps = {
        "checkeuk": (
            (
                "CheckEUK status, warnings, contamination type, and Pfam diagnostics "
                "have no normalized catalog fields."
            ),
        ),
        "checkm1": (
            (
                "CheckM1 marker counts, copy histogram, and strain heterogeneity "
                "have no normalized catalog fields."
            ),
        ),
        "checkm2": (
            (
                "CheckM2 translation table, coding density, model diagnostics, and "
                "notes have no normalized catalog fields."
            ),
        ),
        "gvclass": (
            (
                "GVClass confidence, query-support, and model-reliability labels are "
                "normalized only on existing result rows. Metadata for targets without "
                "those rows remains in the retained native report; model group, "
                "contamination type, species-tree evidence, and marker panels remain "
                "unnormalized."
            ),
        ),
        "gtdbtk": (
            (
                "GTDB-Tk classification method, ANI, aligned fraction, RED, MSA "
                "coverage, and warnings have no normalized catalog fields."
            ),
        ),
        "symclatron": (
            (
                "Symclatron applicability requires downstream taxonomy and "
                "completeness routing; the parser marks output inapplicable."
            ),
        ),
        "quickclade": (
            (
                "QuickClade reference identity and GC, STR, HH, CAGA, and k-mer "
                "distances have no normalized catalog fields."
            ),
        ),
    }
    result = list(gaps[tool])
    if tool == "checkeuk" and detail_path is None:
        result.append(
            "CheckEUK 74-marker inventory is unavailable without "
            "checkeuk_report_detail.tsv."
        )
    return tuple(result)
