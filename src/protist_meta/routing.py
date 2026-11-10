"""Route bins to one gene-calling path from independent native evidence."""

import csv
import json
import math
import shutil
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import Literal, NamedTuple

from protist_meta.sequences import read_fasta

__all__ = ["route_bins"]

type Row = dict[str, str]
type Route = Literal["prokaryotic", "eukaryotic", "viral", "unresolved"]
type CandidateClass = Literal[
    "eukaryotic_candidate",
    "bacterial_archaeal_candidate",
    "viral_candidate",
    "conflicting",
    "unresolved",
]

REPORTS = {
    "quickclade": Path("quickclade.tsv"),
    "checkeuk": Path("result_checkeuk/checkeuk_report.tsv"),
    "gvclass": Path("output/gvclass_summary.tsv"),
    "ssuextract": Path("output/cmsearch_summary.tsv"),
    "viral": Path("genomad/input_summary/input_virus_summary.tsv"),
    "checkm2": Path("output/quality_report.tsv"),
    "tiara": Path("tiara.tsv"),
}
ROUTES: tuple[Route, ...] = (
    "prokaryotic",
    "eukaryotic",
    "viral",
    "unresolved",
)
QUICKCLADE_DOMAINS = {
    "d__Bacteria": "prokaryotic",
    "d__Archaea": "prokaryotic",
    "d__Eukaryota": "eukaryotic",
    "sk__Bacteria": "prokaryotic",
    "sk__Archaea": "prokaryotic",
    "sk__Eukaryota": "eukaryotic",
}
GVCLASS_DOMAINS = {
    "d_BAC": "prokaryotic",
    "d_ARC": "prokaryotic",
    "d_EUK": "eukaryotic",
    "d_NCLDV": "viral",
    "d_PPV": "viral",
    "d_MIRUS": "viral",
    "d_PHAGE": "viral",
    "d_PLASTID": "organelle",
    "d_MITO": "organelle",
    "d_EUK-pEVE": "ambiguous",
}
UNCLASSIFIED_GVCLASS = {"", "d_", "nd"}
TIARA_FIRST_STAGE = {
    "archaea": "prokaryote",
    "bacteria": "prokaryote",
    "prokarya": "prokaryote",
    "eukarya": "eukaryote",
    "organelle": "organelle",
    "unknown": "unclassified",
}
TIARA_ORGANELLE_SECOND_STAGE = frozenset({"mitochondrion", "plastid"})
TIARA_CLASSES = ("eukaryote", "prokaryote", "organelle", "unclassified")
CHECKEUK_DETAIL = Path("result_checkeuk/checkeuk_report_detail.tsv")
CHECKEUK_VOTING_MARKERS = "marker_tax_host_lineage_selected_markers"
MIN_CHECKEUK_VOTING_MARKERS = 5
EUKARYOTIC_SSU_MODEL = "RF01960"
EUKARYOTIC_SSU_DOMAIN = "Eukaryota"
PROVIRUS_SEPARATOR = "|provirus_"
DUPLICATED_CONTIG_MARKER = "_duplicated-yes"
# A Tiara class or whole-contig geNomad calls speak for the bin from this share of
# its base pairs; smaller shares are recorded but never decide or veto a route.
MIN_TIARA_FRACTION = 0.5
MIN_VIRAL_CONTIG_FRACTION = 0.5
# Tiara prokaryote shares are bimodal in bins that carry eukaryote labels: below 0.2
# for real eukaryotes, at least 0.8 for bacteria. A partial share stays unscored.
TIARA_STRONG_PROKARYOTE = 0.8
# Short AT-rich ciliate nanochromosomes look like organelle DNA to Tiara, while the
# organelle bins seen in real samples are smaller than 0.6 Mb.
MAX_TIARA_ORGANELLE_BIN_BP = 2_000_000
MIN_CHECKM2_COMPLETENESS = 10.0

SSU_LINE = "a physically linked eukaryotic 18S locus"
TIARA_LINE = "Tiara"
FALLBACK_REASON = (
    "GVClass supports a eukaryotic route; CheckEUK has too few markers to name "
    "a lineage"
)
CONFLICT_REASON = "supported evidence conflicts across routes"
AMBIGUOUS_REASON = "GVClass reports an ambiguous eukaryote-virus class"
GVCLASS_ORGANELLE_REASON = (
    "GVClass reports an organelle class without a gene-calling route"
)
TIARA_ORGANELLE_REASON = (
    "Tiara classifies most of the bin as organelle sequence without a "
    "gene-calling route"
)
TIARA_CHECKM2_REASON = "Tiara and CheckM2 support a prokaryotic route"
GVCLASS_PROKARYOTIC_REASON = (
    "GVClass supports a prokaryotic route without a "
    "supported eukaryotic or viral conflict"
)
VIRAL_REASON = "GVClass or linked whole-contig geNomad evidence supports a viral route"

EVIDENCE_FIELDS = (
    "sample_id",
    "run_id",
    "bin_id",
    "route",
    "candidate_class",
    "reason",
    "checkeuk_status",
    "checkeuk_lineage",
    "taxonomy",
    "taxonomy_source",
    "quickclade_domain",
    "quickclade_lineage",
    "quickclade_confidence",
    "gvclass_domain",
    "gvclass_taxonomy_majority",
    "gvclass_taxonomy_confidence",
    "linked_eukaryotic_ssu",
    "linked_genomad_viruses",
    "supporting_evidence",
    "conflicting_evidence",
    "bin_bp",
    "tiara_eukaryote_fraction",
    "tiara_prokaryote_fraction",
    "tiara_organelle_fraction",
    "tiara_unclassified_fraction",
    "viral_contig_fraction",
    "checkm2_completeness",
    "checkm2_contamination",
    "duplicated_contig_fraction",
)


@dataclass(frozen=True, slots=True)
class BinComposition:
    """Base-pair shares of one bin that contig-level native calls account for."""

    bin_bp: int
    tiara_eukaryote: float
    tiara_prokaryote: float
    tiara_organelle: float
    tiara_unclassified: float
    viral_contig: float
    duplicated_contig: float

    @property
    def tiara_eukaryote_derived(self) -> float:
        """Share of eukaryote and organelle sequence, which is eukaryote-derived."""
        return self.tiara_eukaryote + self.tiara_organelle


@dataclass(frozen=True, slots=True)
class Checkm2Call:
    """CheckM2 completeness and contamination in percent."""

    completeness: float
    contamination: float


@dataclass(frozen=True, slots=True)
class BinEvidence:
    """Native evidence associated with one physical bin."""

    checkeuk: Row
    checkeuk_voting_markers: int
    quickclade: Row
    gvclass: Row
    eukaryotic_ssu: tuple[Row, ...]
    genomad_viruses: tuple[Row, ...]
    composition: BinComposition
    checkm2: Checkm2Call | None


@dataclass(frozen=True, slots=True)
class NativeEvidence:
    """Every native table, read once and indexed by bin or contig."""

    quickclade: Mapping[str, Row]
    checkeuk: Mapping[str, Row]
    gvclass: Mapping[str, Row]
    voting_markers: Mapping[str, int]
    linked_ssu: Mapping[str, list[Row]]
    linked_viruses: Mapping[str, list[Row]]
    tiara: Mapping[str, str]
    checkm2: Mapping[str, Checkm2Call]


@dataclass(frozen=True, slots=True)
class Decision:
    """One exclusive route and its catalog interpretation."""

    route: Route
    candidate_class: CandidateClass
    reason: str
    supporting: tuple[str, ...]
    conflicting: tuple[str, ...]


class Outcome(NamedTuple):
    """The route, class and reason chosen by one routing rule."""

    route: Route
    candidate_class: CandidateClass
    reason: str


@dataclass(frozen=True, slots=True)
class Signals:
    """Normalized support flags used by the routing policy.

    Attributes:
        quickclade: QuickClade route label; it adds evidence items, never a route.
        gvclass: GVClass route label for its majority domain.
        checkeuk_resolved: CheckEUK reported a lineage.
        checkeuk_fallback: GVClass stands in for a CheckEUK run with too few markers.
        eukaryote_lines: Independent eukaryote lines that hold, in reporting order.
        linked_ssu: A eukaryotic 18S locus sits on a bin contig.
        tiara_eukaryote: Tiara calls most of the bin eukaryotic or organelle.
        tiara_prokaryote: Tiara calls at least 0.8 of the bin prokaryotic.
        tiara_organelle: Tiara calls most of a bin of at most 2 Mb organelle and
            GVClass does not call it prokaryotic.
        checkm2_complete: CheckM2 completeness is at least 10 percent.
        viral: GVClass or whole-contig geNomad calls support a viral route.
    """

    quickclade: str | None
    gvclass: str | None
    checkeuk_resolved: bool
    checkeuk_fallback: bool
    eukaryote_lines: tuple[str, ...]
    linked_ssu: bool
    tiara_eukaryote: bool
    tiara_prokaryote: bool
    tiara_organelle: bool
    checkm2_complete: bool
    viral: bool

    @property
    def eukaryotic(self) -> bool:
        """Whether at least two eukaryote lines or the CheckEUK fallback hold."""
        return len(self.eukaryote_lines) >= 2 or self.checkeuk_fallback

    @property
    def has_eukaryotic_evidence(self) -> bool:
        """Whether any single eukaryote line holds."""
        return bool(self.eukaryote_lines)

    @property
    def has_viral_conflicting_eukaryote_line(self) -> bool:
        """Whether a CheckEUK lineage, GVClass ``d_EUK`` or 18S locus holds.

        Tiara has no virus class and labels giant-virus contigs eukaryotic, so its
        eukaryote line never conflicts with viral support.
        """
        return any(line != TIARA_LINE for line in self.eukaryote_lines)

    @property
    def tiara_checkm2_prokaryotic(self) -> bool:
        """Whether Tiara and CheckM2 together grant a prokaryotic route."""
        return self.tiara_prokaryote and self.checkm2_complete and not self.linked_ssu

    @property
    def prokaryotic(self) -> bool:
        """Whether any prokaryote evidence alone grants the prokaryotic route."""
        return self.gvclass == "prokaryotic" or self.tiara_checkm2_prokaryotic

    @property
    def organelle(self) -> bool:
        """Whether GVClass or Tiara calls the bin organelle."""
        return self.gvclass == "organelle" or self.tiara_organelle


def route_bins(
    sample_path: Path,
    bin_dir: Path,
    evidence_dirs: Mapping[str, Path],
    output_dir: Path,
) -> Path:
    """Copy each bin to one route and write the complete evidence ledger.

    Args:
        sample_path: JSON file with ``sample_id`` and ``run_id``.
        bin_dir: Directory containing the QuickBin ``*.fa`` files.
        evidence_dirs: Native result roots keyed by tool name. The CheckM2
            table may be absent, which means CheckM2 was skipped.
        output_dir: New directory for the evidence table and routed FASTAs.

    Returns:
        The supplied output directory.

    Raises:
        ValueError: If evidence is malformed, incomplete, or scientifically
            unsupported.
    """
    if output_dir.exists():
        raise ValueError(f"routing output already exists: {output_dir}")
    paths = _report_paths(evidence_dirs)
    sample = _sample_identity(sample_path)
    bins, contig_lengths = _bins(bin_dir)
    if not bins:
        _write_outputs(output_dir, {}, {}, [])
        _write_manifest(sample_path, bin_dir, evidence_dirs, output_dir)
        return output_dir
    native = _read_native(paths, set(bins), contig_lengths)
    rows: list[Row] = []
    decisions: dict[str, Decision] = {}
    for bin_id in sorted(bins):
        evidence = _bin_evidence(native, bin_id, contig_lengths[bin_id])
        decision = _decide(evidence)
        decisions[bin_id] = decision
        rows.append(_evidence_row(sample, bin_id, evidence, decision))
    _write_outputs(output_dir, bins, decisions, rows)
    _write_manifest(sample_path, bin_dir, evidence_dirs, output_dir)
    return output_dir


def _report_paths(evidence_dirs: Mapping[str, Path]) -> dict[str, Path]:
    missing = set(REPORTS).difference(evidence_dirs)
    if missing:
        raise ValueError(f"missing evidence directories: {sorted(missing)}")
    return {name: evidence_dirs[name] / report for name, report in REPORTS.items()}


def _read_native(
    paths: Mapping[str, Path],
    bin_ids: set[str],
    contig_lengths: Mapping[str, Mapping[str, int]],
) -> NativeEvidence:
    contig_to_bin = {
        contig_id: bin_id
        for bin_id, lengths in contig_lengths.items()
        for contig_id in lengths
    }
    quickclade = _indexed_report(
        paths["quickclade"],
        ("#QueryName", "Q_Bases", "lineage", "ConfLevel", "Confidence"),
        "#QueryName",
    )
    checkeuk = _indexed_report(
        paths["checkeuk"], ("genome", "status", "lineage"), "genome"
    )
    gvclass = _indexed_report(
        paths["gvclass"],
        ("query", "taxonomy_majority", "taxonomy_confidence", "domain"),
        "query",
    )
    for name, report in (
        ("QuickClade", quickclade),
        ("CheckEUK", checkeuk),
        ("GVClass", gvclass),
    ):
        if set(report) != bin_ids:
            raise ValueError(f"{name} identifiers do not exactly cover physical bins")
    _validate_native_labels(quickclade, checkeuk, gvclass)
    voting_markers = _checkeuk_voting_markers(paths["checkeuk"])
    if set(voting_markers) != bin_ids:
        raise ValueError("CheckEUK detail identifiers do not exactly cover bins")
    return NativeEvidence(
        quickclade,
        checkeuk,
        gvclass,
        voting_markers,
        _linked_ssu(paths["ssuextract"], contig_to_bin),
        _linked_viruses(paths["viral"], contig_to_bin),
        _tiara_classes(paths["tiara"]),
        _checkm2_calls(paths["checkm2"], bin_ids),
    )


def _bin_evidence(
    native: NativeEvidence, bin_id: str, contig_lengths: Mapping[str, int]
) -> BinEvidence:
    viruses = native.linked_viruses.get(bin_id, [])
    viral_contigs = {row["seq_name"] for row in viruses if _is_whole_contig_virus(row)}
    return BinEvidence(
        native.checkeuk[bin_id],
        native.voting_markers[bin_id],
        native.quickclade[bin_id],
        native.gvclass[bin_id],
        tuple(native.linked_ssu.get(bin_id, ())),
        tuple(viruses),
        _composition(contig_lengths, native.tiara, viral_contigs),
        native.checkm2.get(bin_id),
    )


def _checkeuk_voting_markers(report_path: Path) -> dict[str, int]:
    detail = report_path.with_name(CHECKEUK_DETAIL.name)
    required = ("genome", CHECKEUK_VOTING_MARKERS)
    counts: dict[str, int] = {}
    for bin_id, row in _indexed_report(detail, required, "genome").items():
        text = row[CHECKEUK_VOTING_MARKERS].strip()
        try:
            counts[bin_id] = int(text) if text else 0
        except ValueError:
            raise ValueError(
                f"{detail} has a non-integer voting marker count for {bin_id}: {text!r}"
            ) from None
    return counts


def _sample_identity(path: Path) -> dict[str, str]:
    value = json.loads(path.read_text(encoding="utf-8"))
    identity = {name: value.get(name) for name in ("sample_id", "run_id")}
    if not all(isinstance(item, str) and item.strip() for item in identity.values()):
        raise ValueError("sample JSON requires nonempty sample_id and run_id strings")
    return {name: str(item) for name, item in identity.items()}


def _bins(directory: Path) -> tuple[dict[str, Path], dict[str, dict[str, int]]]:
    """Return each bin's FASTA path and the base-pair length of each of its contigs."""
    bins: dict[str, Path] = {}
    contig_lengths: dict[str, dict[str, int]] = {}
    seen: set[str] = set()
    for path in sorted(directory.glob("*.fa")):
        lengths: dict[str, int] = {}
        for record in read_fasta(path):
            if record.identifier in lengths:
                raise ValueError(f"bin FASTA repeats a contig identifier: {path}")
            if record.identifier in seen:
                raise ValueError(f"contig occurs in multiple bins: {record.identifier}")
            lengths[record.identifier] = len(record.sequence)
        if not lengths:
            raise ValueError(f"bin FASTA contains no records: {path}")
        seen.update(lengths)
        bins[path.stem] = path
        contig_lengths[path.stem] = lengths
    return bins, contig_lengths


def _read_report(path: Path, required: tuple[str, ...]) -> list[Row]:
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


def _indexed_report(path: Path, required: tuple[str, ...], key: str) -> dict[str, Row]:
    result: dict[str, Row] = {}
    for row in _read_report(path, required):
        identifier = _bin_id(row[key])
        if identifier in result:
            raise ValueError(f"{path} repeats identifier {identifier}")
        result[identifier] = row
    return result


def _bin_id(raw: str) -> str:
    name = Path(raw).name
    for suffix in (".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if not name:
        raise ValueError("native evidence contains an empty bin identifier")
    return name


def _validate_native_labels(
    quickclade: Mapping[str, Row],
    checkeuk: Mapping[str, Row],
    gvclass: Mapping[str, Row],
) -> None:
    for row in quickclade.values():
        _quickclade_domain(row)
    for row in checkeuk.values():
        if row["status"] not in {"ok", "filtered"}:
            raise ValueError(f"unsupported CheckEUK status: {row['status']}")
    for row in gvclass.values():
        domain = _gvclass_domain(row)
        if domain not in GVCLASS_DOMAINS and domain not in UNCLASSIFIED_GVCLASS:
            raise ValueError(f"unsupported GVClass domain: {domain}")


def _linked_ssu(path: Path, contig_to_bin: Mapping[str, str]) -> dict[str, list[Row]]:
    required = (
        "name",
        "model",
        "contig_name",
        "taxonomy",
        "taxonomy_domain",
    )
    result: dict[str, list[Row]] = {}
    for row in _read_report(path, required):
        bin_id = contig_to_bin.get(row["contig_name"])
        if (
            bin_id is not None
            and row["model"] == EUKARYOTIC_SSU_MODEL
            and _is_eukaryotic_locus(row)
        ):
            result.setdefault(bin_id, []).append(row)
    return result


def _is_eukaryotic_locus(row: Row) -> bool:
    """Accept SSUextract's domain call or either reference-derived lineage.

    SSUextract 1.2.2 leaves ``taxonomy_domain`` as ``Unclassified`` for every
    18S locus, so the lineage columns carry the eukaryotic evidence. Older
    tables lack those columns, which then count as empty.
    """
    return row["taxonomy_domain"] == EUKARYOTIC_SSU_DOMAIN or any(
        row.get(column, "").startswith(EUKARYOTIC_SSU_DOMAIN)
        for column in ("reference_taxonomy", "blast_candidate_taxonomy")
    )


def _linked_viruses(
    path: Path, contig_to_bin: Mapping[str, str]
) -> dict[str, list[Row]]:
    required = ("seq_name", "virus_score", "taxonomy")
    result: dict[str, list[Row]] = {}
    for row in _read_report(path, required):
        contig_id = row["seq_name"].split(PROVIRUS_SEPARATOR, maxsplit=1)[0]
        bin_id = contig_to_bin.get(contig_id)
        if bin_id is not None:
            result.setdefault(bin_id, []).append(row)
    return result


def _tiara_classes(path: Path) -> dict[str, str]:
    """Map each contig Tiara classified to eukaryote, prokaryote, organelle or other.

    Tiara omits contigs below its minimum length, so callers treat a missing
    contig as unclassified.
    """
    required = ("sequence_id", "class_fst_stage", "class_snd_stage")
    classes: dict[str, str] = {}
    for row in _read_report(path, required):
        words = row["sequence_id"].split()
        if not words:
            raise ValueError(f"{path} has an empty sequence_id")
        contig_id = words[0]
        if contig_id in classes:
            raise ValueError(f"{path} repeats identifier {contig_id}")
        classes[contig_id] = _tiara_class(row)
    return classes


def _tiara_class(row: Row) -> str:
    first = row["class_fst_stage"].strip()
    if first not in TIARA_FIRST_STAGE:
        raise ValueError(f"unsupported Tiara class: {first}")
    if row["class_snd_stage"].strip() in TIARA_ORGANELLE_SECOND_STAGE:
        return "organelle"
    return TIARA_FIRST_STAGE[first]


def _checkm2_calls(path: Path, bin_ids: set[str]) -> dict[str, Checkm2Call]:
    """Read CheckM2 scores; a missing table means CheckM2 was skipped."""
    if not path.is_file():
        return {}
    calls: dict[str, Checkm2Call] = {}
    for bin_id, row in _indexed_report(
        path, ("Name", "Completeness", "Contamination"), "Name"
    ).items():
        if bin_id not in bin_ids:
            raise ValueError(
                f"CheckM2 reports {bin_id}, which is not a physical bin: {path}"
            )
        calls[bin_id] = Checkm2Call(
            _percent(row["Completeness"], "CheckM2 Completeness", path),
            _percent(row["Contamination"], "CheckM2 Contamination", path),
        )
    return calls


def _percent(text: str, label: str, path: Path) -> float:
    try:
        value = float(text)
    except ValueError:
        raise ValueError(f"{path} has a non-numeric {label}: {text!r}") from None
    if not math.isfinite(value):
        raise ValueError(f"{path} has a non-finite {label}: {text!r}")
    return value


def _composition(
    contig_lengths: Mapping[str, int],
    tiara_classes: Mapping[str, str],
    viral_contigs: set[str],
) -> BinComposition:
    bin_bp = sum(contig_lengths.values())
    tiara_bp = dict.fromkeys(TIARA_CLASSES, 0)
    viral_bp = 0
    duplicated_bp = 0
    for contig_id, length in contig_lengths.items():
        tiara_bp[tiara_classes.get(contig_id, "unclassified")] += length
        if contig_id in viral_contigs:
            viral_bp += length
        if DUPLICATED_CONTIG_MARKER in contig_id:
            duplicated_bp += length
    return BinComposition(
        bin_bp,
        tiara_bp["eukaryote"] / bin_bp,
        tiara_bp["prokaryote"] / bin_bp,
        tiara_bp["organelle"] / bin_bp,
        tiara_bp["unclassified"] / bin_bp,
        viral_bp / bin_bp,
        duplicated_bp / bin_bp,
    )


def _quickclade_domain(row: Row) -> str:
    first = row["lineage"].strip().split(";", maxsplit=1)[0]
    if not first.startswith(("d__", "sk__")):
        return ""
    if first.partition("__")[2].strip() == "":
        raise ValueError(f"malformed QuickClade domain: {first!r}")
    return first


def _gvclass_domain(row: Row) -> str:
    return row["taxonomy_majority"].strip().split(";", maxsplit=1)[0]


def _decide(evidence: BinEvidence) -> Decision:
    signals = _signals(evidence)
    outcome = _route(signals)
    return Decision(
        outcome.route,
        outcome.candidate_class,
        outcome.reason,
        _supporting(evidence, signals, outcome.route),
        _conflicting(evidence, signals),
    )


def _route(signals: Signals) -> Outcome:
    """Apply the routing rules in order; the first that matches wins."""
    rules: tuple[Callable[[Signals], Outcome | None], ...] = (
        _ambiguous_rule,
        _organelle_rule,
        _tiara_checkm2_rule,
        _eukaryote_rule,
        _viral_rule,
        _gvclass_prokaryote_rule,
    )
    for rule in rules:
        outcome = rule(signals)
        if outcome is not None:
            return outcome
    reason = (
        "eukaryotic evidence does not meet the combined support rule"
        if signals.has_eukaryotic_evidence
        else "no supported routing evidence"
    )
    return Outcome("unresolved", "unresolved", reason)


def _conflict() -> Outcome:
    return Outcome("unresolved", "conflicting", CONFLICT_REASON)


def _ambiguous_rule(signals: Signals) -> Outcome | None:
    if signals.gvclass != "ambiguous":
        return None
    return Outcome("unresolved", "conflicting", AMBIGUOUS_REASON)


def _organelle_rule(signals: Signals) -> Outcome | None:
    if not signals.organelle:
        return None
    reason = (
        GVCLASS_ORGANELLE_REASON
        if signals.gvclass == "organelle"
        else TIARA_ORGANELLE_REASON
    )
    contested = signals.eukaryotic or signals.viral
    return Outcome("unresolved", "conflicting" if contested else "unresolved", reason)


def _tiara_checkm2_rule(signals: Signals) -> Outcome | None:
    if not signals.tiara_checkm2_prokaryotic:
        return None
    if signals.viral:
        return _conflict()
    reason = (
        GVCLASS_PROKARYOTIC_REASON
        if signals.gvclass == "prokaryotic"
        else TIARA_CHECKM2_REASON
    )
    return Outcome("prokaryotic", "bacterial_archaeal_candidate", reason)


def _eukaryote_rule(signals: Signals) -> Outcome | None:
    if not signals.eukaryotic:
        return None
    if signals.viral or signals.gvclass == "prokaryotic" or signals.tiara_prokaryote:
        return _conflict()
    return Outcome("eukaryotic", "eukaryotic_candidate", _eukaryotic_reason(signals))


def _eukaryotic_reason(signals: Signals) -> str:
    if signals.checkeuk_fallback:
        return FALLBACK_REASON
    first, second = signals.eukaryote_lines[:2]
    return f"{first} and {second} support a eukaryotic route"


def _viral_rule(signals: Signals) -> Outcome | None:
    if not signals.viral:
        return None
    if signals.prokaryotic or signals.has_viral_conflicting_eukaryote_line:
        return _conflict()
    return Outcome("viral", "viral_candidate", VIRAL_REASON)


def _gvclass_prokaryote_rule(signals: Signals) -> Outcome | None:
    if signals.gvclass != "prokaryotic":
        return None
    return Outcome(
        "prokaryotic", "bacterial_archaeal_candidate", GVCLASS_PROKARYOTIC_REASON
    )


def _signals(evidence: BinEvidence) -> Signals:
    composition = evidence.composition
    checkeuk_lineage = evidence.checkeuk["lineage"].strip()
    checkeuk_resolved = evidence.checkeuk["status"] == "ok" and bool(checkeuk_lineage)
    gvclass = GVCLASS_DOMAINS.get(_gvclass_domain(evidence.gvclass))
    linked_ssu = bool(evidence.eukaryotic_ssu)
    tiara_eukaryote = composition.tiara_eukaryote_derived >= MIN_TIARA_FRACTION
    lines = tuple(
        name
        for name, holds in (
            ("CheckEUK", checkeuk_resolved),
            ("GVClass", gvclass == "eukaryotic"),
            (SSU_LINE, linked_ssu),
            (TIARA_LINE, tiara_eukaryote),
        )
        if holds
    )
    checkeuk_fallback = (
        evidence.checkeuk["status"] == "ok"
        and not checkeuk_lineage
        and evidence.checkeuk_voting_markers < MIN_CHECKEUK_VOTING_MARKERS
        and gvclass == "eukaryotic"
    )
    checkm2 = evidence.checkm2
    return Signals(
        quickclade=QUICKCLADE_DOMAINS.get(_quickclade_domain(evidence.quickclade)),
        gvclass=gvclass,
        checkeuk_resolved=checkeuk_resolved,
        checkeuk_fallback=checkeuk_fallback,
        eukaryote_lines=lines,
        linked_ssu=linked_ssu,
        tiara_eukaryote=tiara_eukaryote,
        tiara_prokaryote=composition.tiara_prokaryote >= TIARA_STRONG_PROKARYOTE,
        tiara_organelle=(
            composition.tiara_organelle >= MIN_TIARA_FRACTION
            and composition.bin_bp <= MAX_TIARA_ORGANELLE_BIN_BP
            and gvclass != "prokaryotic"
        ),
        checkm2_complete=(
            checkm2 is not None and checkm2.completeness >= MIN_CHECKM2_COMPLETENESS
        ),
        viral=gvclass == "viral"
        or composition.viral_contig >= MIN_VIRAL_CONTIG_FRACTION,
    )


def _supporting(
    evidence: BinEvidence,
    signals: Signals,
    route: Route,
) -> tuple[str, ...]:
    items: list[str] = []
    if signals.gvclass == "organelle":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    if signals.tiara_organelle:
        items.append(f"Tiara:{_fraction(evidence.composition.tiara_organelle)}")
    # Eukaryote labels on a bin that Tiara and CheckM2 call prokaryotic are
    # overridden evidence, so they are listed as conflicts instead.
    if signals.eukaryotic and route != "prokaryotic":
        items.extend(_eukaryotic_support(evidence, signals))
    if signals.prokaryotic:
        items.extend(_prokaryotic_support(evidence, signals))
    if signals.viral:
        items.extend(_viral_support(evidence, signals))
    return tuple(items)


def _eukaryotic_support(evidence: BinEvidence, signals: Signals) -> list[str]:
    items: list[str] = []
    if signals.checkeuk_resolved:
        items.append(f"CheckEUK:{evidence.checkeuk['lineage']}")
    if signals.gvclass == "eukaryotic":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    items.extend(
        f"SSU:{row['name']}@{row['contig_name']}" for row in evidence.eukaryotic_ssu
    )
    if signals.tiara_eukaryote:
        derived = evidence.composition.tiara_eukaryote_derived
        items.append(f"Tiara:{_fraction(derived)}")
    return items


def _prokaryotic_support(evidence: BinEvidence, signals: Signals) -> list[str]:
    items: list[str] = []
    if signals.gvclass == "prokaryotic":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
        if signals.quickclade == "prokaryotic":
            items.append(f"QuickClade:{evidence.quickclade['lineage']}")
    if signals.tiara_checkm2_prokaryotic and evidence.checkm2 is not None:
        items.append(f"Tiara:{_fraction(evidence.composition.tiara_prokaryote)}")
        items.append(f"CheckM2:{_checkm2_text(evidence.checkm2)}")
    return items


def _viral_support(evidence: BinEvidence, signals: Signals) -> list[str]:
    items: list[str] = []
    if signals.gvclass == "viral":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    if evidence.composition.viral_contig >= MIN_VIRAL_CONTIG_FRACTION:
        items.extend(
            f"geNomad:{row['seq_name']}"
            for row in evidence.genomad_viruses
            if _is_whole_contig_virus(row)
        )
    return items


def _is_whole_contig_virus(row: Row) -> bool:
    return PROVIRUS_SEPARATOR not in row["seq_name"]


def _conflicting(
    evidence: BinEvidence,
    signals: Signals,
) -> tuple[str, ...]:
    items = _quickclade_conflicts(evidence, signals)
    items.extend(_cross_route_conflicts(evidence, signals))
    if signals.tiara_checkm2_prokaryotic:
        items.extend(_eukaryote_labels(evidence, signals))
    if 0 < evidence.composition.viral_contig < MIN_VIRAL_CONTIG_FRACTION:
        items.append(f"geNomad-minority:{_fraction(evidence.composition.viral_contig)}")
    return tuple(items)


def _quickclade_conflicts(evidence: BinEvidence, signals: Signals) -> list[str]:
    items: list[str] = []
    quickclade_domain = _quickclade_domain(evidence.quickclade)
    if quickclade_domain and quickclade_domain not in QUICKCLADE_DOMAINS:
        items.append(f"QuickClade-unresolved:{evidence.quickclade['lineage']}")
    if signals.eukaryotic and signals.quickclade == "prokaryotic":
        items.append(f"QuickClade:{evidence.quickclade['lineage']}")
    return items


def _cross_route_conflicts(evidence: BinEvidence, signals: Signals) -> list[str]:
    items: list[str] = []
    if signals.viral and signals.prokaryotic:
        items.append("prokaryotic evidence conflicts with viral evidence")
    if signals.viral and signals.has_viral_conflicting_eukaryote_line:
        items.append("eukaryotic evidence conflicts with viral evidence")
    if signals.eukaryotic and signals.gvclass == "prokaryotic":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    if (
        signals.eukaryotic
        and signals.tiara_prokaryote
        and not signals.tiara_checkm2_prokaryotic
    ):
        items.append(f"Tiara:{_fraction(evidence.composition.tiara_prokaryote)}")
    if signals.gvclass == "ambiguous":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    if signals.organelle and (signals.eukaryotic or signals.viral):
        items.append("organelle evidence conflicts with supported routing evidence")
    return items


def _eukaryote_labels(evidence: BinEvidence, signals: Signals) -> list[str]:
    """List CheckEUK and GVClass eukaryote labels that a prokaryotic route overrides."""
    items: list[str] = []
    if signals.checkeuk_resolved:
        items.append(f"CheckEUK:{evidence.checkeuk['lineage']}")
    if signals.gvclass == "eukaryotic":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    return items


def _fraction(value: float) -> str:
    return f"{value:.4f}"


def _checkm2_text(call: Checkm2Call) -> str:
    return f"{call.completeness:.2f}/{call.contamination:.2f}"


def _evidence_row(
    sample: Mapping[str, str],
    bin_id: str,
    evidence: BinEvidence,
    decision: Decision,
) -> Row:
    ssu = [
        {
            "name": row["name"],
            "contig_name": row["contig_name"],
            "taxonomy": row["taxonomy"],
            "taxonomy_domain": row["taxonomy_domain"],
            "reference_taxonomy": row.get("reference_taxonomy", ""),
            "blast_candidate_taxonomy": row.get("blast_candidate_taxonomy", ""),
        }
        for row in evidence.eukaryotic_ssu
    ]
    viruses = [
        {
            "seq_name": row["seq_name"],
            "virus_score": row["virus_score"],
            "taxonomy": row["taxonomy"],
        }
        for row in evidence.genomad_viruses
    ]
    taxonomy, taxonomy_source = _taxonomy(evidence, decision)
    composition = evidence.composition
    checkm2 = evidence.checkm2
    return {
        **sample,
        "bin_id": bin_id,
        "route": decision.route,
        "candidate_class": decision.candidate_class,
        "reason": decision.reason,
        "checkeuk_status": evidence.checkeuk["status"],
        "checkeuk_lineage": evidence.checkeuk["lineage"],
        "taxonomy": taxonomy,
        "taxonomy_source": taxonomy_source,
        "quickclade_domain": _quickclade_domain(evidence.quickclade),
        "quickclade_lineage": evidence.quickclade["lineage"],
        "quickclade_confidence": evidence.quickclade["Confidence"],
        "gvclass_domain": _gvclass_domain(evidence.gvclass),
        "gvclass_taxonomy_majority": evidence.gvclass["taxonomy_majority"],
        "gvclass_taxonomy_confidence": evidence.gvclass["taxonomy_confidence"],
        "linked_eukaryotic_ssu": _json(ssu),
        "linked_genomad_viruses": _json(viruses),
        "supporting_evidence": _json(decision.supporting),
        "conflicting_evidence": _json(decision.conflicting),
        "bin_bp": str(composition.bin_bp),
        "tiara_eukaryote_fraction": _fraction(composition.tiara_eukaryote),
        "tiara_prokaryote_fraction": _fraction(composition.tiara_prokaryote),
        "tiara_organelle_fraction": _fraction(composition.tiara_organelle),
        "tiara_unclassified_fraction": _fraction(composition.tiara_unclassified),
        "viral_contig_fraction": _fraction(composition.viral_contig),
        "checkm2_completeness": ""
        if checkm2 is None
        else f"{checkm2.completeness:.2f}",
        "checkm2_contamination": ""
        if checkm2 is None
        else f"{checkm2.contamination:.2f}",
        "duplicated_contig_fraction": _fraction(composition.duplicated_contig),
    }


def _taxonomy(evidence: BinEvidence, decision: Decision) -> tuple[str, str]:
    lineage = evidence.checkeuk["lineage"].strip()
    if evidence.checkeuk["status"] == "ok" and lineage:
        return lineage, "checkeuk"
    if decision.reason == FALLBACK_REASON:
        return evidence.gvclass["taxonomy_majority"], "gvclass_fallback"
    return "", ""


def _json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _write_outputs(
    output_dir: Path,
    bins: Mapping[str, Path],
    decisions: Mapping[str, Decision],
    rows: list[Row],
) -> None:
    output_dir.mkdir(parents=True)
    for route in ROUTES:
        (output_dir / route).mkdir()
    for bin_id, source in bins.items():
        destination = output_dir / decisions[bin_id].route / f"{bin_id}.fna"
        shutil.copyfile(source, destination)
    with (output_dir / "evidence.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=EVIDENCE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_manifest(
    sample_path: Path,
    bin_dir: Path,
    evidence_dirs: Mapping[str, Path],
    output_dir: Path,
) -> None:
    command = [
        "protist-meta",
        "route",
        "--sample-json",
        str(sample_path),
        "--bins-dir",
        str(bin_dir),
    ]
    for name in REPORTS:
        command.extend((f"--{name}-dir", str(evidence_dirs[name])))
    command.extend(("--output-dir", str(output_dir)))
    manifest = {
        "stage": "routing",
        "status": "completed",
        "tool": "protist-meta",
        "tool_version": version("protist-meta"),
        "command": command,
        "outputs": {
            "evidence": "evidence.tsv",
            **{route: route for route in ROUTES},
        },
    }
    (output_dir / "stage.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
