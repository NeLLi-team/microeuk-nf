"""Route bins to one gene-calling path from independent native evidence."""

import csv
import json
import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import Literal

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
    "d_EUK": "eukaryotic",
    "d_NCLDV": "viral",
    "d_PPV": "viral",
    "d_MIRUS": "viral",
    "d_EUK-pEVE": "ambiguous",
}
UNCLASSIFIED_GVCLASS = {"", "d_", "nd"}
EUKARYOTIC_SSU_MODEL = "RF01960"
PROVIRUS_SEPARATOR = "|provirus_"

EVIDENCE_FIELDS = (
    "sample_id",
    "run_id",
    "bin_id",
    "route",
    "candidate_class",
    "reason",
    "checkeuk_status",
    "checkeuk_lineage",
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
)


@dataclass(frozen=True, slots=True)
class BinEvidence:
    """Native evidence associated with one physical bin."""

    checkeuk: Row
    quickclade: Row
    gvclass: Row
    eukaryotic_ssu: tuple[Row, ...]
    genomad_viruses: tuple[Row, ...]


@dataclass(frozen=True, slots=True)
class Decision:
    """One exclusive route and its catalog interpretation."""

    route: Route
    candidate_class: CandidateClass
    reason: str
    supporting: tuple[str, ...]
    conflicting: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class Signals:
    """Normalized support flags used by the routing policy."""

    checkeuk_resolved: bool
    quickclade: str | None
    gvclass: str | None
    eukaryotic: bool
    prokaryotic: bool
    viral: bool


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
        evidence_dirs: Native result roots keyed by tool name.
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
    bins, contig_to_bin = _bins(bin_dir)
    if not bins:
        _write_outputs(output_dir, {}, {}, [])
        _write_manifest(sample_path, bin_dir, evidence_dirs, output_dir)
        return output_dir
    bin_ids = set(bins)
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
    linked_ssu = _linked_ssu(paths["ssuextract"], contig_to_bin)
    linked_viruses = _linked_viruses(paths["viral"], contig_to_bin)
    rows: list[Row] = []
    decisions: dict[str, Decision] = {}
    for bin_id in sorted(bins):
        evidence = BinEvidence(
            checkeuk[bin_id],
            quickclade[bin_id],
            gvclass[bin_id],
            tuple(linked_ssu.get(bin_id, ())),
            tuple(linked_viruses.get(bin_id, ())),
        )
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


def _sample_identity(path: Path) -> dict[str, str]:
    value = json.loads(path.read_text(encoding="utf-8"))
    identity = {name: value.get(name) for name in ("sample_id", "run_id")}
    if not all(isinstance(item, str) and item.strip() for item in identity.values()):
        raise ValueError("sample JSON requires nonempty sample_id and run_id strings")
    return {name: str(item) for name, item in identity.items()}


def _bins(directory: Path) -> tuple[dict[str, Path], dict[str, str]]:
    bins: dict[str, Path] = {}
    contig_to_bin: dict[str, str] = {}
    for path in sorted(directory.glob("*.fa")):
        bin_id = path.stem
        identifiers = [record.identifier for record in read_fasta(path)]
        if not identifiers:
            raise ValueError(f"bin FASTA contains no records: {path}")
        if len(identifiers) != len(set(identifiers)):
            raise ValueError(f"bin FASTA repeats a contig identifier: {path}")
        for contig_id in identifiers:
            if contig_id in contig_to_bin:
                raise ValueError(f"contig occurs in multiple bins: {contig_id}")
            contig_to_bin[contig_id] = bin_id
        bins[bin_id] = path
    return bins, contig_to_bin


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
    required = ("name", "model", "contig_name", "taxonomy", "taxonomy_domain")
    result: dict[str, list[Row]] = {}
    for row in _read_report(path, required):
        bin_id = contig_to_bin.get(row["contig_name"])
        if (
            bin_id is not None
            and row["model"] == EUKARYOTIC_SSU_MODEL
            and row["taxonomy_domain"] == "Eukaryota"
        ):
            result.setdefault(bin_id, []).append(row)
    return result


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
    supporting = _supporting(evidence, signals)
    conflicting = _conflicting(evidence, signals)
    route: Route
    candidate_class: CandidateClass
    if signals.gvclass == "ambiguous":
        route = "unresolved"
        candidate_class = "conflicting"
        reason = "GVClass reports an ambiguous eukaryote-virus class"
    elif signals.eukaryotic and (signals.viral or signals.gvclass == "prokaryotic"):
        route = "unresolved"
        candidate_class = "conflicting"
        reason = "supported evidence conflicts across routes"
    elif signals.eukaryotic:
        route = "eukaryotic"
        candidate_class = "eukaryotic_candidate"
        reason = (
            "CheckEUK and GVClass support a eukaryotic route"
            if signals.gvclass == "eukaryotic"
            else (
                "CheckEUK and a physically linked eukaryotic 18S locus support "
                "a eukaryotic route"
            )
        )
    elif signals.viral and (
        signals.prokaryotic or signals.checkeuk_resolved or evidence.eukaryotic_ssu
    ):
        route = "unresolved"
        candidate_class = "conflicting"
        reason = "supported evidence conflicts across routes"
    elif signals.viral:
        route = "viral"
        candidate_class = "viral_candidate"
        reason = (
            "GVClass or linked whole-contig geNomad evidence supports a viral route"
        )
    elif signals.prokaryotic:
        route = "prokaryotic"
        candidate_class = "bacterial_archaeal_candidate"
        reason = (
            "GVClass or QuickClade supports a prokaryotic route without a "
            "supported eukaryotic or viral conflict"
        )
    else:
        route = "unresolved"
        candidate_class = "unresolved"
        reason = (
            "eukaryotic evidence does not meet the combined support rule"
            if (
                signals.checkeuk_resolved
                or signals.gvclass == "eukaryotic"
                or evidence.eukaryotic_ssu
            )
            else "no supported routing evidence"
        )
    return Decision(route, candidate_class, reason, supporting, conflicting)


def _signals(evidence: BinEvidence) -> Signals:
    checkeuk_lineage = evidence.checkeuk["lineage"].strip()
    checkeuk_resolved = evidence.checkeuk["status"] == "ok" and bool(checkeuk_lineage)
    quickclade = QUICKCLADE_DOMAINS.get(_quickclade_domain(evidence.quickclade))
    gvclass_domain = _gvclass_domain(evidence.gvclass)
    gvclass = GVCLASS_DOMAINS.get(gvclass_domain)
    linked_ssu = bool(evidence.eukaryotic_ssu)
    linked_virus = any(_is_whole_contig_virus(row) for row in evidence.genomad_viruses)
    euk_supported = checkeuk_resolved and (gvclass == "eukaryotic" or linked_ssu)
    prok_supported = gvclass == "prokaryotic" or quickclade == "prokaryotic"
    virus_supported = gvclass == "viral" or linked_virus
    return Signals(
        checkeuk_resolved,
        quickclade,
        gvclass,
        euk_supported,
        prok_supported,
        virus_supported,
    )


def _supporting(
    evidence: BinEvidence,
    signals: Signals,
) -> tuple[str, ...]:
    items: list[str] = []
    if signals.eukaryotic:
        items.append(f"CheckEUK:{evidence.checkeuk['lineage']}")
        if _gvclass_domain(evidence.gvclass) == "d_EUK":
            items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
        items.extend(
            f"SSU:{row['name']}@{row['contig_name']}" for row in evidence.eukaryotic_ssu
        )
    if signals.prokaryotic:
        if GVCLASS_DOMAINS.get(_gvclass_domain(evidence.gvclass)) == "prokaryotic":
            items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
        if signals.quickclade == "prokaryotic":
            items.append(f"QuickClade:{evidence.quickclade['lineage']}")
    if signals.viral:
        if GVCLASS_DOMAINS.get(_gvclass_domain(evidence.gvclass)) == "viral":
            items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
        items.extend(
            f"geNomad:{row['seq_name']}"
            for row in evidence.genomad_viruses
            if _is_whole_contig_virus(row)
        )
    return tuple(items)


def _is_whole_contig_virus(row: Row) -> bool:
    return PROVIRUS_SEPARATOR not in row["seq_name"]


def _conflicting(
    evidence: BinEvidence,
    signals: Signals,
) -> tuple[str, ...]:
    items: list[str] = []
    quickclade_domain = _quickclade_domain(evidence.quickclade)
    if quickclade_domain and quickclade_domain not in QUICKCLADE_DOMAINS:
        items.append(f"QuickClade-unresolved:{evidence.quickclade['lineage']}")
    if signals.eukaryotic and signals.quickclade == "prokaryotic":
        items.append(f"QuickClade:{evidence.quickclade['lineage']}")
    if signals.viral and signals.prokaryotic:
        items.append("prokaryotic evidence conflicts with viral evidence")
    if signals.viral and (
        evidence.checkeuk["lineage"].strip() or evidence.eukaryotic_ssu
    ):
        items.append("eukaryotic evidence conflicts with viral evidence")
    if signals.eukaryotic and signals.gvclass == "prokaryotic":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    if signals.gvclass == "ambiguous":
        items.append(f"GVClass:{evidence.gvclass['taxonomy_majority']}")
    return tuple(items)


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
    return {
        **sample,
        "bin_id": bin_id,
        "route": decision.route,
        "candidate_class": decision.candidate_class,
        "reason": decision.reason,
        "checkeuk_status": evidence.checkeuk["status"],
        "checkeuk_lineage": evidence.checkeuk["lineage"],
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
    }


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
