"""Export portable native evidence without changing biological classifications."""

import csv
import hashlib
import io
import json
import re
import sqlite3
from pathlib import Path

__all__ = ["export_visualization"]

type Row = dict[str, object]

TABLES = (
    "read_stats",
    "assemblies",
    "bins",
    "contigs",
    "memberships",
    "qc",
    "taxonomy",
    "ssu",
    "viruses",
    "phenotypes",
    "stages",
    "artifacts",
)
METHODS = (
    "routing",
    "checkeuk",
    "checkeuk_detail",
    "gvclass",
    "checkm2",
    "checkm1",
    "gtdbtk",
    "quickclade",
)
VERSION_FIELDS = ("tool_name", "tool_version", "database_name", "database_version")


def export_visualization(database: Path, output_dir: Path) -> list[Path]:
    """Write one digest-bound JSON per native assembly, plus a file index.

    Small registered native reports are read only after digest verification.
    Missing completed-stage QC or taxonomy reports stop the export, as do
    changed artifacts. Other missing evidence and coverage remain unavailable.
    """
    database = database.resolve(strict=True)
    catalog_digest = _digest(database)
    bundles = []
    with sqlite3.connect(f"{database.as_uri()}?mode=ro", uri=True) as connection:
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only = ON")
        metadata = [
            dict(row) for row in connection.execute("SELECT * FROM catalog_metadata")
        ]
        if len(metadata) != 1:
            raise ValueError("expected exactly one catalog metadata row")
        samples = {
            row["sample_id"]: dict(row)
            for row in connection.execute("SELECT * FROM samples")
        }
        for run in connection.execute("SELECT * FROM runs ORDER BY sample_id, run_id"):
            tables = _read_tables(connection, run)
            bundles.extend(
                _bundle(
                    samples[run["sample_id"]],
                    dict(run),
                    metadata[0],
                    assembly,
                    tables,
                    {
                        "source_id": "catalog",
                        "path": str(database),
                        "sha256": catalog_digest,
                        "kind": "catalog",
                    },
                )
                for assembly in tables["assemblies"] or [None]
            )
    if _digest(database) != catalog_digest:
        raise ValueError("catalog SHA256 changed during export")
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    index = []
    for bundle in bundles:
        path = output_dir / f"{bundle['bundle_id']}.json"
        _write_json(path, bundle)
        paths.append(path)
        index.append(
            {
                "bundle_id": bundle["bundle_id"],
                "path": path.name,
                "sha256": _digest(path),
                "sample_id": bundle["sample_id"],
                "run_id": bundle["run_id"],
                "assembly_id": bundle["assembly_id"],
            }
        )
    _write_json(output_dir / "index.json", {"schema_version": 1, "bundles": index})
    return paths


def _read_tables(
    connection: sqlite3.Connection, run: sqlite3.Row
) -> dict[str, list[Row]]:
    # Table names are this module's fixed allowlist; external identifiers are bound.
    return {
        table: [
            dict(row)
            for row in connection.execute(
                f'SELECT * FROM "{table}" '  # noqa: S608  # Table is from the fixed module allowlist.
                "WHERE sample_id=? AND run_id=? ORDER BY rowid",
                (run["sample_id"], run["run_id"]),
            )
        ]
        for table in TABLES
    }


def _bundle(
    sample: Row,
    run: Row,
    metadata: Row,
    assembly: Row | None,
    tables: dict[str, list[Row]],
    catalog: Row,
) -> Row:
    assembly_id = assembly["assembly_id"] if assembly else None
    scoped = _scope_tables(tables, assembly_id)
    stages = {row["stage_id"]: row for row in tables["stages"]}
    sources = [catalog]
    artifacts = {row["artifact_id"]: row for row in tables["artifacts"]}
    native = _native_evidence(artifacts, stages, sources)
    _scope_evidence(native, {row["bin_id"] for row in scoped["bins"]})
    scoped["qc"] = [
        {
            **row,
            **_version(stages[row["stage_id"]]),
            "status": stages[row["stage_id"]]["status"],
        }
        for row in scoped["qc"]
    ]
    if assembly:
        assembly = {**assembly, "inventory_status": _inventory(assembly, scoped)}
    mapping = _mapping(scoped, artifacts, stages, sources)
    if assembly and assembly["inventory_status"] == "partial":
        mapping["status"] = "partial"
        mapping["status_reason"] = "partial_assembly_inventory"
        mapping["without_qualifying_alignment_reads"] = None
    referenced = {row.get("artifact_id") for rows in scoped.values() for row in rows}
    for identifier in sorted(referenced - {None}):
        if identifier in artifacts:
            _receipt(artifacts[identifier], sources, "catalog_record")
    bundle = {
        "schema_version": 1,
        "sample_id": run["sample_id"],
        "run_id": run["run_id"],
        "assembly_id": assembly_id,
        "sample": sample,
        "run": {**run, **metadata, "stages": tables["stages"]},
        "assembly": assembly,
        **scoped,
        "evidence": native,
        "mapping": mapping,
        "sources": sources,
    }
    digest = hashlib.sha256(_encode(bundle)).hexdigest()
    return {**bundle, "bundle_id": f"bundle_{digest}", "content_sha256": digest}


def _scope_tables(
    tables: dict[str, list[Row]], assembly_id: object
) -> dict[str, list[Row]]:
    scoped = {
        name: tables[name]
        for name in TABLES
        if name not in {"assemblies", "stages", "artifacts"}
    }
    for name in ("bins", "contigs"):
        scoped[name] = [
            row for row in tables[name] if row["assembly_id"] == assembly_id
        ]
    bins = {row["bin_id"] for row in scoped["bins"]}
    contigs = {row["contig_id"] for row in scoped["contigs"]}
    scoped["memberships"] = [
        row for row in tables["memberships"] if row["bin_id"] in bins
    ]
    seen = set()
    for row in scoped["memberships"]:
        if row["contig_id"] not in contigs or row["contig_id"] in seen:
            raise ValueError(
                "memberships must be disjoint and refer to assembly contigs"
            )
        seen.add(row["contig_id"])
    for name in ("qc", "taxonomy"):
        scoped[name] = [
            row
            for row in tables[name]
            if (row["target_kind"] == "bin" and row["target_id"] in bins)
            or (row["target_kind"] == "contig" and row["target_id"] in contigs)
            or (row["target_kind"] == "assembly" and row["target_id"] == assembly_id)
        ]
    for name in ("ssu", "viruses"):
        scoped[name] = [row for row in tables[name] if row["contig_id"] in contigs]
    scoped["phenotypes"] = [
        row for row in tables["phenotypes"] if row["bin_id"] in bins
    ]
    return scoped


def _inventory(assembly: Row, scoped: dict[str, list[Row]]) -> str:
    count = len(scoped["contigs"])
    span = sum(row["length_bp"] for row in scoped["contigs"])
    if count > assembly["contig_count"] or span > assembly["assembly_length_bp"]:
        raise ValueError("contig inventory exceeds assembly totals")
    if count == assembly["contig_count"] and span == assembly["assembly_length_bp"]:
        return "complete"
    return "partial"


def _version(stage: Row) -> Row:
    return {key: stage.get(key) for key in VERSION_FIELDS}


def _native_evidence(
    artifacts: dict[str, Row], stages: dict[str, Row], sources: list[Row]
) -> Row:
    evidence = {}
    for method in METHODS:
        stage_id = (
            "binning" if method == "quickclade" else method.removesuffix("_detail")
        )
        stage = stages.get(stage_id, {})
        evidence[method] = {
            **_version(stage),
            "status": stage.get("status", "unassessed"),
            "stage_status": stage.get("status", "unassessed"),
            "status_reason": stage.get("status_reason"),
            "rows": [],
            "source_ids": [],
        }
    seen = set()
    for artifact in artifacts.values():
        if artifact["kind"] not in {"qc_report", "taxonomy_report"}:
            continue
        receipt = _receipt(artifact, sources, "unavailable")
        content = _read_artifact(artifact, receipt)
        stage = stages.get(artifact["producer_stage_id"], {})
        if content is None and stage.get("status") == "completed":
            raise FileNotFoundError(
                "completed stage requires registered evidence: "
                f"{artifact['artifact_id']}"
            )
        method = _evidence_method(artifact)
        if method is None or (method, artifact["sha256"]) in seen:
            continue
        seen.add((method, artifact["sha256"]))
        report = evidence[method]
        report["source_ids"].append(receipt["source_id"])
        if content is None:
            report["status"] = "unavailable"
            report["status_reason"] = "registered_artifact_unavailable"
            continue
        report["rows"].extend(
            {**row, "source_id": receipt["source_id"]} for row in _parse_table(content)
        )
    return evidence


def _scope_evidence(evidence: Row, bins: set[str]) -> None:
    keys = {
        "routing": "bin_id",
        "checkeuk": "genome",
        "checkeuk_detail": "genome",
        "gvclass": "query",
        "checkm2": "Name",
        "checkm1": "Bin Id",
        "gtdbtk": "user_genome",
        "quickclade": "#QueryName",
    }
    for method, report in evidence.items():
        key = keys[method]
        if any(key not in row for row in report["rows"]):
            raise ValueError(f"native {method} report lacks identifier column {key}")
        report["rows"] = [
            row for row in report["rows"] if _native_bin_id(row[key]) in bins
        ]


def _native_bin_id(value: str) -> str:
    name = Path(value).name.removesuffix(".gz")
    return name.removesuffix(".fasta").removesuffix(".fna").removesuffix(".fa")


def _evidence_method(artifact: Row) -> str | None:
    identifier = str(artifact["artifact_id"])
    if identifier == "checkeuk.detail":
        return "checkeuk_detail"
    if identifier == "quickclade":
        return "quickclade"
    stage = str(artifact["producer_stage_id"])
    if stage in METHODS and artifact["kind"] in {"qc_report", "taxonomy_report"}:
        return stage
    return None


def _mapping(
    scoped: dict[str, list[Row]],
    artifacts: dict[str, Row],
    stages: dict[str, Row],
    sources: list[Row],
) -> Row:
    retained = [
        row
        for row in scoped["read_stats"]
        if row["stage_id"] == "read_qc" and row["retained_reads"] is not None
    ]
    denominator = retained[0]["retained_reads"] if len(retained) == 1 else None
    mapping = {
        "status": "unavailable",
        "input_reads": denominator,
        "assembly_mapped_reads": None,
        "bin_mapped_reads": None,
        "unbinned_mapped_reads": None,
        "without_qualifying_alignment_reads": None,
        "read_unit": "unpaired reads",
        "denominator_source": retained[0]["record_id"] if len(retained) == 1 else None,
        "per_contig": [],
        "per_bin": [],
        "source_ids": [],
        **_mapping_policy(stages),
    }
    identifiers = {row["coverage_artifact_id"] for row in scoped["contigs"]}
    if not identifiers or None in identifiers:
        mapping["status_reason"] = "coverage_not_registered"
        return mapping
    coverage = []
    for identifier in sorted(identifiers):
        artifact = artifacts[identifier]
        receipt = _receipt(artifact, sources, "unavailable")
        mapping["source_ids"].append(receipt["source_id"])
        content = _read_artifact(artifact, receipt)
        if content is None:
            mapping["status_reason"] = "coverage_artifact_unavailable"
            return mapping
        coverage.extend(_parse_table(content))
    if mapping["minimum_mapq"] is None or mapping["excluded_sam_flags"] is None:
        mapping["status_reason"] = "primary_alignment_policy_unavailable"
        return mapping
    per_contig = _coverage_counts(coverage, scoped["contigs"])
    by_contig = {row["contig_id"]: row for row in per_contig}
    by_bin = {
        row["bin_id"]: {"bin_id": row["bin_id"], "numreads": 0, "covbases": 0}
        for row in scoped["bins"]
    }
    for membership in scoped["memberships"]:
        for field in ("numreads", "covbases"):
            by_bin[membership["bin_id"]][field] += by_contig[membership["contig_id"]][
                field
            ]
    aligned = sum(row["numreads"] for row in per_contig)
    binned = sum(row["numreads"] for row in by_bin.values())
    if denominator is not None and aligned > denominator:
        raise ValueError("assembly primary reads exceed retained-read denominator")
    return {
        **mapping,
        "status": "available" if denominator is not None else "partial",
        "status_reason": None if denominator is not None else "denominator_unavailable",
        "assembly_mapped_reads": aligned,
        "bin_mapped_reads": binned,
        "unbinned_mapped_reads": aligned - binned,
        "without_qualifying_alignment_reads": denominator - aligned
        if denominator is not None
        else None,
        "per_contig": per_contig,
        "per_bin": list(by_bin.values()),
    }


def _mapping_policy(stages: dict[str, Row]) -> Row:
    command = str(stages.get("binning", {}).get("command", ""))
    view = re.search(r"samtools\s+view\s+([^;|\n]+)", command)
    flags = re.search(r"-F\s+(\d+)", view[1]) if view else None
    mapq = re.search(r"-q\s+(\d+)", view[1]) if view else None
    coverage = re.search(r"samtools\s+coverage\s+([^;|\n]+)", command)
    # Native samtools coverage defaults also exclude QC-failed and duplicate reads.
    mask = int(flags[1]) | 1796 if flags and coverage else None
    primary = mask is not None and mask & 2308 == 2308
    return {
        "minimum_mapq": int(mapq[1]) if mapq and primary else None,
        "excluded_sam_flags": mask if primary else None,
        "method": (
            "qualifying primary alignments from samtools coverage numreads; "
            "disjoint catalog membership"
        ),
        "mapping_command": command or None,
        "qc_commands": [
            {"stage_id": key, "command": stage["command"]}
            for key, stage in stages.items()
            if key == "read_qc"
        ],
        "coverage_policy": (
            "samtools coverage default excludes UNMAP,SECONDARY,QCFAIL,DUP; -d 0"
        )
        if coverage
        else None,
    }


def _coverage_counts(coverage: list[dict[str, str]], contigs: list[Row]) -> list[Row]:
    lengths = {row["contig_id"]: row["length_bp"] for row in contigs}
    rows = []
    seen = set()
    for row in coverage:
        identifier = row["#rname"]
        if identifier in seen or identifier not in lengths:
            raise ValueError("coverage must contain each assembly contig exactly once")
        seen.add(identifier)
        count, bases = int(row["numreads"]), int(row["covbases"])
        if count < 0 or not 0 <= bases <= lengths[identifier]:
            raise ValueError("invalid coverage read count or covered bases")
        if int(row["startpos"]) != 1 or int(row["endpos"]) != lengths[identifier]:
            raise ValueError("coverage interval must span the complete contig")
        rows.append({"contig_id": identifier, "numreads": count, "covbases": bases})
    if seen != set(lengths):
        raise ValueError("coverage does not cover exactly the assembly contigs")
    return rows


def _receipt(artifact: Row, sources: list[Row], status: str) -> Row:
    source_id = f"artifact:{artifact['artifact_id']}"
    for receipt in sources:
        if receipt["source_id"] == source_id:
            return receipt
    receipt = {**artifact, "source_id": source_id, "verification": status}
    sources.append(receipt)
    return receipt


def _read_artifact(artifact: Row, receipt: Row) -> str | None:
    path = Path(str(artifact["path"]))
    if not path.is_file():
        return None
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != artifact["sha256"]:
        raise ValueError(f"artifact SHA256 mismatch: {artifact['artifact_id']}")
    receipt["verification"] = "verified"
    return content.decode("utf-8")


def _parse_table(content: str) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(content), delimiter="\t"))


def _digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _encode(payload: Row) -> bytes:
    return (
        json.dumps(
            payload,
            sort_keys=True,
            ensure_ascii=False,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _write_json(path: Path, payload: Row) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_bytes(_encode(payload))
    temporary.replace(path)
