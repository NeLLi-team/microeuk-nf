"""Normalize core workflow outputs into the catalog's typed record bundle."""

import csv
import hashlib
import json
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from protist_meta.annotation_parsers import parse_eggnog, parse_interpro
from protist_meta.gene_inputs import lift_prodigal_gene, load_prodigal_origins
from protist_meta.gene_parsers import parse_braker3, parse_prodigal_gv
from protist_meta.quality_parsers import ParsedRows, parse_quality
from protist_meta.sequence_parsers import (
    ParsedSequenceRows,
    parse_checkv,
    parse_genomad,
    parse_ssu,
)
from protist_meta.sequences import read_fasta

__all__ = ["collect_records"]

COLLECTIONS = (
    "samples",
    "runs",
    "artifacts",
    "stages",
    "read_stats",
    "assemblies",
    "contigs",
    "bins",
    "memberships",
    "qc",
    "taxonomy",
    "ssu",
    "viruses",
    "phenotypes",
    "genes",
    "annotations",
)
FULL_STAGES = (
    "checkeuk",
    "checkm2",
    "checkm1",
    "gvclass",
    "gtdbtk",
    "ssuextract",
    "genomad",
    "checkv",
    "symclatron",
    "routing",
    "repeat_masking",
    "rna_alignment",
    "prodigal_gv",
    "braker3",
    "eggnog_mapper",
    "interproscan",
)
PATH_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
CHECKM2_NO_ANNOTATION_NOTE = (
    "CheckM2 ran gene calling and DIAMOND, but DIAMOND returned no annotations; "
    "CheckM2 quality prediction is unavailable."
)
RESULT_ID_FIELDS = {
    "read_stats": "record_id",
    "assemblies": "assembly_id",
    "contigs": "contig_id",
    "bins": "bin_id",
    "memberships": "membership_id",
    "qc": "record_id",
    "taxonomy": "record_id",
    "ssu": "record_id",
    "viruses": "record_id",
    "phenotypes": "record_id",
    "genes": "gene_id",
    "annotations": "annotation_id",
}
type Row = dict[str, object]
type Collections = dict[str, list[Row]]


class ToolManifest(BaseModel):
    """One scientific tool within a grouped Nextflow process."""

    model_config = ConfigDict(extra="forbid")

    version: str
    output: str
    status: Literal["completed", "failed", "skipped", "pending"] = "completed"
    database_name: str | None = None
    database_version: str | None = None
    reason: str | None = None


class StageManifest(BaseModel):
    """External stage metadata emitted by a completed Nextflow process."""

    model_config = ConfigDict(extra="forbid")

    stage: str
    status: Literal["completed", "failed", "skipped", "pending"]
    tool: str
    tool_version: str
    database_name: str | None = None
    database_version: str | None = None
    database_sha256: str | None = None
    command: str | list[str]
    outputs: dict[str, str] = Field(default_factory=dict)
    tools: dict[str, ToolManifest] = Field(default_factory=dict)
    reason: str | None = None


@dataclass(frozen=True, slots=True)
class CollectionContext:
    """Shared identity for records in one sample's workflow execution."""

    sample_id: str
    run_id: str
    records: Collections
    viruses_by_native: dict[str, str] = field(default_factory=dict)

    def scope(self) -> Row:
        """Return the composite scope required on run-linked records."""
        return {"sample_id": self.sample_id, "run_id": self.run_id}


def collect_records(sample_path: Path, stages: list[Path], output: Path) -> Path:
    """Parse every supplied stage; reject stages lacking an implemented parser."""
    sample = json.loads(sample_path.read_text(encoding="utf-8"))
    records: Collections = {name: [] for name in COLLECTIONS}
    context = CollectionContext(sample["sample_id"], sample["run_id"], records)
    records["samples"].append({"sample_id": context.sample_id})
    records["runs"].append({**context.scope(), "platform": sample["platform"].lower()})
    if sample.get("genetic_code", "auto") != "auto":
        records["samples"][0]["genetic_code"] = int(sample["genetic_code"])
    for name in ("basecaller", "basecaller_model", "library_prep"):
        if sample.get(name):
            records["runs"][0][name] = sample[name]
    raw_id = _artifact(context, Path(sample["reads"]), "raw_reads", "dna_reads", None)
    _collect_external_evidence(context, sample)
    for directory in sorted(stages):
        _collect_stage(context, directory, raw_id)
    present = {row["stage_id"] for row in records["stages"]}
    if sample.get("run_mode") == "core":
        _record_unrequested(context, present)
    bundle: Row = {
        "schema_version": "1.0.0",
        "workflow_name": "protist-meta-nf",
        "workflow_version": sample["workflow_version"],
        "source_revision": sample["source_revision"],
        **records,
    }
    output.write_text(json.dumps(bundle, indent=2) + "\n", encoding="utf-8")
    return output


def _collect_external_evidence(context: CollectionContext, sample: Row) -> None:
    for column, identifier, kind in (
        ("assembly", "supplied_assembly", "assembly"),
        ("rna_reads", "rna_reads", "rna_reads"),
        ("protein_reference", "protein_reference", "protein_fasta"),
        ("softmasked", "supplied_softmasked", "assembly"),
    ):
        if sample.get(column):
            _artifact(context, Path(str(sample[column])), identifier, kind, None)
    if sample.get("rna_platform"):
        technologies = {
            "ont_cdna": "ont_cdna",
            "ont_direct": "ont_direct_rna",
            "pacbio_isoseq": "pacbio_isoseq",
            "illumina": "illumina_short",
        }
        context.records["runs"][0]["rna_technology"] = technologies[
            str(sample["rna_platform"]).lower()
        ]


def _collect_stage(context: CollectionContext, directory: Path, raw_id: str) -> None:
    metadata = directory / "stage.json"
    manifest = StageManifest.model_validate_json(metadata.read_text(encoding="utf-8"))
    if manifest.stage == "braker3":
        _collect_manifest(context, directory, manifest, raw_id)
        if manifest.status in {"completed", "pending"}:
            _collect_braker_children(context, directory, manifest, raw_id)
        return
    if not manifest.tools:
        _collect_manifest(context, directory, manifest, raw_id)
        return
    if manifest.stage == "functional_annotation":
        _collect_manifest(context, directory, manifest, raw_id)
        if manifest.status == "skipped" and manifest.reason == "no_called_proteins":
            _collect_annotation_inputs(context, directory, manifest)
    for name, tool in manifest.tools.items():
        status = tool.status if manifest.status == "completed" else manifest.status
        reason = tool.reason or manifest.reason
        if status != "completed" and not reason:
            reason = f"{name} had no applicable input in {manifest.stage}."
        child = StageManifest(
            stage=name,
            status=status,
            tool=name,
            tool_version=tool.version,
            database_name=tool.database_name,
            database_version=tool.database_version,
            command=manifest.command,
            outputs={**manifest.outputs, "result": tool.output},
            reason=reason,
        )
        _collect_manifest(context, directory, child, raw_id)
        if (
            child.stage == "checkm2"
            and child.status == "skipped"
            and child.reason == "no_diamond_annotations"
        ):
            _collect_skipped_checkm2_evidence(context, directory, child)


def _collect_skipped_checkm2_evidence(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    """Retain native evidence for an eligible CheckM2 no-result run."""
    evidence = manifest.model_copy(
        update={"outputs": {"result": manifest.outputs["result"]}}
    )
    _collect_artifact_outputs(context, directory, evidence, {})
    context.records["stages"][-1]["notes"] = CHECKM2_NO_ANNOTATION_NOTE


def _collect_manifest(
    context: CollectionContext, directory: Path, manifest: StageManifest, raw_id: str
) -> None:
    metadata = directory / "stage.json"
    _artifact(context, metadata, f"{manifest.stage}.manifest", "other", manifest.stage)
    stage: Row = {
        **context.scope(),
        "stage_id": manifest.stage,
        "name": manifest.stage,
        "status": manifest.status,
        "tool_name": manifest.tool,
        "tool_version": manifest.tool_version.strip(),
        "database_name": manifest.database_name,
        "database_version": manifest.database_version,
        "database_sha256": manifest.database_sha256,
        "command": _command(directory, manifest),
        "input_artifact_ids": _stage_inputs(context, manifest.stage, raw_id),
        "expected_result_keys": [],
    }
    if manifest.reason:
        stage["status_reason"] = manifest.reason
    context.records["stages"].append(stage)
    if manifest.status != "completed":
        return
    _collect_completed_stage(context, directory, manifest)
    stage["expected_result_keys"] = [
        f"{collection}:{_record_id(collection, row)}"
        for collection, rows in context.records.items()
        if collection not in {"samples", "runs", "artifacts", "stages"}
        for row in rows
        if row["stage_id"] == manifest.stage
    ]


def _command(directory: Path, manifest: StageManifest) -> str:
    """Return the producer script when Nextflow provenance is available."""
    script = directory.resolve().parent / ".command.sh"
    if script.is_file():
        command = script.read_text(encoding="utf-8").strip()
        if not command:
            raise ValueError(f"producer command script is empty: {script}")
        return command
    if isinstance(manifest.command, list):
        return json.dumps(manifest.command)
    return manifest.command


def _collect_completed_stage(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    collectors: dict[str, Callable[[CollectionContext, Path, StageManifest], None]] = {
        "preflight": _collect_preflight,
        "read_qc": _collect_read_manifest,
        "qc": _collect_read_manifest,
        "assembly": _collect_assembly_manifest,
        "binning": _collect_binning_manifest,
        "checkeuk": _collect_quality,
        "checkm2": _collect_quality,
        "checkm1": _collect_quality,
        "gvclass": _collect_quality,
        "gtdbtk": _collect_quality,
        "symclatron": _collect_quality,
        "routing": _collect_routing_manifest,
        "repeat_masking": _collect_repeat_masking,
        "rna_alignment": _collect_rna_alignment,
        "prodigal_gv": _collect_prodigal_genes,
        "braker3": _collect_no_results,
        "functional_annotation": _collect_annotation_inputs,
        "eggnog_mapper": _collect_annotations,
        "interproscan": _collect_annotations,
        "ssu": _collect_sequences,
        "genomad": _collect_sequences,
        "checkv": _collect_sequences,
    }
    if manifest.stage.startswith("braker3."):
        _collect_braker_genes(context, directory, manifest)
        return
    try:
        collector = collectors[manifest.stage]
    except KeyError as error:
        message = f"no validated parser for stage {manifest.stage}"
        raise NotImplementedError(message) from error
    collector(context, directory, manifest)


def _collect_no_results(
    _context: CollectionContext,
    _directory: Path,
    _manifest: StageManifest,
) -> None:
    """Accept a completed stage with provenance artifacts but no result rows."""


def _collect_preflight(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    """Record declared preflight receipts without duplicating its manifest."""
    receipts = manifest.model_copy(
        update={
            "outputs": {
                label: path
                for label, path in manifest.outputs.items()
                if label != "report"
            }
        }
    )
    _collect_artifact_outputs(context, directory, receipts, {})


def _collect_read_manifest(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    _collect_reads(context, directory, manifest.stage)


def _collect_assembly_manifest(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    _collect_assembly(context, directory, manifest.stage)


def _collect_binning_manifest(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    _collect_binning(context, directory, manifest)


def _collect_routing_manifest(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    _collect_routing(context, directory, manifest.stage)


def _collect_repeat_masking(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    _collect_artifact_outputs(
        context,
        directory,
        manifest,
        {"masked_bins": "assembly"},
    )


def _collect_rna_alignment(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    _collect_artifact_outputs(
        context,
        directory,
        manifest,
        {
            "whole_assembly_bam": "alignment",
            "bin_bams": "alignment",
        },
    )


def _collect_quality(
    context: CollectionContext, directory: Path, manifest: StageManifest
) -> None:
    path = directory / manifest.outputs["result"]
    paths = sorted(path.glob("gtdbtk.*.summary.tsv")) if path.is_dir() else [path]
    if not paths:
        raise FileNotFoundError(f"no native quality tables at {path}")
    observed: list[str] = []
    gaps: set[str] = set()
    for report_path in paths:
        detail = None
        if manifest.stage == "checkeuk":
            detail = report_path.with_name("checkeuk_report_detail.tsv")
            _artifact(context, detail, "checkeuk.detail", "qc_report", manifest.stage)
        parsed = parse_quality(manifest.stage, report_path, detail_path=detail)
        observed.extend(parsed.target_ids)
        gaps.update(parsed.evidence_gaps)
        _append_quality(context, report_path, manifest.stage, parsed)
    expected = _quality_bins(context, manifest.stage)
    if len(observed) != len(set(observed)) or set(observed) != expected:
        raise ValueError(
            f"{manifest.stage} output does not cover exactly its input bins"
        )
    context.records["stages"][-1]["notes"] = "\n".join(sorted(gaps)) or None


def _append_quality(
    context: CollectionContext, path: Path, stage: str, parsed: ParsedRows
) -> None:
    tables = (
        ("qc", parsed.qc, "qc_report"),
        ("taxonomy", parsed.taxonomy, "taxonomy_report"),
        ("phenotypes", parsed.phenotypes, "phenotype_report"),
    )
    _append_results(context, path, stage, tables)


def _append_results(
    context: CollectionContext,
    path: Path,
    stage: str,
    tables: tuple[tuple[str, list[Row], str], ...],
) -> None:
    if not any(rows for _, rows, _ in tables):
        _artifact(context, path, f"{stage}.{path.stem}.empty", "other", stage)
    for collection, rows, artifact_kind in tables:
        if not rows:
            continue
        artifact_id = _artifact(
            context, path, f"{stage}.{path.stem}.{collection}", artifact_kind, stage
        )
        context.records[collection].extend(
            {
                **context.scope(),
                "stage_id": stage,
                "artifact_id": artifact_id,
                **row,
            }
            for row in rows
        )


def _collect_sequences(
    context: CollectionContext, directory: Path, manifest: StageManifest
) -> None:
    if manifest.stage == "ssu":
        path = directory / manifest.outputs["cmsearch_summary"]
        parsed = parse_ssu(path)
    elif manifest.stage == "genomad":
        path = directory / manifest.outputs["result"] / "input_virus_summary.tsv"
        lengths = {
            str(row["contig_id"]): int(row["length_bp"])
            for row in context.records["contigs"]
        }
        parsed = parse_genomad(path, lengths)
        viral_fasta = path.with_name("input_virus.fna")
        if viral_fasta.exists():
            _artifact(
                context, viral_fasta, "genomad.viral_fasta", "assembly", manifest.stage
            )
        context.viruses_by_native.update(
            zip(
                parsed.target_ids,
                (str(row["record_id"]) for row in parsed.viruses),
                strict=True,
            )
        )
    else:
        path = directory / manifest.outputs["result"] / "quality_summary.tsv"
        parsed = parse_checkv(path, context.viruses_by_native)
    _append_sequences(context, path, manifest.stage, parsed)
    context.records["stages"][-1]["notes"] = "\n".join(parsed.evidence_gaps) or None


def _append_sequences(
    context: CollectionContext, path: Path, stage: str, parsed: ParsedSequenceRows
) -> None:
    _append_results(
        context,
        path,
        stage,
        (
            ("ssu", parsed.ssu, "ssu_report"),
            ("viruses", parsed.viruses, "viral_report"),
            ("qc", parsed.qc, "qc_report"),
            ("taxonomy", parsed.taxonomy, "taxonomy_report"),
        ),
    )


def _quality_bins(context: CollectionContext, stage: str) -> set[str]:
    all_bins = stage in {"checkeuk", "checkm2", "gvclass"}
    return {
        str(row["bin_id"])
        for row in context.records["bins"]
        if all_bins or row["candidate_class"] == "bacterial_archaeal_candidate"
    }


def _collect_routing(context: CollectionContext, directory: Path, stage: str) -> None:
    path = directory / "evidence.tsv"
    _artifact(context, path, "routing.evidence", "taxonomy_report", stage)
    rows = _table(path)
    routed = {row["bin_id"]: row for row in rows}
    expected = {str(row["bin_id"]) for row in context.records["bins"]}
    if len(routed) != len(rows) or set(routed) != expected:
        raise ValueError("routing evidence does not cover exactly the input bins")
    for bin_record in context.records["bins"]:
        decision = routed[str(bin_record["bin_id"])]
        if bin_record["candidate_class"] != "chaff":
            bin_record["candidate_class"] = decision["candidate_class"]


def _collect_artifact_outputs(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
    kinds: dict[str, str],
) -> None:
    """Record every file below declared artifact-only stage outputs."""
    for label, relative in sorted(manifest.outputs.items()):
        output_path = directory / relative
        if not output_path.exists():
            raise FileNotFoundError(f"missing {manifest.stage} output: {output_path}")
        if output_path.is_dir():
            files = sorted(path for path in output_path.rglob("*") if path.is_file())
        else:
            files = [output_path]
            index = Path(f"{output_path}.bai")
            if index.is_file():
                files.append(index)
        for path in files:
            _artifact(
                context,
                path,
                _output_artifact_id(manifest.stage, label, path, directory),
                kinds.get(label, "other"),
                manifest.stage,
            )


def _collect_prodigal_genes(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    """Parse each included Prodigal-GV source and lift staged coordinates."""
    sources_path = directory / manifest.outputs["sources"]
    origins_path = directory / manifest.outputs["origins"]
    source_rows = _required_table(
        sources_path,
        ("source_id", "category", "bin_id", "fasta", "origin_map", "status", "reason"),
    )
    _artifact(context, sources_path, "prodigal_gv.sources", "other", manifest.stage)
    _artifact(context, origins_path, "prodigal_gv.origins", "other", manifest.stage)
    indexed: dict[str, dict[str, str]] = {}
    for row in source_rows:
        source_id = _path_identifier(row["source_id"], "Prodigal-GV source")
        if source_id in indexed:
            raise ValueError(f"duplicate Prodigal-GV source: {source_id}")
        if row["status"] not in {"included", "excluded"}:
            raise ValueError(f"invalid Prodigal-GV source status: {row['status']!r}")
        indexed[source_id] = row

    expected = {
        source_id for source_id, row in indexed.items() if row["status"] == "included"
    }
    output_roots = {
        "gff": directory / manifest.outputs["gff"],
        "proteins": directory / manifest.outputs["proteins"],
        "genes": directory / manifest.outputs["genes"],
    }
    for label, root in output_roots.items():
        suffix = {"gff": ".gff", "proteins": ".faa", "genes": ".fna"}[label]
        observed = {path.stem for path in root.glob(f"*{suffix}") if path.is_file()}
        if observed != expected:
            raise ValueError(
                f"Prodigal-GV {label}/source mismatch; "
                f"missing={sorted(expected - observed)}; "
                f"extra={sorted(observed - expected)}"
            )

    source_root = sources_path.parent
    sample_code = context.records["samples"][0].get("genetic_code")
    for source_id in sorted(expected):
        source = indexed[source_id]
        bin_id = source["bin_id"] or None
        if bin_id is not None:
            _path_identifier(bin_id, "Prodigal-GV bin")
        input_fasta = _relative_path(source_root, source["fasta"])
        origin_map = _relative_path(source_root, source["origin_map"])
        gff = output_roots["gff"] / f"{source_id}.gff"
        proteins = output_roots["proteins"] / f"{source_id}.faa"
        nucleotide = output_roots["genes"] / f"{source_id}.fna"
        prefix = f"prodigal_gv.{source_id}"
        _artifact(context, input_fasta, f"{prefix}.input", "assembly", manifest.stage)
        _artifact(context, origin_map, f"{prefix}.origins", "other", manifest.stage)
        artifact_id = _artifact(
            context, gff, f"{prefix}.gff", "gene_annotation", manifest.stage
        )
        _artifact(
            context, proteins, f"{prefix}.proteins", "protein_fasta", manifest.stage
        )
        _artifact(
            context,
            nucleotide,
            f"{prefix}.nucleotide",
            "gene_annotation",
            manifest.stage,
        )
        genetic_code = (
            int(sample_code)
            if source["category"] == "prokaryotic_bin" and sample_code is not None
            else None
        )
        genes = parse_prodigal_gv(
            gff,
            proteins,
            source_id=source_id,
            bin_id=bin_id,
            genetic_code=genetic_code,
        )
        origins = load_prodigal_origins(origin_map)
        context.records["genes"].extend(
            {
                **context.scope(),
                "stage_id": manifest.stage,
                "artifact_id": artifact_id,
                **lift_prodigal_gene(gene, origins),
            }
            for gene in genes
        )


def _collect_braker_children(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
    raw_id: str,
) -> None:
    """Retain a pending parent and completed per-bin BRAKER3 child stages."""
    expected = {
        str(row["bin_id"])
        for row in context.records["bins"]
        if row["candidate_class"] == "eukaryotic_candidate"
    }
    annotations = directory / manifest.outputs["annotations"]
    proteins = directory / manifest.outputs["proteins"]
    gff_bins = {
        path.parent.name
        for path in annotations.glob("*/braker.gff3")
        if path.is_file() and path.stat().st_size > 0
    }
    protein_bins = {
        path.stem
        for path in proteins.glob("*.faa")
        if path.is_file() and path.stat().st_size > 0
    }
    if gff_bins != protein_bins:
        raise ValueError(
            "BRAKER3 GFF/protein bin mismatch; "
            f"gff_only={sorted(gff_bins - protein_bins)}; "
            f"protein_only={sorted(protein_bins - gff_bins)}"
        )
    completed = gff_bins
    pending_path = directory / "pending.tsv"
    pending = _braker_pending(pending_path) if pending_path.exists() else {}
    if pending:
        _artifact(context, pending_path, "braker3.pending", "other", manifest.stage)
    overlap = completed.intersection(pending)
    if overlap:
        message = f"BRAKER3 bins are both completed and pending: {sorted(overlap)}"
        raise ValueError(message)
    if completed.union(pending) != expected:
        raise ValueError(
            "BRAKER3 completed/pending bins do not cover routed eukaryotic bins; "
            f"missing={sorted(expected - completed - set(pending))}; "
            f"extra={sorted((completed | set(pending)) - expected)}"
        )
    if manifest.status == "completed" and pending:
        raise ValueError("completed BRAKER3 parent contains pending bins")
    if manifest.status == "pending" and not pending:
        raise ValueError("pending BRAKER3 parent lacks pending.tsv entries")

    for bin_id in sorted(completed):
        _path_identifier(bin_id, "BRAKER3 bin")
        child = StageManifest(
            stage=f"braker3.{bin_id}",
            status="completed",
            tool=manifest.tool,
            tool_version=manifest.tool_version,
            database_name=manifest.database_name,
            database_version=manifest.database_version,
            database_sha256=manifest.database_sha256,
            command=manifest.command,
            outputs={
                "gff": str(
                    Path(manifest.outputs["annotations"]) / bin_id / "braker.gff3"
                ),
                "proteins": str(Path(manifest.outputs["proteins"]) / f"{bin_id}.faa"),
            },
        )
        _collect_manifest(context, directory, child, raw_id)


def _collect_braker_genes(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    """Parse one successfully called BRAKER3 bin into a child stage."""
    bin_id = manifest.stage.removeprefix("braker3.")
    gff = directory / manifest.outputs["gff"]
    proteins = directory / manifest.outputs["proteins"]
    artifact_id = _artifact(
        context, gff, f"{manifest.stage}.gff", "gene_annotation", manifest.stage
    )
    _artifact(
        context,
        proteins,
        f"{manifest.stage}.proteins",
        "protein_fasta",
        manifest.stage,
    )
    genes = parse_braker3(gff, proteins, source_id="braker3", bin_id=bin_id)
    context.records["genes"].extend(
        {
            **context.scope(),
            "stage_id": manifest.stage,
            "artifact_id": artifact_id,
            **gene,
        }
        for gene in genes
    )


def _collect_annotation_inputs(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    """Record the merged protein FASTA and its namespacing ledger."""
    _artifact(
        context,
        directory / manifest.outputs["proteins"],
        "functional_annotation.proteins",
        "protein_fasta",
        manifest.stage,
    )
    _artifact(
        context,
        directory / manifest.outputs["protein_map"],
        "functional_annotation.protein_map",
        "other",
        manifest.stage,
    )


def _collect_annotations(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    """Parse a grouped annotation tool using the shared protein ID map."""
    protein_ids = _annotation_protein_ids(
        context, directory / manifest.outputs["protein_map"]
    )
    database_version = manifest.database_version
    if database_version is None:
        raise ValueError(f"{manifest.stage} lacks a database version")
    result = directory / manifest.outputs["result"]
    if manifest.stage == "eggnog_mapper":
        path = result
        annotations = parse_eggnog(path, protein_ids, database_version)
    else:
        paths = sorted(path for path in result.glob("*.tsv") if path.is_file())
        if len(paths) != 1:
            raise ValueError(
                f"InterProScan output requires one TSV, observed {len(paths)}: {result}"
            )
        path = paths[0]
        annotations = parse_interpro(path, protein_ids, database_version)
    artifact_id = _artifact(
        context,
        path,
        f"{manifest.stage}.result",
        "functional_annotation",
        manifest.stage,
    )
    context.records["annotations"].extend(
        {
            **context.scope(),
            "stage_id": manifest.stage,
            "artifact_id": artifact_id,
            **annotation,
        }
        for annotation in annotations
    )


def _collect_reads(context: CollectionContext, directory: Path, stage: str) -> None:
    raw = _table(directory / "raw_stats.tsv")[0]
    retained = _table(directory / "filtered_stats.tsv")[0]
    artifact_id = _artifact(
        context, directory / "filtered_stats.tsv", "read_stats", "qc_report", stage
    )
    _artifact(
        context, directory / "raw_stats.tsv", "raw_read_stats", "qc_report", stage
    )
    _artifact(
        context, directory / "filtered.fastq.gz", "filtered_reads", "dna_reads", stage
    )
    before_reads, after_reads = int(raw["num_seqs"]), int(retained["num_seqs"])
    before_bases, after_bases = int(raw["sum_len"]), int(retained["sum_len"])
    context.records["read_stats"].append(
        {
            **context.scope(),
            "record_id": "reads",
            "stage_id": stage,
            "artifact_id": artifact_id,
            "total_reads": before_reads,
            "retained_reads": after_reads,
            "rejected_reads": before_reads - after_reads,
            "total_bases": before_bases,
            "retained_bases": after_bases,
            "rejected_bases": before_bases - after_bases,
        }
    )


def _collect_assembly(context: CollectionContext, directory: Path, stage: str) -> None:
    path = directory / "assembly_primary.fa"
    artifact_id = _artifact(context, path, "assembly", "assembly", stage)
    lengths: list[int] = []
    gc_count = 0
    for sequence in read_fasta(path):
        length = len(sequence.sequence)
        gc = sequence.sequence.upper().count("G") + sequence.sequence.upper().count("C")
        lengths.append(length)
        gc_count += gc
        context.records["contigs"].append(
            {
                **context.scope(),
                "contig_id": sequence.identifier,
                "assembly_id": "assembly",
                "stage_id": stage,
                "artifact_id": artifact_id,
                "length_bp": length,
                "gc_fraction": gc / length,
            }
        )
    if not lengths:
        raise ValueError("assembly contains no contigs")
    total = sum(lengths)
    context.records["assemblies"].append(
        {
            **context.scope(),
            "assembly_id": "assembly",
            "stage_id": stage,
            "artifact_id": artifact_id,
            "assembly_length_bp": total,
            "contig_count": len(lengths),
            "n50_bp": _n50(lengths),
            "gc_fraction": gc_count / total,
        }
    )


def _collect_binning(
    context: CollectionContext,
    directory: Path,
    manifest: StageManifest,
) -> None:
    stage = manifest.stage
    coverage_path = directory / manifest.outputs["coverage"]
    coverage = {row["#rname"]: float(row["meandepth"]) for row in _table(coverage_path)}
    coverage_id = _artifact(context, coverage_path, "coverage", "depth_table", stage)
    expected_contigs = {str(row["contig_id"]) for row in context.records["contigs"]}
    if set(coverage) != expected_contigs:
        raise ValueError("coverage does not cover exactly the assembly contigs")
    for contig in context.records["contigs"]:
        contig["mean_depth_x"] = coverage.get(str(contig["contig_id"]))
        contig["coverage_artifact_id"] = coverage_id
    for path in sorted((directory / manifest.outputs["bins"]).glob("*.fa")):
        bin_id = path.stem
        artifact_id = _artifact(context, path, f"bin.{bin_id}", "bin_fasta", stage)
        context.records["bins"].append(
            {
                **context.scope(),
                "bin_id": bin_id,
                "assembly_id": "assembly",
                "stage_id": stage,
                "artifact_id": artifact_id,
                "candidate_class": "unresolved",
            }
        )
        for sequence in read_fasta(path):
            context.records["memberships"].append(
                {
                    **context.scope(),
                    "membership_id": f"{bin_id}:{sequence.identifier}",
                    "bin_id": bin_id,
                    "contig_id": sequence.identifier,
                    "stage_id": stage,
                    "artifact_id": artifact_id,
                }
            )
    quickclade = manifest.outputs.get("quickclade")
    if quickclade is None:
        if context.records["bins"]:
            raise ValueError("binning omitted QuickClade output for accepted bins")
        return
    _collect_quickclade(context, directory / quickclade, stage)


def _collect_quickclade(context: CollectionContext, path: Path, stage: str) -> None:
    artifact_id = _artifact(context, path, "quickclade", "taxonomy_report", stage)
    rows = _table(path)
    observed = [Path(row["#QueryName"]).stem for row in rows]
    expected = {str(row["bin_id"]) for row in context.records["bins"]}
    if len(observed) != len(set(observed)) or set(observed) != expected:
        raise ValueError("QuickClade output does not cover exactly the accepted bins")
    for row, bin_id in zip(rows, observed, strict=True):
        context.records["taxonomy"].append(
            {
                **context.scope(),
                "record_id": f"quickclade:{bin_id}",
                "stage_id": stage,
                "artifact_id": artifact_id,
                "target_kind": "bin",
                "target_id": bin_id,
                "rank": "lineage",
                "taxon_name": row["lineage"],
                "taxon_id": row["R_TaxID"],
            }
        )


def _artifact(
    context: CollectionContext,
    path: Path,
    identifier: str,
    kind: str,
    stage: str | None,
) -> str:
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    context.records["artifacts"].append(
        {
            **context.scope(),
            "artifact_id": identifier,
            "kind": kind,
            "path": str(path.resolve()),
            "sha256": digest,
            "size_bytes": path.stat().st_size,
            "producer_stage_id": stage,
        }
    )
    return identifier


def _record_id(collection: str, row: Row) -> object:
    identifier = RESULT_ID_FIELDS[collection]
    try:
        return row[identifier]
    except KeyError as error:
        raise ValueError(f"{collection} result lacks {identifier}") from error


def _table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _required_table(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
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


def _braker_pending(path: Path) -> dict[str, str]:
    pending: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not all(fields):
                raise ValueError(f"malformed BRAKER3 pending row {line_number}: {path}")
            bin_id = _path_identifier(fields[0], "BRAKER3 pending bin")
            if bin_id in pending:
                raise ValueError(f"duplicate BRAKER3 pending bin: {bin_id}")
            pending[bin_id] = fields[1]
    return pending


def _annotation_protein_ids(
    context: CollectionContext,
    path: Path,
) -> set[str]:
    rows = _required_table(
        path,
        (
            "normalized_id",
            "native_id",
            "native_gene_id",
            "caller",
            "source_id",
            "bin_id",
        ),
    )
    protein_ids = {row["normalized_id"] for row in rows}
    if len(protein_ids) != len(rows):
        raise ValueError(f"protein map repeats a normalized ID: {path}")
    gene_ids = {str(row["gene_id"]) for row in context.records["genes"]}
    if protein_ids != gene_ids:
        raise ValueError(
            "protein map does not cover exactly the collected genes; "
            f"missing={sorted(gene_ids - protein_ids)}; "
            f"extra={sorted(protein_ids - gene_ids)}"
        )
    return protein_ids


def _path_identifier(value: str, label: str) -> str:
    if PATH_IDENTIFIER.fullmatch(value) is None:
        raise ValueError(f"invalid {label} identifier: {value!r}")
    return value


def _relative_path(root: Path, relative: str) -> Path:
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"output path escapes its stage directory: {relative}")
    return path


def _output_artifact_id(stage: str, label: str, path: Path, root: Path) -> str:
    relative = path.relative_to(root).as_posix()
    digest = hashlib.sha256(relative.encode()).hexdigest()[:12]
    return f"{stage}.{label}.{digest}"


def _n50(lengths: list[int]) -> int:
    midpoint = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= midpoint:
            return length
    raise ValueError("cannot calculate N50 from an empty assembly")


def _record_unrequested(context: CollectionContext, present: set[object]) -> None:
    for name in FULL_STAGES:
        if name not in present:
            context.records["stages"].append(
                {
                    **context.scope(),
                    "stage_id": name,
                    "name": name,
                    "status": "skipped",
                    "status_reason": "Analysis excluded by core validation mode.",
                    "tool_name": name,
                    "tool_version": "not_run",
                    "input_artifact_ids": [],
                    "expected_result_keys": [],
                }
            )


def _stage_inputs(context: CollectionContext, stage: str, raw_id: str) -> list[str]:
    if stage.startswith("braker3."):
        return _braker_inputs(context, stage.removeprefix("braker3."))
    builders: dict[str, Callable[[CollectionContext, str], list[str]]] = {
        "assembly": _assembly_inputs,
        "checkeuk": _quality_inputs,
        "checkm2": _quality_inputs,
        "checkm1": _quality_inputs,
        "gvclass": _quality_inputs,
        "gtdbtk": _quality_inputs,
        "symclatron": _quality_inputs,
        "routing": _routing_inputs,
        "repeat_masking": _repeat_masking_inputs,
        "rna_alignment": _rna_alignment_inputs,
        "prodigal_gv": _prodigal_inputs,
        "braker3": _braker_parent_inputs,
        "functional_annotation": _functional_annotation_inputs,
        "eggnog_mapper": _annotation_inputs,
        "interproscan": _annotation_inputs,
        "checkv": _checkv_inputs,
    }
    if builder := builders.get(stage):
        return builder(context, stage)
    return {
        "read_qc": [raw_id],
        "qc": [raw_id],
        "binning": ["filtered_reads", "assembly"],
        "ssu": ["assembly"],
        "genomad": ["assembly"],
    }.get(stage, [])


def _assembly_inputs(context: CollectionContext, _stage: str) -> list[str]:
    supplied = any(
        row["artifact_id"] == "supplied_assembly"
        for row in context.records["artifacts"]
    )
    return ["supplied_assembly" if supplied else "filtered_reads"]


def _quality_inputs(context: CollectionContext, stage: str) -> list[str]:
    return [f"bin.{name}" for name in sorted(_quality_bins(context, stage))]


def _routing_inputs(context: CollectionContext, _stage: str) -> list[str]:
    return _artifact_ids(
        context,
        lambda row: (
            row["kind"]
            in {"bin_fasta", "taxonomy_report", "ssu_report", "viral_report"}
        ),
    )


def _repeat_masking_inputs(context: CollectionContext, _stage: str) -> list[str]:
    eukaryotic_bins = {
        str(row["bin_id"])
        for row in context.records["bins"]
        if row["candidate_class"] == "eukaryotic_candidate"
    }
    return _artifact_ids(
        context,
        lambda row: (
            row["artifact_id"] == "supplied_softmasked"
            or (
                row["kind"] == "bin_fasta"
                and Path(str(row["path"])).stem in eukaryotic_bins
            )
        ),
    )


def _rna_alignment_inputs(context: CollectionContext, _stage: str) -> list[str]:
    return _artifact_ids(
        context,
        lambda row: (
            row["artifact_id"] in {"assembly", "rna_reads"}
            or (
                row["producer_stage_id"] == "repeat_masking"
                and row["kind"] == "assembly"
            )
        ),
    )


def _prodigal_inputs(context: CollectionContext, _stage: str) -> list[str]:
    selected_bins = {
        str(row["bin_id"])
        for row in context.records["bins"]
        if row["candidate_class"] in {"bacterial_archaeal_candidate", "viral_candidate"}
    }
    return _artifact_ids(
        context,
        lambda row: (
            row["artifact_id"] in {"routing.evidence", "genomad.viral_fasta"}
            or (
                row["kind"] == "bin_fasta"
                and Path(str(row["path"])).stem in selected_bins
            )
        ),
    )


def _braker_parent_inputs(context: CollectionContext, _stage: str) -> list[str]:
    return _braker_inputs(context, None)


def _functional_annotation_inputs(
    context: CollectionContext,
    _stage: str,
) -> list[str]:
    return _artifact_ids(context, _is_gene_protein)


def _annotation_inputs(context: CollectionContext, _stage: str) -> list[str]:
    return _artifact_ids(
        context,
        lambda row: (
            row["artifact_id"]
            in {
                "functional_annotation.proteins",
                "functional_annotation.protein_map",
            }
        ),
    )


def _checkv_inputs(context: CollectionContext, _stage: str) -> list[str]:
    return _artifact_ids(
        context,
        lambda row: row["artifact_id"] == "genomad.viral_fasta",
    )


def _artifact_ids(
    context: CollectionContext,
    predicate: Callable[[Row], bool],
) -> list[str]:
    selected = {
        str(row["artifact_id"])
        for row in context.records["artifacts"]
        if predicate(row)
    }
    return sorted(selected)


def _is_gene_protein(row: Row) -> bool:
    return (
        row["kind"] == "protein_fasta"
        and row["producer_stage_id"] is not None
        and (
            row["producer_stage_id"] == "prodigal_gv"
            or str(row["producer_stage_id"]).startswith("braker3.")
        )
    )


def _braker_inputs(context: CollectionContext, bin_id: str | None) -> list[str]:
    def selected(row: Row) -> bool:
        if row["artifact_id"] == "protein_reference":
            return True
        producer = row["producer_stage_id"]
        if producer not in {"repeat_masking", "rna_alignment"}:
            return False
        if bin_id is None:
            return True
        name = Path(str(row["path"])).name
        return name.startswith(f"{bin_id}.")

    return _artifact_ids(context, selected)
