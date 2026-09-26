"""Validate workflow bundles and publish normalized SQLite catalogs."""

from __future__ import annotations

import sqlite3
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from pydantic import BaseModel, ValidationError

from protist_meta.models import (
    Artifact,
    MetadataBundle,
    ResultRecord,
    ScopedRecord,
    Stage,
)

__all__ = [
    "CatalogBuildError",
    "CatalogValidationError",
    "build_catalog",
    "validate_bundle",
]

type ScopedKey = tuple[str, str, str]
type RunKey = tuple[str, str]
type StageKey = tuple[str, str, str]

COLLECTION_ID_FIELDS = {
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
RESULT_ARTIFACT_KINDS = {
    "read_stats": "qc_report",
    "assemblies": "assembly",
    "contigs": "assembly",
    "bins": "bin_fasta",
    "memberships": "bin_fasta",
    "qc": "qc_report",
    "taxonomy": "taxonomy_report",
    "ssu": "ssu_report",
    "viruses": "viral_report",
    "phenotypes": "phenotype_report",
    "genes": "gene_annotation",
    "annotations": "functional_annotation",
}

TABLE_COLUMNS = {
    "samples": ("sample_id", "description", "genetic_code"),
    "runs": (
        "sample_id",
        "run_id",
        "platform",
        "rna_technology",
        "basecaller",
        "basecaller_model",
        "library_prep",
    ),
    "stages": (
        "sample_id",
        "run_id",
        "stage_id",
        "name",
        "status",
        "status_reason",
        "notes",
        "tool_name",
        "tool_version",
        "command",
        "database_name",
        "database_version",
        "database_sha256",
    ),
    "artifacts": (
        "sample_id",
        "run_id",
        "artifact_id",
        "kind",
        "path",
        "sha256",
        "size_bytes",
        "producer_stage_id",
    ),
    "read_stats": (
        "sample_id",
        "run_id",
        "record_id",
        "stage_id",
        "artifact_id",
        "total_reads",
        "total_bases",
        "retained_reads",
        "retained_bases",
        "rejected_reads",
        "rejected_bases",
    ),
    "assemblies": (
        "sample_id",
        "run_id",
        "assembly_id",
        "stage_id",
        "artifact_id",
        "assembly_length_bp",
        "contig_count",
        "n50_bp",
        "gc_fraction",
    ),
    "contigs": (
        "sample_id",
        "run_id",
        "contig_id",
        "stage_id",
        "artifact_id",
        "assembly_id",
        "length_bp",
        "gc_fraction",
        "mean_depth_x",
        "coverage_artifact_id",
        "is_circular",
    ),
    "bins": (
        "sample_id",
        "run_id",
        "bin_id",
        "stage_id",
        "artifact_id",
        "assembly_id",
        "candidate_class",
    ),
    "memberships": (
        "sample_id",
        "run_id",
        "membership_id",
        "stage_id",
        "artifact_id",
        "bin_id",
        "contig_id",
    ),
    "qc": (
        "sample_id",
        "run_id",
        "record_id",
        "stage_id",
        "artifact_id",
        "target_kind",
        "target_id",
        "completeness_percent",
        "contamination_percent",
        "completeness_basis",
        "contamination_basis",
        "completeness_query_support",
        "completeness_model_reliability",
        "contamination_query_support",
        "contamination_model_reliability",
        "checkeuk_pfam_ridge_completeness_percent",
        "checkeuk_bom_completeness_percent",
        "checkeuk_inventory74_present",
        "checkeuk_expected_markers_present",
        "checkeuk_expected_markers_total",
        "checkeuk_pfam1085_distinct",
    ),
    "taxonomy": (
        "sample_id",
        "run_id",
        "record_id",
        "stage_id",
        "artifact_id",
        "target_kind",
        "target_id",
        "rank",
        "taxon_name",
        "taxon_id",
        "score",
        "lineage_confidence",
    ),
    "ssu": (
        "sample_id",
        "run_id",
        "record_id",
        "stage_id",
        "artifact_id",
        "contig_id",
        "ssu_type",
        "start",
        "end",
        "strand",
        "hit_accession",
        "reference_identifiers",
        "reference_source",
        "reference_versions",
        "reference_taxonomy",
        "taxon_name",
        "taxonomy_assignment_method",
        "identity_percent",
        "query_coverage_percent",
    ),
    "viruses": (
        "sample_id",
        "run_id",
        "record_id",
        "stage_id",
        "artifact_id",
        "contig_id",
        "start",
        "end",
        "strand",
        "score",
        "viral_taxon",
    ),
    "phenotypes": (
        "sample_id",
        "run_id",
        "record_id",
        "stage_id",
        "artifact_id",
        "bin_id",
        "raw_class",
        "raw_score",
        "threshold",
        "passes_threshold",
        "thresholded_class",
        "is_applicable",
        "applicability_reason",
    ),
    "genes": (
        "sample_id",
        "run_id",
        "gene_id",
        "stage_id",
        "artifact_id",
        "contig_id",
        "bin_id",
        "start",
        "end",
        "strand",
        "phase",
        "genetic_code",
    ),
    "annotations": (
        "sample_id",
        "run_id",
        "annotation_id",
        "stage_id",
        "artifact_id",
        "gene_id",
        "database_name",
        "database_version",
        "accession",
        "function_name",
        "start",
        "end",
        "score",
        "evalue",
    ),
}

SCHEMA_SQL = """
CREATE TABLE catalog_metadata (
    schema_version TEXT PRIMARY KEY,
    workflow_name TEXT NOT NULL,
    workflow_version TEXT NOT NULL,
    source_revision TEXT NOT NULL
) STRICT;

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    description TEXT,
    genetic_code INTEGER CHECK (genetic_code >= 1)
) STRICT;

CREATE TABLE runs (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    platform TEXT NOT NULL CHECK (platform IN ('ont', 'pacbio_hifi')),
    rna_technology TEXT CHECK (
        rna_technology IN (
            'illumina_short', 'pacbio_isoseq', 'ont_cdna', 'ont_direct_rna'
        )
    ),
    basecaller TEXT,
    basecaller_model TEXT,
    library_prep TEXT,
    PRIMARY KEY (sample_id, run_id),
    FOREIGN KEY (sample_id) REFERENCES samples (sample_id),
    CHECK (platform = 'ont' OR (basecaller IS NULL AND basecaller_model IS NULL))
) STRICT;

CREATE TABLE stages (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    name TEXT NOT NULL,
    status TEXT NOT NULL CHECK (
        status IN ('completed', 'failed', 'skipped', 'pending')
    ),
    status_reason TEXT,
    notes TEXT,
    tool_name TEXT NOT NULL,
    tool_version TEXT NOT NULL,
    command TEXT,
    database_name TEXT,
    database_version TEXT,
    database_sha256 TEXT CHECK (
        database_sha256 IS NULL OR length(database_sha256) = 64
    ),
    PRIMARY KEY (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id) REFERENCES runs (sample_id, run_id),
    CHECK (
        (status = 'completed' AND status_reason IS NULL AND command IS NOT NULL)
        OR (status <> 'completed' AND status_reason IS NOT NULL)
    ),
    CHECK (
        (database_name IS NULL AND database_version IS NULL)
        OR (database_name IS NOT NULL AND database_version IS NOT NULL)
    )
) STRICT;

CREATE TABLE artifacts (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    kind TEXT NOT NULL CHECK (
        kind IN (
            'dna_reads', 'rna_reads', 'assembly', 'alignment', 'depth_table',
            'bin_fasta', 'qc_report', 'taxonomy_report', 'ssu_report',
            'viral_report', 'phenotype_report', 'gene_annotation',
            'protein_fasta', 'functional_annotation', 'workflow_report', 'other'
        )
    ),
    path TEXT NOT NULL CHECK (
        length(path) > 0
        AND instr(path, char(10)) = 0
        AND instr(path, char(13)) = 0
        AND instr(path, char(34)) = 0
        AND instr(path, char(39)) = 0
    ),
    sha256 TEXT NOT NULL CHECK (length(sha256) = 64),
    size_bytes INTEGER NOT NULL CHECK (size_bytes >= 0),
    producer_stage_id TEXT,
    PRIMARY KEY (sample_id, run_id, artifact_id),
    UNIQUE (sample_id, run_id, artifact_id, producer_stage_id),
    FOREIGN KEY (sample_id, run_id) REFERENCES runs (sample_id, run_id),
    FOREIGN KEY (sample_id, run_id, producer_stage_id)
        REFERENCES stages (sample_id, run_id, stage_id)
) STRICT;

CREATE TABLE stage_inputs (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    position INTEGER NOT NULL CHECK (position >= 0),
    artifact_id TEXT NOT NULL,
    PRIMARY KEY (sample_id, run_id, stage_id, position),
    UNIQUE (sample_id, run_id, stage_id, artifact_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id)
        REFERENCES artifacts (sample_id, run_id, artifact_id)
) STRICT;

CREATE TABLE stage_expected_keys (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    collection TEXT NOT NULL CHECK (
        collection IN (
            'read_stats', 'assemblies', 'contigs', 'bins', 'memberships',
            'qc', 'taxonomy', 'ssu', 'viruses', 'phenotypes', 'genes',
            'annotations'
        )
    ),
    record_id TEXT NOT NULL,
    PRIMARY KEY (sample_id, run_id, stage_id, collection, record_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id)
) STRICT;

CREATE TABLE read_stats (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    record_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    total_reads INTEGER NOT NULL CHECK (total_reads >= 0),
    total_bases INTEGER NOT NULL CHECK (total_bases >= 0),
    retained_reads INTEGER NOT NULL CHECK (retained_reads >= 0),
    retained_bases INTEGER NOT NULL CHECK (retained_bases >= 0),
    rejected_reads INTEGER NOT NULL CHECK (rejected_reads >= 0),
    rejected_bases INTEGER NOT NULL CHECK (rejected_bases >= 0),
    PRIMARY KEY (sample_id, run_id, record_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    CHECK (retained_reads + rejected_reads = total_reads),
    CHECK (retained_bases + rejected_bases = total_bases)
) STRICT;

CREATE TABLE assemblies (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    assembly_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    assembly_length_bp INTEGER NOT NULL CHECK (assembly_length_bp >= 0),
    contig_count INTEGER NOT NULL CHECK (contig_count >= 0),
    n50_bp INTEGER CHECK (n50_bp >= 1),
    gc_fraction REAL NOT NULL CHECK (gc_fraction BETWEEN 0 AND 1),
    PRIMARY KEY (sample_id, run_id, assembly_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    CHECK (
        (contig_count = 0 AND n50_bp IS NULL AND assembly_length_bp = 0)
        OR (contig_count > 0 AND n50_bp IS NOT NULL)
    )
) STRICT;

CREATE TABLE contigs (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    contig_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    assembly_id TEXT NOT NULL,
    length_bp INTEGER NOT NULL CHECK (length_bp >= 1),
    gc_fraction REAL NOT NULL CHECK (gc_fraction BETWEEN 0 AND 1),
    mean_depth_x REAL CHECK (mean_depth_x >= 0),
    coverage_artifact_id TEXT,
    is_circular INTEGER CHECK (is_circular IN (0, 1)),
    PRIMARY KEY (sample_id, run_id, contig_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, assembly_id)
        REFERENCES assemblies (sample_id, run_id, assembly_id),
    FOREIGN KEY (sample_id, run_id, coverage_artifact_id)
        REFERENCES artifacts (sample_id, run_id, artifact_id),
    CHECK (mean_depth_x IS NULL OR coverage_artifact_id IS NOT NULL)
) STRICT;

CREATE TABLE bins (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    bin_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    assembly_id TEXT NOT NULL,
    candidate_class TEXT NOT NULL CHECK (
        candidate_class IN (
            'eukaryotic_candidate', 'bacterial_archaeal_candidate',
            'viral_candidate', 'conflicting', 'unresolved', 'chaff'
        )
    ),
    PRIMARY KEY (sample_id, run_id, bin_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, assembly_id)
        REFERENCES assemblies (sample_id, run_id, assembly_id)
) STRICT;

CREATE TABLE memberships (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    membership_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    bin_id TEXT NOT NULL,
    contig_id TEXT NOT NULL,
    PRIMARY KEY (sample_id, run_id, membership_id),
    UNIQUE (sample_id, run_id, contig_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, bin_id)
        REFERENCES bins (sample_id, run_id, bin_id),
    FOREIGN KEY (sample_id, run_id, contig_id)
        REFERENCES contigs (sample_id, run_id, contig_id)
) STRICT;

CREATE TABLE qc (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    record_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    target_kind TEXT NOT NULL CHECK (
        target_kind IN ('assembly', 'contig', 'bin', 'virus', 'gene')
    ),
    target_id TEXT NOT NULL,
    completeness_percent REAL CHECK (completeness_percent BETWEEN 0 AND 100),
    contamination_percent REAL CHECK (
        contamination_percent BETWEEN 0 AND 1.7976931348623157e+308
    ),
    completeness_basis TEXT,
    contamination_basis TEXT,
    completeness_query_support TEXT,
    completeness_model_reliability TEXT,
    contamination_query_support TEXT,
    contamination_model_reliability TEXT,
    checkeuk_pfam_ridge_completeness_percent REAL CHECK (
        checkeuk_pfam_ridge_completeness_percent BETWEEN 0 AND 100
    ),
    checkeuk_bom_completeness_percent REAL CHECK (
        checkeuk_bom_completeness_percent BETWEEN 0 AND 100
    ),
    checkeuk_inventory74_present INTEGER CHECK (
        checkeuk_inventory74_present BETWEEN 0 AND 74
    ),
    checkeuk_expected_markers_present INTEGER CHECK (
        checkeuk_expected_markers_present BETWEEN 0 AND 74
    ),
    checkeuk_expected_markers_total INTEGER CHECK (
        checkeuk_expected_markers_total BETWEEN 1 AND 74
    ),
    checkeuk_pfam1085_distinct INTEGER CHECK (
        checkeuk_pfam1085_distinct BETWEEN 0 AND 1085
    ),
    PRIMARY KEY (sample_id, run_id, record_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    CHECK (
        (checkeuk_expected_markers_present IS NULL
            AND checkeuk_expected_markers_total IS NULL)
        OR (checkeuk_expected_markers_present IS NOT NULL
            AND checkeuk_expected_markers_total IS NOT NULL
            AND checkeuk_expected_markers_present <= checkeuk_expected_markers_total)
    )
) STRICT;

CREATE TABLE taxonomy (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    record_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    target_kind TEXT NOT NULL CHECK (
        target_kind IN ('assembly', 'contig', 'bin', 'virus', 'gene')
    ),
    target_id TEXT NOT NULL,
    rank TEXT NOT NULL,
    taxon_name TEXT NOT NULL,
    taxon_id TEXT,
    score REAL,
    lineage_confidence TEXT,
    PRIMARY KEY (sample_id, run_id, record_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        )
) STRICT;

CREATE TABLE ssu (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    record_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    contig_id TEXT NOT NULL,
    ssu_type TEXT NOT NULL CHECK (
        ssu_type IN (
            'rrna_16s', 'rrna_18s', 'rrna_23s', 'rrna_28s', 'rrna_5s',
            'rrna_5_8s'
        )
    ),
    start INTEGER NOT NULL CHECK (start >= 1),
    end INTEGER NOT NULL CHECK (end >= start),
    strand TEXT NOT NULL CHECK (strand IN ('plus', 'minus', 'unknown')),
    hit_accession TEXT,
    reference_identifiers TEXT,
    reference_source TEXT,
    reference_versions TEXT,
    reference_taxonomy TEXT,
    taxon_name TEXT,
    taxonomy_assignment_method TEXT,
    identity_percent REAL CHECK (identity_percent BETWEEN 0 AND 100),
    query_coverage_percent REAL CHECK (query_coverage_percent BETWEEN 0 AND 100),
    PRIMARY KEY (sample_id, run_id, record_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, contig_id)
        REFERENCES contigs (sample_id, run_id, contig_id)
) STRICT;

CREATE TABLE viruses (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    record_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    contig_id TEXT NOT NULL,
    start INTEGER NOT NULL CHECK (start >= 1),
    end INTEGER NOT NULL CHECK (end >= start),
    strand TEXT NOT NULL CHECK (strand IN ('plus', 'minus', 'unknown')),
    score REAL,
    viral_taxon TEXT,
    PRIMARY KEY (sample_id, run_id, record_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, contig_id)
        REFERENCES contigs (sample_id, run_id, contig_id)
) STRICT;

CREATE TABLE phenotypes (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    record_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    bin_id TEXT NOT NULL,
    raw_class TEXT NOT NULL,
    raw_score REAL NOT NULL CHECK (raw_score BETWEEN 0 AND 1),
    threshold REAL NOT NULL CHECK (threshold BETWEEN 0 AND 1),
    passes_threshold INTEGER NOT NULL CHECK (passes_threshold IN (0, 1)),
    thresholded_class TEXT,
    is_applicable INTEGER NOT NULL CHECK (is_applicable IN (0, 1)),
    applicability_reason TEXT,
    PRIMARY KEY (sample_id, run_id, record_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, bin_id)
        REFERENCES bins (sample_id, run_id, bin_id),
    CHECK (is_applicable = 1 OR applicability_reason IS NOT NULL)
) STRICT;

CREATE TABLE genes (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    gene_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    contig_id TEXT NOT NULL,
    bin_id TEXT,
    start INTEGER NOT NULL CHECK (start >= 1),
    end INTEGER NOT NULL CHECK (end >= start),
    strand TEXT NOT NULL CHECK (strand IN ('plus', 'minus', 'unknown')),
    phase INTEGER CHECK (phase BETWEEN 0 AND 2),
    genetic_code INTEGER NOT NULL CHECK (genetic_code >= 1),
    PRIMARY KEY (sample_id, run_id, gene_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, contig_id)
        REFERENCES contigs (sample_id, run_id, contig_id),
    FOREIGN KEY (sample_id, run_id, bin_id)
        REFERENCES bins (sample_id, run_id, bin_id)
) STRICT;

CREATE TABLE annotations (
    sample_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    annotation_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    gene_id TEXT NOT NULL,
    database_name TEXT NOT NULL,
    database_version TEXT NOT NULL,
    accession TEXT NOT NULL,
    function_name TEXT,
    start INTEGER CHECK (start >= 1),
    end INTEGER CHECK (end >= start),
    score REAL,
    evalue REAL CHECK (evalue >= 0),
    PRIMARY KEY (sample_id, run_id, annotation_id),
    FOREIGN KEY (sample_id, run_id, stage_id)
        REFERENCES stages (sample_id, run_id, stage_id),
    FOREIGN KEY (sample_id, run_id, artifact_id, stage_id)
        REFERENCES artifacts (
            sample_id, run_id, artifact_id, producer_stage_id
        ),
    FOREIGN KEY (sample_id, run_id, gene_id)
        REFERENCES genes (sample_id, run_id, gene_id),
    CHECK ((start IS NULL) = (end IS NULL))
) STRICT;
"""


class CatalogValidationError(ValueError):
    """Raised when a bundle fails shape or relationship validation."""


class CatalogBuildError(RuntimeError):
    """Raised when SQLite publication fails after bundle validation."""


@dataclass(frozen=True, slots=True)
class _ValidationIndex:
    sample_ids: set[str]
    run_keys: set[RunKey]
    stage_keys: set[ScopedKey]
    artifact_keys: set[ScopedKey]
    stages_by_key: dict[StageKey, Stage]
    artifacts_by_key: dict[ScopedKey, Artifact]


def validate_bundle(bundle_path: Path) -> MetadataBundle:
    """Load and validate a complete metadata bundle without writing outputs."""
    try:
        source = bundle_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise CatalogValidationError(
            f"cannot read metadata bundle: {bundle_path}"
        ) from error

    try:
        bundle = MetadataBundle.model_validate_json(source)
    except ValidationError as error:
        details = "; ".join(
            f"{_format_location(item['loc'])}: {item['type']}"
            for item in error.errors(include_input=False)
        )
        raise CatalogValidationError(
            f"bundle model validation failed: {details}"
        ) from error

    _validate_relationships(bundle)
    return bundle


def build_catalog(bundle_path: Path, output: Path) -> Path:
    """Build, check, and atomically publish a normalized SQLite catalog.

    Validation completes before the output directory is created. A failed build
    removes only its temporary database and preserves any catalog at ``output``.
    """
    bundle = validate_bundle(bundle_path)
    output_path = output.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary_path = Path(handle.name)

    try:
        _write_database(bundle, temporary_path)
        temporary_path.replace(output_path)
    except CatalogBuildError:
        temporary_path.unlink(missing_ok=True)
        raise
    except (OSError, sqlite3.Error) as error:
        temporary_path.unlink(missing_ok=True)
        raise CatalogBuildError(f"catalog publication failed: {output_path}") from error
    return output_path


def _format_location(location: tuple[int | str, ...]) -> str:
    return ".".join(str(part) for part in location) or "bundle"


def _text(value: object) -> str:
    if isinstance(value, Enum):
        return str(value.value)
    return str(value)


def _scope(record: ScopedRecord) -> RunKey:
    return (record.sample_id, record.run_id)


def _key(record: ScopedRecord, identifier: str) -> ScopedKey:
    value = getattr(record, identifier)
    if not isinstance(value, str):
        raise CatalogValidationError(f"{identifier} is not a string")
    return (record.sample_id, record.run_id, value)


def _stage_key(record: ResultRecord) -> StageKey:
    return (record.sample_id, record.run_id, record.stage_id)


def _require_unique(
    records: Iterable[ScopedRecord],
    identifier: str,
    label: str,
) -> set[ScopedKey]:
    keys: set[ScopedKey] = set()
    for record in records:
        key = _key(record, identifier)
        if key in keys:
            raise CatalogValidationError(f"duplicate {label} key: {key}")
        keys.add(key)
    return keys


def _require_reference(key: ScopedKey, valid_keys: set[ScopedKey], label: str) -> None:
    if key not in valid_keys:
        raise CatalogValidationError(f"missing {label} reference: {key}")


def _validate_relationships(bundle: MetadataBundle) -> None:
    sample_ids = _validate_samples_and_runs(bundle)
    run_keys = {(run.sample_id, run.run_id) for run in bundle.runs}
    stage_keys = _require_unique(bundle.stages, "stage_id", "stage")
    artifact_keys = _require_unique(bundle.artifacts, "artifact_id", "artifact")
    index = _ValidationIndex(
        sample_ids,
        run_keys,
        stage_keys,
        artifact_keys,
        {
            (stage.sample_id, stage.run_id, stage.stage_id): stage
            for stage in bundle.stages
        },
        {
            (artifact.sample_id, artifact.run_id, artifact.artifact_id): artifact
            for artifact in bundle.artifacts
        },
    )
    _validate_artifacts(bundle, index)
    _validate_stages(bundle, index)
    result_keys = _validate_result_links(bundle, index)
    _validate_stage_coverage(index.stages_by_key, result_keys)
    _validate_domain_relationships(bundle, index.stages_by_key)


def _validate_samples_and_runs(bundle: MetadataBundle) -> set[str]:
    sample_ids: set[str] = set()
    for sample in bundle.samples:
        if sample.sample_id in sample_ids:
            raise CatalogValidationError(f"duplicate sample key: {sample.sample_id}")
        sample_ids.add(sample.sample_id)

    run_keys: set[RunKey] = set()
    for run in bundle.runs:
        if run.sample_id not in sample_ids:
            raise CatalogValidationError(
                f"run references missing sample: {_scope(run)}"
            )
        run_key = _scope(run)
        if run_key in run_keys:
            raise CatalogValidationError(f"duplicate run key: {run_key}")
        run_keys.add(run_key)
        if _text(run.platform) != "ont" and (
            run.basecaller is not None or run.basecaller_model is not None
        ):
            raise CatalogValidationError(
                f"non-ONT run contains ONT basecaller provenance: {run_key}"
            )
    return sample_ids


def _validate_artifacts(
    bundle: MetadataBundle,
    index: _ValidationIndex,
) -> None:
    fingerprints: dict[str, tuple[str, int]] = {}
    for artifact in bundle.artifacts:
        if _scope(artifact) not in index.run_keys:
            raise CatalogValidationError(
                f"artifact references missing run: {_scope(artifact)}"
            )
        fingerprint = (artifact.sha256, artifact.size_bytes)
        if artifact.path in fingerprints and fingerprints[artifact.path] != fingerprint:
            raise CatalogValidationError(
                f"conflicting artifact content identity: {artifact.path}"
            )
        fingerprints[artifact.path] = fingerprint
        if artifact.producer_stage_id is not None:
            producer_key = (*_scope(artifact), artifact.producer_stage_id)
            _require_reference(producer_key, index.stage_keys, "producer stage")


def _validate_stages(bundle: MetadataBundle, index: _ValidationIndex) -> None:
    for stage in bundle.stages:
        if (
            stage.sample_id not in index.sample_ids
            or _scope(stage) not in index.run_keys
        ):
            raise CatalogValidationError(
                f"stage references missing run: {_scope(stage)}"
            )
        _validate_stage_metadata(stage)
        if len(stage.input_artifact_ids) != len(set(stage.input_artifact_ids)):
            raise CatalogValidationError(
                f"stage repeats an input artifact: {_key(stage, 'stage_id')}"
            )
        for artifact_id in stage.input_artifact_ids:
            artifact_key = (*_scope(stage), artifact_id)
            _require_reference(
                artifact_key,
                index.artifact_keys,
                "stage input artifact",
            )
            artifact = index.artifacts_by_key[artifact_key]
            if artifact.producer_stage_id == stage.stage_id:
                raise CatalogValidationError(
                    f"stage consumes its own output artifact: {_key(stage, 'stage_id')}"
                )


def _validate_stage_metadata(stage: Stage) -> None:
    identity = _key(stage, "stage_id")
    if (stage.database_name is None) != (stage.database_version is None):
        raise CatalogValidationError(
            f"stage database name/version must appear together: {identity}"
        )
    if stage.database_sha256 is not None and stage.database_name is None:
        raise CatalogValidationError(
            f"stage database checksum lacks database identity: {identity}"
        )
    status = _text(stage.status)
    if status == "completed" and stage.status_reason is not None:
        raise CatalogValidationError(f"completed stage has a status reason: {identity}")
    if status == "completed" and stage.command is None:
        raise CatalogValidationError(
            f"completed stage lacks its exact command: {identity}"
        )
    if status != "completed" and stage.status_reason is None:
        raise CatalogValidationError(
            f"non-completed stage lacks a status reason: {identity}"
        )


def _validate_result_links(
    bundle: MetadataBundle,
    index: _ValidationIndex,
) -> dict[StageKey, set[str]]:
    result_keys: dict[StageKey, set[str]] = defaultdict(set)
    for collection, identifier in COLLECTION_ID_FIELDS.items():
        records: Sequence[ResultRecord] = getattr(bundle, collection)
        _require_unique(records, identifier, collection)
        for record in records:
            if _scope(record) not in index.run_keys:
                raise CatalogValidationError(
                    f"{collection} record references missing run: {_scope(record)}"
                )
            stage_key = _stage_key(record)
            _require_reference(stage_key, index.stage_keys, "result stage")
            artifact_key = (*_scope(record), record.artifact_id)
            _require_reference(artifact_key, index.artifact_keys, "result artifact")
            coverage_artifact_id = getattr(record, "coverage_artifact_id", None)
            if coverage_artifact_id is not None:
                _require_reference(
                    (*_scope(record), coverage_artifact_id),
                    index.artifact_keys,
                    "coverage artifact",
                )
            stage = index.stages_by_key[stage_key]
            if _text(stage.status) != "completed":
                raise CatalogValidationError(
                    f"{collection} record belongs to non-completed stage: {stage_key}"
                )
            artifact = index.artifacts_by_key[artifact_key]
            if artifact.producer_stage_id != record.stage_id:
                raise CatalogValidationError(
                    f"{collection} artifact was not produced by its stage: "
                    f"{artifact_key}"
                )
            if _text(artifact.kind) != RESULT_ARTIFACT_KINDS[collection]:
                raise CatalogValidationError(
                    f"{collection} record references the wrong artifact kind: "
                    f"{artifact_key}"
                )
            record_id = getattr(record, identifier)
            result_keys[stage_key].add(f"{collection}:{record_id}")
    return result_keys


def _validate_stage_coverage(
    stages_by_key: dict[StageKey, object],
    result_keys: dict[StageKey, set[str]],
) -> None:
    for stage_key, stage in stages_by_key.items():
        expected = stage.expected_result_keys
        if len(expected) != len(set(expected)):
            raise CatalogValidationError(
                f"stage repeats an expected result key: {stage_key}"
            )
        if _text(stage.status) != "completed":
            if result_keys.get(stage_key):
                raise CatalogValidationError(
                    f"non-completed stage has result rows: {stage_key}"
                )
            continue
        actual = result_keys.get(stage_key, set())
        if set(expected) != actual:
            missing = sorted(set(expected) - actual)
            unexpected = sorted(actual - set(expected))
            raise CatalogValidationError(
                f"completed stage key coverage differs for {stage_key}; "
                f"missing={missing}, unexpected={unexpected}"
            )


def _validate_domain_relationships(
    bundle: MetadataBundle,
    stages_by_key: dict[StageKey, object],
) -> None:
    assembly_keys = {
        (record.sample_id, record.run_id, record.assembly_id)
        for record in bundle.assemblies
    }
    contigs_by_key = {
        (record.sample_id, record.run_id, record.contig_id): record
        for record in bundle.contigs
    }
    bin_keys = {
        (record.sample_id, record.run_id, record.bin_id) for record in bundle.bins
    }
    virus_keys = {
        (record.sample_id, record.run_id, record.record_id) for record in bundle.viruses
    }
    gene_keys = {
        (record.sample_id, record.run_id, record.gene_id) for record in bundle.genes
    }
    _validate_read_statistics(bundle)
    _validate_assemblies(bundle, assembly_keys)
    membership_pairs = _validate_bins(bundle, assembly_keys, contigs_by_key, bin_keys)
    target_keys = {
        "assembly": assembly_keys,
        "contig": set(contigs_by_key),
        "bin": bin_keys,
        "virus": virus_keys,
        "gene": gene_keys,
    }
    _validate_qc(bundle, stages_by_key, target_keys)
    _validate_taxonomy(bundle, target_keys)
    _validate_features(bundle, contigs_by_key, bin_keys, gene_keys, membership_pairs)


def _validate_read_statistics(bundle: MetadataBundle) -> None:
    for record in bundle.read_stats:
        if record.retained_reads + record.rejected_reads != record.total_reads:
            raise CatalogValidationError(
                f"retained and rejected read counts do not sum to total: "
                f"{_key(record, 'record_id')}"
            )
        if record.retained_bases + record.rejected_bases != record.total_bases:
            raise CatalogValidationError(
                f"retained and rejected base counts do not sum to total: "
                f"{_key(record, 'record_id')}"
            )


def _validate_assemblies(
    bundle: MetadataBundle,
    assembly_keys: set[ScopedKey],
) -> None:
    lengths_by_assembly: dict[ScopedKey, list[int]] = defaultdict(list)
    for contig in bundle.contigs:
        assembly_key = (*_scope(contig), contig.assembly_id)
        _require_reference(assembly_key, assembly_keys, "contig assembly")
        if contig.mean_depth_x is not None and contig.coverage_artifact_id is None:
            raise CatalogValidationError(
                f"contig depth lacks a coverage artifact: {_key(contig, 'contig_id')}"
            )
        lengths_by_assembly[assembly_key].append(contig.length_bp)

    for assembly in bundle.assemblies:
        assembly_key = (*_scope(assembly), assembly.assembly_id)
        lengths = lengths_by_assembly.get(assembly_key, [])
        if assembly.contig_count != len(lengths):
            raise CatalogValidationError(
                f"assembly contig count differs from contig rows: {assembly_key}"
            )
        if assembly.assembly_length_bp != sum(lengths):
            raise CatalogValidationError(
                f"assembly length differs from contig rows: {assembly_key}"
            )
        if lengths and assembly.n50_bp is None:
            raise CatalogValidationError(
                f"non-empty assembly lacks N50: {assembly_key}"
            )
        if not lengths and assembly.n50_bp is not None:
            raise CatalogValidationError(f"empty assembly has N50: {assembly_key}")
        if assembly.n50_bp is not None and assembly.n50_bp > max(lengths):
            raise CatalogValidationError(
                f"assembly N50 exceeds longest contig: {assembly_key}"
            )
        if lengths and assembly.n50_bp != _n50(lengths):
            raise CatalogValidationError(
                f"assembly N50 differs from contig rows: {assembly_key}"
            )


def _n50(lengths: Sequence[int]) -> int:
    half_length = (sum(lengths) + 1) // 2
    cumulative_length = 0
    for length in sorted(lengths, reverse=True):
        cumulative_length += length
        if cumulative_length >= half_length:
            return length
    raise ValueError("cannot calculate N50 from an empty sequence")


def _validate_bins(
    bundle: MetadataBundle,
    assembly_keys: set[ScopedKey],
    contigs_by_key: dict[ScopedKey, object],
    bin_keys: set[ScopedKey],
) -> set[tuple[str, str, str, str]]:
    for bin_record in bundle.bins:
        _require_reference(
            (*_scope(bin_record), bin_record.assembly_id),
            assembly_keys,
            "bin assembly",
        )

    membership_pairs: set[tuple[str, str, str, str]] = set()
    assigned_contigs: set[ScopedKey] = set()
    for membership in bundle.memberships:
        bin_key = (*_scope(membership), membership.bin_id)
        contig_key = (*_scope(membership), membership.contig_id)
        _require_reference(bin_key, bin_keys, "membership bin")
        if contig_key not in contigs_by_key:
            raise CatalogValidationError(
                f"missing membership contig reference: {contig_key}"
            )
        if contig_key in assigned_contigs:
            raise CatalogValidationError(
                f"contig has multiple bin memberships: {contig_key}"
            )
        assigned_contigs.add(contig_key)
        membership_pairs.add(
            (*_scope(membership), membership.bin_id, membership.contig_id)
        )
    return membership_pairs


def _validate_qc(
    bundle: MetadataBundle,
    stages_by_key: dict[StageKey, object],
    target_keys: dict[str, set[ScopedKey]],
) -> None:
    special_fields = (
        "checkeuk_pfam_ridge_completeness_percent",
        "checkeuk_bom_completeness_percent",
        "checkeuk_inventory74_present",
        "checkeuk_expected_markers_present",
        "checkeuk_expected_markers_total",
        "checkeuk_pfam1085_distinct",
    )
    for record in bundle.qc:
        _validate_target(record, target_keys)
        stage = stages_by_key[_stage_key(record)]
        is_checkeuk = _text(stage.tool_name).casefold() == "checkeuk"
        special_values = [getattr(record, field) for field in special_fields]
        if not is_checkeuk and any(value is not None for value in special_values):
            raise CatalogValidationError(
                f"non-CheckEUK QC contains CheckEUK fields: {_key(record, 'record_id')}"
            )
        if (
            record.completeness_percent is None
            and record.contamination_percent is None
            and all(value is None for value in special_values)
        ):
            raise CatalogValidationError(
                f"QC record contains no measurements: {_key(record, 'record_id')}"
            )
        if (
            record.completeness_percent is not None
            and record.completeness_basis is None
        ):
            raise CatalogValidationError(
                f"QC completeness lacks its basis: {_key(record, 'record_id')}"
            )
        if (
            record.contamination_percent is not None
            and record.contamination_basis is None
        ):
            raise CatalogValidationError(
                f"QC contamination lacks its basis: {_key(record, 'record_id')}"
            )
        present = record.checkeuk_expected_markers_present
        total = record.checkeuk_expected_markers_total
        if (present is None) != (total is None):
            raise CatalogValidationError(
                f"CheckEUK marker numerator/denominator must appear together: "
                f"{_key(record, 'record_id')}"
            )
        if present is not None and total is not None and present > total:
            raise CatalogValidationError(
                f"CheckEUK marker numerator exceeds denominator: "
                f"{_key(record, 'record_id')}"
            )


def _validate_taxonomy(
    bundle: MetadataBundle,
    target_keys: dict[str, set[ScopedKey]],
) -> None:
    for record in bundle.taxonomy:
        _validate_target(record, target_keys)


def _validate_target(
    record: ResultRecord,
    target_keys: dict[str, set[ScopedKey]],
) -> None:
    target_kind = _text(record.target_kind)
    target_id = record.target_id
    target_key = (*_scope(record), target_id)
    _require_reference(target_key, target_keys[target_kind], f"{target_kind} target")


def _validate_features(
    bundle: MetadataBundle,
    contigs_by_key: dict[ScopedKey, object],
    bin_keys: set[ScopedKey],
    gene_keys: set[ScopedKey],
    membership_pairs: set[tuple[str, str, str, str]],
) -> None:
    _validate_coordinates(bundle, contigs_by_key)
    _validate_genes(bundle, bin_keys, membership_pairs)
    _validate_phenotypes(bundle, bin_keys)
    _validate_annotations(bundle, gene_keys)


def _validate_coordinates(
    bundle: MetadataBundle,
    contigs_by_key: dict[ScopedKey, object],
) -> None:
    for record in [*bundle.ssu, *bundle.viruses, *bundle.genes]:
        contig_key = (*_scope(record), record.contig_id)
        if contig_key not in contigs_by_key:
            raise CatalogValidationError(f"missing feature contig: {contig_key}")
        if record.start > record.end:
            raise CatalogValidationError(f"feature start exceeds end: {contig_key}")
        contig = contigs_by_key[contig_key]
        if record.end > contig.length_bp:
            raise CatalogValidationError(f"feature exceeds contig length: {contig_key}")


def _validate_genes(
    bundle: MetadataBundle,
    bin_keys: set[ScopedKey],
    membership_pairs: set[tuple[str, str, str, str]],
) -> None:
    samples_by_id = {sample.sample_id: sample for sample in bundle.samples}
    prokaryotic_bins = {
        _key(bin_record, "bin_id")
        for bin_record in bundle.bins
        if bin_record.candidate_class == "bacterial_archaeal_candidate"
    }
    for gene in bundle.genes:
        if gene.bin_id is not None:
            bin_key = (*_scope(gene), gene.bin_id)
            _require_reference(bin_key, bin_keys, "gene bin")
            pair = (*_scope(gene), gene.bin_id, gene.contig_id)
            if pair not in membership_pairs:
                raise CatalogValidationError(
                    f"gene contig is not assigned to its bin: {_key(gene, 'gene_id')}"
                )
        sample_code = samples_by_id[gene.sample_id].genetic_code
        fixed_prokaryotic = (*_scope(gene), gene.bin_id) in prokaryotic_bins
        if (
            fixed_prokaryotic
            and sample_code is not None
            and gene.genetic_code != sample_code
        ):
            raise CatalogValidationError(
                f"gene code differs from prokaryotic override: {_key(gene, 'gene_id')}"
            )


def _validate_phenotypes(
    bundle: MetadataBundle,
    bin_keys: set[ScopedKey],
) -> None:
    for phenotype in bundle.phenotypes:
        _require_reference(
            (*_scope(phenotype), phenotype.bin_id),
            bin_keys,
            "phenotype bin",
        )
        if phenotype.passes_threshold != (phenotype.raw_score >= phenotype.threshold):
            raise CatalogValidationError(
                f"phenotype threshold flag disagrees with score: "
                f"{_key(phenotype, 'record_id')}"
            )
        if not phenotype.is_applicable and phenotype.applicability_reason is None:
            raise CatalogValidationError(
                f"inapplicable phenotype lacks a reason: {_key(phenotype, 'record_id')}"
            )


def _validate_annotations(
    bundle: MetadataBundle,
    gene_keys: set[ScopedKey],
) -> None:
    for annotation in bundle.annotations:
        _require_reference(
            (*_scope(annotation), annotation.gene_id),
            gene_keys,
            "annotation gene",
        )
        if (annotation.start is None) != (annotation.end is None):
            raise CatalogValidationError(
                f"annotation coordinates must appear together: "
                f"{_key(annotation, 'annotation_id')}"
            )
        if annotation.start is not None and annotation.start > annotation.end:
            raise CatalogValidationError(
                f"annotation start exceeds end: {_key(annotation, 'annotation_id')}"
            )


def _write_database(bundle: MetadataBundle, path: Path) -> None:
    with sqlite3.connect(path) as connection:
        connection.execute("PRAGMA foreign_keys = ON")
        connection.executescript(f"BEGIN IMMEDIATE;\n{SCHEMA_SQL}")
        _insert_bundle(connection, bundle)
        connection.commit()
        integrity = connection.execute("PRAGMA integrity_check").fetchone()
        if integrity != ("ok",):
            raise CatalogBuildError("SQLite integrity_check did not return ok")
        foreign_keys = connection.execute("PRAGMA foreign_key_check").fetchall()
        if foreign_keys:
            raise CatalogBuildError("SQLite foreign_key_check found violations")


def _insert_bundle(connection: sqlite3.Connection, bundle: MetadataBundle) -> None:
    connection.execute(
        "INSERT INTO catalog_metadata VALUES (?, ?, ?, ?)",
        (
            bundle.schema_version,
            bundle.workflow_name,
            bundle.workflow_version,
            bundle.source_revision,
        ),
    )
    _insert_models(connection, "samples", bundle.samples)
    _insert_models(connection, "runs", bundle.runs)
    _insert_models(connection, "stages", bundle.stages)
    _insert_models(connection, "artifacts", bundle.artifacts)
    _insert_stage_lists(connection, bundle)
    for table in COLLECTION_ID_FIELDS:
        _insert_models(connection, table, getattr(bundle, table))


def _insert_models(
    connection: sqlite3.Connection,
    table: str,
    models: Sequence[BaseModel],
) -> None:
    if not models:
        return
    columns = TABLE_COLUMNS[table]
    quoted_columns = ", ".join(f'"{column}"' for column in columns)
    placeholders = ", ".join("?" for _ in columns)
    statement = f'INSERT INTO "{table}" ({quoted_columns}) VALUES ({placeholders})'  # noqa: S608  # Identifiers come only from TABLE_COLUMNS; values are bound below.
    rows = []
    for model in models:
        values = model.model_dump(mode="json")
        rows.append(tuple(values.get(column) for column in columns))
    connection.executemany(statement, rows)


def _insert_stage_lists(
    connection: sqlite3.Connection,
    bundle: MetadataBundle,
) -> None:
    input_rows = []
    expected_rows = []
    for stage in bundle.stages:
        prefix = (stage.sample_id, stage.run_id, stage.stage_id)
        input_rows.extend(
            (*prefix, position, artifact_id)
            for position, artifact_id in enumerate(stage.input_artifact_ids)
        )
        expected_rows.extend(
            (*prefix, *result_key.split(":", maxsplit=1))
            for result_key in stage.expected_result_keys
        )
    connection.executemany(
        "INSERT INTO stage_inputs VALUES (?, ?, ?, ?, ?)",
        input_rows,
    )
    connection.executemany(
        "INSERT INTO stage_expected_keys VALUES (?, ?, ?, ?, ?)",
        expected_rows,
    )
