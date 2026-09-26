"""Pydantic models generated from the authoritative LinkML catalog schema."""

from __future__ import annotations

from enum import StrEnum
from typing import Annotated

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
)


class ConfiguredBaseModel(BaseModel):
    """Reject unknown fields and implicit type conversions."""

    model_config = ConfigDict(
        strict=True,
        extra="forbid",
        validate_assignment=True,
        validate_default=True,
        use_enum_values=True,
    )


class Platform(StrEnum):
    """Allowed Platform values."""

    ont = "ont"
    pacbio_hifi = "pacbio_hifi"


class RnaTechnology(StrEnum):
    """Allowed RnaTechnology values."""

    illumina_short = "illumina_short"
    pacbio_isoseq = "pacbio_isoseq"
    ont_cdna = "ont_cdna"
    ont_direct_rna = "ont_direct_rna"


class ArtifactKind(StrEnum):
    """Allowed ArtifactKind values."""

    dna_reads = "dna_reads"
    rna_reads = "rna_reads"
    assembly = "assembly"
    alignment = "alignment"
    depth_table = "depth_table"
    bin_fasta = "bin_fasta"
    qc_report = "qc_report"
    taxonomy_report = "taxonomy_report"
    ssu_report = "ssu_report"
    viral_report = "viral_report"
    phenotype_report = "phenotype_report"
    gene_annotation = "gene_annotation"
    protein_fasta = "protein_fasta"
    functional_annotation = "functional_annotation"
    workflow_report = "workflow_report"
    other = "other"


class StageStatus(StrEnum):
    """Allowed StageStatus values."""

    completed = "completed"
    failed = "failed"
    skipped = "skipped"
    pending = "pending"


class TargetKind(StrEnum):
    """Allowed TargetKind values."""

    assembly = "assembly"
    contig = "contig"
    bin = "bin"
    virus = "virus"
    gene = "gene"


class CandidateClass(StrEnum):
    """Allowed CandidateClass values."""

    eukaryotic_candidate = "eukaryotic_candidate"
    bacterial_archaeal_candidate = "bacterial_archaeal_candidate"
    viral_candidate = "viral_candidate"
    conflicting = "conflicting"
    unresolved = "unresolved"
    chaff = "chaff"


class Strand(StrEnum):
    """Allowed Strand values."""

    plus = "plus"
    minus = "minus"
    unknown = "unknown"


class SsuType(StrEnum):
    """Allowed SsuType values."""

    rrna_16s = "rrna_16s"
    rrna_18s = "rrna_18s"
    rrna_23s = "rrna_23s"
    rrna_28s = "rrna_28s"
    rrna_5s = "rrna_5s"
    rrna_5_8s = "rrna_5_8s"


class ScopedRecord(ConfiguredBaseModel):
    """Validated ScopedRecord record."""

    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class ResultRecord(ScopedRecord):
    """Validated ResultRecord record."""

    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Sample(ConfiguredBaseModel):
    """Validated Sample record."""

    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    description: str | None = Field(
        default=None,
        description="""Short operator-supplied description.""",
    )
    genetic_code: int | None = Field(
        default=None,
        description="""Optional translation-table constraint for routed
    prokaryotic bins only. Eukaryotic and viral calls retain
    their own translation tables.""",
        ge=1,
    )


class Run(ScopedRecord):
    """Validated Run record."""

    platform: Platform = Field(
        default=...,
        description="""DNA sequencing platform for this run.""",
    )
    rna_technology: RnaTechnology | None = Field(
        default=None,
        description="""Technology of optional, sample-bound RNA evidence.""",
    )
    basecaller: str | None = Field(
        default=None,
        description="""Basecaller name recorded for basecalled ONT reads, when
    available.""",
        pattern=r"^\S+$",
    )
    basecaller_model: str | None = Field(
        default=None,
        description="""Exact basecaller model recorded for basecalled ONT
    reads.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    library_prep: str | None = Field(
        default=None,
        description="""Sequencing library preparation recorded for this run.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Artifact(ScopedRecord):
    """Validated Artifact record."""

    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    kind: ArtifactKind = Field(
        default=...,
        description="""Artifact role.""",
    )
    path: str = Field(
        default=...,
        description="""Workflow-relative or absolute artifact path recorded at
    publication.""",
        pattern=r"^[^\x22\x27\r\n]+$",
    )
    sha256: str = Field(
        default=...,
        description="""Lowercase SHA-256 digest of artifact bytes.""",
        pattern=r"^[0-9a-f]{64}$",
    )
    size_bytes: int = Field(
        default=...,
        description="""Artifact size in bytes.""",
        ge=0,
    )
    producer_stage_id: str | None = Field(
        default=None,
        description="""Stage that produced the artifact; absent only for
    external inputs.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Stage(ScopedRecord):
    """Validated Stage record."""

    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    name: str = Field(
        default=...,
        description="""Stable human-readable name.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    status: StageStatus = Field(
        default=...,
        description="""Terminal workflow stage state.""",
    )
    status_reason: str | None = Field(
        default=None,
        description="""Required reason for failed, skipped, or pending stages.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    notes: str | None = Field(
        default=None,
        description="""Interpretation limits and native evidence retained
    outside normalized fields.""",
    )
    tool_name: str = Field(
        default=...,
        description="""Executable or workflow component that ran the stage.""",
        pattern=r"^\S+$",
    )
    tool_version: str = Field(
        default=...,
        description="""Exact tool version or immutable source revision.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    command: str | None = Field(
        default=None,
        description="""Exact executed command for a completed stage, with
    secrets excluded.""",
        pattern=r"^\S(?:[\s\S]*\S)?$",
    )
    database_name: str | None = Field(
        default=None,
        description="""Reference database name used by the stage, when
    applicable.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    database_version: str | None = Field(
        default=None,
        description="""Exact reference database release or immutable digest.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    database_sha256: str | None = Field(
        default=None,
        description="""SHA-256 digest of the database manifest, when recorded.""",
        pattern=r"^[0-9a-f]{64}$",
    )
    input_artifact_ids: list[str] = Field(
        default=...,
        description="""Exact input artifacts consumed by this stage in the same
    sample and run.""",
    )
    expected_result_keys: list[
        Annotated[
            str,
            Field(
                pattern=r"^(read_stats|assemblies|contigs|bins|memberships|qc|taxonomy|ssu|viruses|phenotypes|genes|annotations):[A-Za-z0-9][A-Za-z0-9_.:-]*$"
            ),
        ]
    ] = Field(
        default=...,
        description="""Exact collection-qualified result keys expected when
    this stage is completed, for example assemblies:asm1 or
    contigs:ctg1.""",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class ReadStatistics(ResultRecord):
    """Validated ReadStatistics record."""

    record_id: str = Field(
        default=...,
        description="""Stable record identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    total_reads: int = Field(
        default=...,
        description="""Number of reads presented to filtering or summarization.""",
        ge=0,
    )
    total_bases: int = Field(
        default=...,
        description="""Total read bases presented, in base pairs.""",
        ge=0,
    )
    retained_reads: int = Field(
        default=...,
        description="""Number of reads retained after filtering.""",
        ge=0,
    )
    retained_bases: int = Field(
        default=...,
        description="""Number of retained read bases, in base pairs.""",
        ge=0,
    )
    rejected_reads: int = Field(
        default=...,
        description="""Number of reads rejected by filtering.""",
        ge=0,
    )
    rejected_bases: int = Field(
        default=...,
        description="""Number of rejected read bases, in base pairs.""",
        ge=0,
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Assembly(ResultRecord):
    """Validated Assembly record."""

    assembly_id: str = Field(
        default=...,
        description="""Stable assembly identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    assembly_length_bp: int = Field(
        default=...,
        description="""Total assembled sequence length in base pairs.""",
        ge=0,
    )
    contig_count: int = Field(
        default=...,
        description="""Number of contigs in the assembly.""",
        ge=0,
    )
    n50_bp: int | None = Field(
        default=None,
        description="""Assembly N50 contig length in base pairs; absent for an
    empty assembly.""",
        ge=1,
    )
    gc_fraction: float = Field(
        default=...,
        description="""Fraction of G and C bases on a zero-to-one scale.""",
        ge=0,
        le=1,
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Contig(ResultRecord):
    """Validated Contig record."""

    contig_id: str = Field(
        default=...,
        description="""Stable contig identifier scoped by sample_id and run_id.""",
        pattern=r"^\S+$",
    )
    assembly_id: str = Field(
        default=...,
        description="""Stable assembly identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    length_bp: int = Field(
        default=...,
        description="""Sequence length in base pairs.""",
        ge=1,
    )
    gc_fraction: float = Field(
        default=...,
        description="""Fraction of G and C bases on a zero-to-one scale.""",
        ge=0,
        le=1,
    )
    mean_depth_x: float | None = Field(
        default=None,
        description="""Mean mapped DNA read depth in fold coverage.""",
        ge=0,
    )
    coverage_artifact_id: str | None = Field(
        default=None,
        description="""Artifact containing the mapping or depth evidence for
    mean_depth_x.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    is_circular: bool | None = Field(
        default=None,
        description="""Whether the producing tool called the sequence circular.""",
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Bin(ResultRecord):
    """Validated Bin record."""

    bin_id: str = Field(
        default=...,
        description="""Stable bin identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    assembly_id: str = Field(
        default=...,
        description="""Stable assembly identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    candidate_class: CandidateClass = Field(
        default=...,
        description="""Evidence-routing class; this field is not a MAG
    acceptance decision.""",
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Membership(ResultRecord):
    """Validated Membership record."""

    membership_id: str = Field(
        default=...,
        description="""Stable bin-membership identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    bin_id: str = Field(
        default=...,
        description="""Stable bin identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    contig_id: str = Field(
        default=...,
        description="""Stable contig identifier scoped by sample_id and run_id.""",
        pattern=r"^\S+$",
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class QualityAssessment(ResultRecord):
    """Validated QualityAssessment record."""

    record_id: str = Field(
        default=...,
        description="""Stable record identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    target_kind: TargetKind = Field(
        default=...,
        description="""Entity type evaluated by this record.""",
    )
    target_id: str = Field(
        default=...,
        description="""Identifier of the evaluated entity in the same sample
    and run.""",
        pattern=r"^\S+$",
    )
    completeness_percent: float | None = Field(
        default=None,
        description="""Tool-reported generic completeness in percent, with
    tool-specific meaning.""",
        ge=0,
        le=100,
    )
    contamination_percent: float | None = Field(
        default=None,
        description="""Tool-reported contamination in percent, with tool-
    specific meaning.""",
        ge=0,
        le=1.7976931348623157e308,
    )
    completeness_basis: str | None = Field(
        default=None,
        description="""Exact method or model defining completeness_percent.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    contamination_basis: str | None = Field(
        default=None,
        description="""Exact method or denominator defining
    contamination_percent.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    completeness_query_support: str | None = Field(
        default=None,
        description="""Verbatim GVClass query-support label for the
    completeness estimate, not a calibrated probability.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    completeness_model_reliability: str | None = Field(
        default=None,
        description="""Verbatim GVClass model-reliability label for the
    completeness estimate, not a calibrated probability.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    contamination_query_support: str | None = Field(
        default=None,
        description="""Verbatim GVClass query-support label for the
    contamination estimate, not a calibrated probability.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    contamination_model_reliability: str | None = Field(
        default=None,
        description="""Verbatim GVClass model-reliability label for the
    contamination estimate, not a calibrated probability.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    checkeuk_pfam_ridge_completeness_percent: float | None = Field(
        default=None,
        description="""CheckEUK Pfam-union ridge completeness in percent.""",
        ge=0,
        le=100,
    )
    checkeuk_bom_completeness_percent: float | None = Field(
        default=None,
        description="""CheckEUK best-of-models completeness estimate in
    percent.""",
        ge=0,
        le=100,
    )
    checkeuk_inventory74_present: int | None = Field(
        default=None,
        description="""Count of present families in the fixed 74-marker
    inventory.""",
        ge=0,
        le=74,
    )
    checkeuk_expected_markers_present: int | None = Field(
        default=None,
        description="""Present markers in the CheckEUK clade-scored expected
    set.""",
        ge=0,
        le=74,
    )
    checkeuk_expected_markers_total: int | None = Field(
        default=None,
        description="""Denominator of the CheckEUK clade-scored expected marker
    set.""",
        ge=1,
        le=74,
    )
    checkeuk_pfam1085_distinct: int | None = Field(
        default=None,
        description="""Distinct recovered families in the global 1,085-family
    Pfam panel.""",
        ge=0,
        le=1085,
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class TaxonomyEvidence(ResultRecord):
    """Validated TaxonomyEvidence record."""

    record_id: str = Field(
        default=...,
        description="""Stable record identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    target_kind: TargetKind = Field(
        default=...,
        description="""Entity type evaluated by this record.""",
    )
    target_id: str = Field(
        default=...,
        description="""Identifier of the evaluated entity in the same sample
    and run.""",
        pattern=r"^\S+$",
    )
    rank: str = Field(
        default=...,
        description="""Taxonomic rank reported by the classification database.""",
        pattern=r"^\S+$",
    )
    taxon_name: str = Field(
        default=...,
        description="""Taxon label reported at rank.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    taxon_id: str | None = Field(
        default=None,
        description="""Database-specific taxon identifier, when available.""",
        pattern=r"^\S+$",
    )
    score: float | None = Field(
        default=None,
        description="""Raw tool score in the tool's declared units.""",
    )
    lineage_confidence: str | None = Field(
        default=None,
        description="""Verbatim GVClass label qualifying the full native
    lineage, not a calibrated per-rank probability.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class SsuLocus(ResultRecord):
    """Validated SsuLocus record."""

    record_id: str = Field(
        default=...,
        description="""Stable record identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    contig_id: str = Field(
        default=...,
        description="""Stable contig identifier scoped by sample_id and run_id.""",
        pattern=r"^\S+$",
    )
    ssu_type: SsuType = Field(
        default=...,
        description="""Ribosomal RNA locus type.""",
    )
    start: int = Field(
        default=...,
        description="""One-based inclusive feature start coordinate on the
    referenced sequence.""",
        ge=1,
    )
    end: int = Field(
        default=...,
        description="""One-based inclusive feature end coordinate on the
    referenced sequence.""",
        ge=1,
    )
    strand: Strand = Field(
        default=...,
        description="""Feature strand relative to the referenced sequence.""",
    )
    hit_accession: str | None = Field(
        default=None,
        description="""Best database hit accession, when a hit is available.""",
        pattern=r"^\S+$",
    )
    reference_identifiers: str | None = Field(
        default=None,
        description="""Native identifiers for candidate reference hits,
    retained as raw delimited text.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    reference_source: str | None = Field(
        default=None,
        description="""Native per-row source for candidate reference evidence,
    distinct from stage database provenance.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    reference_versions: str | None = Field(
        default=None,
        description="""Native per-row versions for candidate reference
    evidence, retained as raw delimited text.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    reference_taxonomy: str | None = Field(
        default=None,
        description="""Native candidate reference taxonomy, retained as raw
    decoded JSON array text.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    taxon_name: str | None = Field(
        default=None,
        description="""Taxon label reported at rank.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    taxonomy_assignment_method: str | None = Field(
        default=None,
        description="""Native per-row method explaining how the final
    taxon_name was assigned.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    identity_percent: float | None = Field(
        default=None,
        description="""Alignment identity in percent.""",
        ge=0,
        le=100,
    )
    query_coverage_percent: float | None = Field(
        default=None,
        description="""Fraction of the query aligned, expressed as percent.""",
        ge=0,
        le=100,
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class VirusCall(ResultRecord):
    """Validated VirusCall record."""

    record_id: str = Field(
        default=...,
        description="""Stable record identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    contig_id: str = Field(
        default=...,
        description="""Stable contig identifier scoped by sample_id and run_id.""",
        pattern=r"^\S+$",
    )
    start: int = Field(
        default=...,
        description="""One-based inclusive feature start coordinate on the
    referenced sequence.""",
        ge=1,
    )
    end: int = Field(
        default=...,
        description="""One-based inclusive feature end coordinate on the
    referenced sequence.""",
        ge=1,
    )
    strand: Strand = Field(
        default=...,
        description="""Feature strand relative to the referenced sequence.""",
    )
    score: float | None = Field(
        default=None,
        description="""Raw tool score in the tool's declared units.""",
    )
    viral_taxon: str | None = Field(
        default=None,
        description="""Viral lineage or label assigned by the caller.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Phenotype(ResultRecord):
    """Validated Phenotype record."""

    record_id: str = Field(
        default=...,
        description="""Stable record identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    bin_id: str = Field(
        default=...,
        description="""Stable bin identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    raw_class: str = Field(
        default=...,
        description="""Unmodified phenotype class emitted by the predictor.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    raw_score: float = Field(
        default=...,
        description="""Unmodified phenotype confidence score on a zero-to-one
    scale.""",
        ge=0,
        le=1,
    )
    threshold: float = Field(
        default=...,
        description="""Decision threshold applied to raw_score on a zero-to-one
    scale.""",
        ge=0,
        le=1,
    )
    passes_threshold: bool = Field(
        default=...,
        description="""Whether raw_score meets the recorded decision threshold.""",
    )
    thresholded_class: str | None = Field(
        default=None,
        description="""Predictor class after applying the recorded threshold,
    if any.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    is_applicable: bool = Field(
        default=...,
        description="""Whether workflow evidence permits biological
    interpretation for this target.""",
    )
    applicability_reason: str | None = Field(
        default=None,
        description="""Required explanation when a phenotype prediction is not
    applicable.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class Gene(ResultRecord):
    """Validated Gene record."""

    gene_id: str = Field(
        default=...,
        description="""Stable gene identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    contig_id: str = Field(
        default=...,
        description="""Stable contig identifier scoped by sample_id and run_id.""",
        pattern=r"^\S+$",
    )
    bin_id: str | None = Field(
        default=None,
        description="""Stable bin identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    start: int = Field(
        default=...,
        description="""One-based inclusive feature start coordinate on the
    referenced sequence.""",
        ge=1,
    )
    end: int = Field(
        default=...,
        description="""One-based inclusive feature end coordinate on the
    referenced sequence.""",
        ge=1,
    )
    strand: Strand = Field(
        default=...,
        description="""Feature strand relative to the referenced sequence.""",
    )
    phase: int | None = Field(
        default=None,
        description="""Coding-sequence phase as zero, one, or two bases.""",
        ge=0,
        le=2,
    )
    genetic_code: int = Field(
        default=...,
        description="""NCBI translation table number used or explicitly
    requested for gene calling.""",
        ge=1,
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class FunctionalAnnotation(ResultRecord):
    """Validated FunctionalAnnotation record."""

    annotation_id: str = Field(
        default=...,
        description="""Stable functional-annotation identifier scoped by
    sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    gene_id: str = Field(
        default=...,
        description="""Stable gene identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$",
    )
    database_name: str = Field(
        default=...,
        description="""Reference database name used by the stage, when
    applicable.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    database_version: str = Field(
        default=...,
        description="""Exact reference database release or immutable digest.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    accession: str = Field(
        default=...,
        description="""Functional database accession.""",
        pattern=r"^\S+$",
    )
    function_name: str | None = Field(
        default=None,
        description="""Function or domain name supplied by the source database.""",
        pattern=r"^\S(?:.*\S)?$",
    )
    start: int | None = Field(
        default=None,
        description="""One-based inclusive feature start coordinate on the
    referenced sequence.""",
        ge=1,
    )
    end: int | None = Field(
        default=None,
        description="""One-based inclusive feature end coordinate on the
    referenced sequence.""",
        ge=1,
    )
    score: float | None = Field(
        default=None,
        description="""Raw tool score in the tool's declared units.""",
    )
    evalue: float | None = Field(
        default=None,
        description="""Reported expectation value; lower values indicate
    stronger matches.""",
        ge=0,
    )
    stage_id: str = Field(
        default=...,
        description="""Stable stage identifier scoped by sample_id and run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    artifact_id: str = Field(
        default=...,
        description="""Stable artifact identifier scoped by sample_id and
    run_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    sample_id: str = Field(
        default=...,
        description="""Stable sample identifier unique within the catalog.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )
    run_id: str = Field(
        default=...,
        description="""Stable workflow run identifier scoped by sample_id.""",
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    )


class MetadataBundle(ConfiguredBaseModel):
    """Validated MetadataBundle record."""

    schema_version: str = Field(
        default=...,
        description="""Version of this complete bundle contract.""",
        pattern=r"^1\.0\.0$",
    )
    workflow_name: str = Field(
        default=...,
        description="""Canonical workflow name.""",
        pattern=r"^protist-meta-nf$",
    )
    workflow_version: str = Field(
        default=...,
        description="""Workflow version label that produced this bundle.""",
        pattern=r"^\S+$",
    )
    source_revision: str = Field(
        default=...,
        description="""Immutable source revision or working-tree content hash
    for this bundle.""",
        pattern=r"^\S+$",
    )
    samples: list[Sample] = Field(
        default=...,
        description="""Samples represented in the catalog.""",
    )
    runs: list[Run] = Field(
        default=...,
        description="""Sample-scoped workflow runs.""",
    )
    artifacts: list[Artifact] = Field(
        default=...,
        description="""Checksummed input and output files.""",
    )
    stages: list[Stage] = Field(
        default=...,
        description="""Stage states, provenance, and expected key coverage.""",
    )
    read_stats: list[ReadStatistics] = Field(
        default=...,
        description="""Read filtering and retention statistics.""",
    )
    assemblies: list[Assembly] = Field(
        default=...,
        description="""Assembly summaries.""",
    )
    contigs: list[Contig] = Field(
        default=...,
        description="""Assembly contigs.""",
    )
    bins: list[Bin] = Field(
        default=...,
        description="""QuickBin bins.""",
    )
    memberships: list[Membership] = Field(
        default=...,
        description="""Contig-to-bin memberships.""",
    )
    qc: list[QualityAssessment] = Field(
        default=...,
        description="""Independent tool-specific QC assessments.""",
    )
    taxonomy: list[TaxonomyEvidence] = Field(
        default=...,
        description="""Independent taxonomy evidence.""",
    )
    ssu: list[SsuLocus] = Field(
        default=...,
        description="""Assembly-wide ribosomal RNA loci.""",
    )
    viruses: list[VirusCall] = Field(
        default=...,
        description="""Assembly-wide viral calls.""",
    )
    phenotypes: list[Phenotype] = Field(
        default=...,
        description="""Raw and interpreted phenotype predictions.""",
    )
    genes: list[Gene] = Field(
        default=...,
        description="""Predicted genes.""",
    )
    annotations: list[FunctionalAnnotation] = Field(
        default=...,
        description="""Functional database hits.""",
    )


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
ScopedRecord.model_rebuild()
ResultRecord.model_rebuild()
Sample.model_rebuild()
Run.model_rebuild()
Artifact.model_rebuild()
Stage.model_rebuild()
ReadStatistics.model_rebuild()
Assembly.model_rebuild()
Contig.model_rebuild()
Bin.model_rebuild()
Membership.model_rebuild()
QualityAssessment.model_rebuild()
TaxonomyEvidence.model_rebuild()
SsuLocus.model_rebuild()
VirusCall.model_rebuild()
Phenotype.model_rebuild()
Gene.model_rebuild()
FunctionalAnnotation.model_rebuild()
MetadataBundle.model_rebuild()
