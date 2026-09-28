"""Full collector integration tests for gene calls and functional annotations."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from protist_meta.catalog import build_catalog, validate_bundle
from protist_meta.collect import collect_records
from protist_meta.report import build_report

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


def test_collects_lifted_genes_partial_braker_and_annotations(
    tmp_path: Path,
) -> None:
    sample, stages = _full_fixture(tmp_path)
    bundle_path = collect_records(sample, stages, tmp_path / "bundle.json")

    bundle = validate_bundle(bundle_path)
    database = build_catalog(bundle_path, tmp_path / "catalog.sqlite")

    stages_by_id = {stage.stage_id: stage for stage in bundle.stages}
    assert stages_by_id["braker3"].status == "pending"
    assert stages_by_id["braker3.euk1"].status == "completed"
    assert "braker3.euk2" not in stages_by_id
    assert stages_by_id["braker3.euk1"].expected_result_keys == [
        "genes:braker3:braker3:euk1:g1.t1"
    ]

    genes = {gene.gene_id: gene for gene in bundle.genes}
    lifted = genes["prodigal-gv:genomad-unbinned:viral_fragment_1"]
    assert lifted.contig_id == "virus_parent"
    assert (lifted.start, lifted.end, lifted.genetic_code) == (104, 400, 15)
    assert set(stages_by_id["prodigal_gv"].expected_result_keys) == {
        "genes:prodigal-gv:genomad-unbinned:viral_fragment_1",
        "genes:prodigal-gv:prok-source:prok:prok_contig_1",
    }

    assert {annotation.stage_id for annotation in bundle.annotations} == {
        "eggnog_mapper",
        "interproscan",
    }
    assert {annotation.gene_id for annotation in bundle.annotations} == {
        "braker3:braker3:euk1:g1.t1",
        "prodigal-gv:prok-source:prok:prok_contig_1",
    }
    parent_outputs = {
        "functional_annotation.proteins",
        "functional_annotation.protein_map",
    }
    assert set(stages_by_id["eggnog_mapper"].input_artifact_ids) == parent_outputs
    assert set(stages_by_id["interproscan"].input_artifact_ids) == parent_outputs

    artifact_kinds = {str(artifact.kind) for artifact in bundle.artifacts}
    assert {"alignment", "gene_annotation", "functional_annotation"}.issubset(
        artifact_kinds
    )
    with sqlite3.connect(database) as connection:
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute("SELECT COUNT(*) FROM genes").fetchone() == (3,)
        assert connection.execute("SELECT COUNT(*) FROM annotations").fetchone() == (2,)


@pytest.mark.parametrize(
    ("skip_genes", "skip_annotation", "stage_count", "gene_count", "reason"),
    [
        (
            False,
            True,
            8,
            3,
            "Functional annotation disabled by user (--skip-annotation).",
        ),
        (True, False, 4, 0, "Gene calling disabled by user (--skip-gene-calling)."),
        (True, True, 4, 0, "Gene calling disabled by user (--skip-gene-calling)."),
    ],
)
def test_disabled_stages_survive_collection_catalog_and_report(
    tmp_path: Path,
    skip_genes: bool,
    skip_annotation: bool,
    stage_count: int,
    gene_count: int,
    reason: str,
) -> None:
    """Explicit skips preserve upstream results without fabricated products."""
    sample, stages = _full_fixture(tmp_path)
    metadata = json.loads(sample.read_text(encoding="utf-8"))
    metadata.update(skip_gene_calling=skip_genes, skip_annotation=skip_annotation)
    sample.write_text(json.dumps(metadata), encoding="utf-8")

    bundle_path = collect_records(
        sample, stages[:stage_count], tmp_path / "bundle.json"
    )
    bundle = validate_bundle(bundle_path)
    database = build_catalog(bundle_path, tmp_path / "catalog.sqlite")
    notebook = build_report(database, tmp_path / "report")

    skipped = {
        stage.stage_id: stage for stage in bundle.stages if stage.status == "skipped"
    }
    assert {"functional_annotation", "eggnog_mapper", "interproscan"} <= skipped.keys()
    assert {stage.status_reason for stage in skipped.values()} == {reason}
    assert {stage.tool_version for stage in skipped.values()} == {"not_run"}
    assert all(stage.expected_result_keys == [] for stage in skipped.values())
    assert all(stage.input_artifact_ids == [] for stage in skipped.values())
    assert (
        not {artifact.producer_stage_id for artifact in bundle.artifacts}
        & skipped.keys()
    )
    assert len(bundle.genes) == gene_count
    assert bundle.annotations == []
    assert len(bundle.bins) == 3
    assert len(bundle.contigs) == 4
    assert ("braker3.euk1" in skipped) == skip_genes
    assert ("braker3.euk2" in skipped) == skip_genes
    assert notebook.is_file()
    assert reason in (tmp_path / "report/index.html").read_text(encoding="utf-8")
    with sqlite3.connect(database) as connection:
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute("SELECT COUNT(*) FROM genes").fetchone() == (
            gene_count,
        )
        assert connection.execute("SELECT COUNT(*) FROM annotations").fetchone() == (0,)


@pytest.mark.parametrize(
    ("flag", "stage_name"),
    [
        ("skip_gene_calling", "repeat_masking"),
        ("skip_gene_calling", "braker3.euk1"),
        ("skip_gene_calling", "functional_annotation"),
        ("skip_annotation", "functional_annotation"),
        ("skip_annotation", "eggnog_mapper"),
    ],
)
def test_disabled_stages_reject_supplied_manifests(
    tmp_path: Path, flag: str, stage_name: str
) -> None:
    """A disabled manifest cannot silently populate supposedly skipped results."""
    sample, stages = _full_fixture(tmp_path)
    metadata = json.loads(sample.read_text(encoding="utf-8"))
    metadata[flag] = True
    sample.write_text(json.dumps(metadata), encoding="utf-8")
    stage = _stage(
        tmp_path / "stages",
        "disabled",
        {
            "stage": stage_name,
            "status": "completed",
            "tool": "fixture",
            "tool_version": "fixture",
            "command": "fixture",
        },
    )

    with pytest.raises(ValueError, match=f"supplied stage {stage_name} is disabled"):
        collect_records(sample, [*stages[:4], stage], tmp_path / "bundle.json")


@pytest.mark.parametrize("flag", ["skip_gene_calling", "skip_annotation"])
def test_skip_flags_require_json_booleans(tmp_path: Path, flag: str) -> None:
    """String booleans cannot invert requested execution at collection."""
    sample, stages = _full_fixture(tmp_path)
    metadata = json.loads(sample.read_text(encoding="utf-8"))
    metadata[flag] = "false"
    sample.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(TypeError, match=f"{flag} must be a boolean"):
        collect_records(sample, stages, tmp_path / "bundle.json")


def test_pending_braker_requires_exact_bin_coverage(tmp_path: Path) -> None:
    sample, stages = _full_fixture(tmp_path)
    (tmp_path / "stages/06_braker3/pending.tsv").unlink()

    with pytest.raises(ValueError, match="do not cover routed eukaryotic bins"):
        collect_records(sample, stages, tmp_path / "bundle.json")


def test_collects_checkm2_contamination_above_100(tmp_path: Path) -> None:
    """Native estimates above 100 survive collection and the catalog unchanged."""
    sample, stages = _full_fixture(tmp_path)
    screen = _stage(
        tmp_path / "stages",
        "03b_checkm2",
        {
            "stage": "checkm2",
            "status": "completed",
            "tool": "CheckM2",
            "tool_version": "1.1.0",
            "command": "checkm2 predict",
            "outputs": {"result": "quality_report.tsv"},
        },
    )
    _write(
        screen / "quality_report.tsv",
        "Name\tCompleteness\tContamination\tCompleteness_Model_Used\n"
        "euk1\t90\t103.62\tNeural Network (Specific Model)\n"
        "euk2\t80\t101.71\tNeural Network (Specific Model)\n"
        "prok\t70\t101.77\tNeural Network (Specific Model)\n",
    )
    bundle_path = collect_records(
        sample, [*stages[:3], screen], tmp_path / "bundle.json"
    )
    database = build_catalog(bundle_path, tmp_path / "catalog.sqlite")

    with sqlite3.connect(database) as connection:
        assert connection.execute(
            "SELECT target_id, contamination_percent FROM qc ORDER BY target_id"
        ).fetchall() == [("euk1", 103.62), ("euk2", 101.71), ("prok", 101.77)]
        with pytest.raises(sqlite3.IntegrityError, match="CHECK constraint"):
            connection.execute(
                "UPDATE qc SET contamination_percent = ?", (float("inf"),)
            )


def test_collects_nextflow_command_for_parent_and_child_stages(
    tmp_path: Path,
) -> None:
    """Resolved producer scripts override summary text for grouped stages."""
    sample, stages = _full_fixture(tmp_path)
    commands = {
        "06_braker3": "#!/usr/bin/env bash\nset -euo pipefail\nbraker.pl --gff3",
        "08_functional_annotation": (
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            "emapper.py --output annotations && interproscan.sh"
        ),
    }
    for stage_name, command in commands.items():
        original = tmp_path / "stages" / stage_name
        producer = tmp_path / "work" / stage_name
        producer.mkdir(parents=True)
        produced_stage = producer / stage_name
        original.rename(produced_stage)
        (producer / ".command.sh").write_text(f"\n{command}\n", encoding="utf-8")
        original.symlink_to(produced_stage, target_is_directory=True)

    bundle_path = collect_records(sample, stages, tmp_path / "bundle.json")
    bundle = validate_bundle(bundle_path)
    stages_by_id = {stage.stage_id: stage for stage in bundle.stages}

    assert stages_by_id["braker3"].command == commands["06_braker3"]
    assert stages_by_id["braker3.euk1"].command == commands["06_braker3"]
    assert stages_by_id["eggnog_mapper"].command == commands["08_functional_annotation"]
    assert stages_by_id["interproscan"].command == commands["08_functional_annotation"]


def test_completed_zero_bin_stage_omits_quickclade(tmp_path: Path) -> None:
    sample, stages = _full_fixture(tmp_path)
    binning = tmp_path / "stages/03_binning"
    for path in (binning / "bins").glob("*.fa"):
        path.unlink()
    (binning / "quickclade.tsv").unlink()
    manifest_path = binning / "stage.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"].pop("quickclade")
    manifest["outputs"].update(
        {
            "quickbin_report": "quickbin_report.tsv",
            "quickbin_coverage": "quickbin_coverage.tsv",
        }
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _write(binning / "quickbin_report.tsv", "# no accepted bins\n")
    _write(binning / "quickbin_coverage.tsv", "# no accepted bins\n")

    screen = _stage(
        tmp_path / "stages",
        "03b_all_bin_screen",
        {
            "stage": "all_bin_screen",
            "status": "skipped",
            "tool": "multi-tool",
            "tool_version": "not_run",
            "command": "skip screens when QuickBin accepts no bins",
            "reason": "no_bins",
            "outputs": {},
            "tools": {
                name: {
                    "version": "not_run",
                    "output": name,
                    "status": "skipped",
                    "reason": "no_bins",
                }
                for name in ("checkeuk", "checkm2", "gvclass")
            },
        },
    )
    for name in ("checkeuk", "checkm2", "gvclass"):
        (screen / name).mkdir()
    routing = _stage(
        tmp_path / "stages",
        "03c_routing",
        {
            "stage": "routing",
            "status": "completed",
            "tool": "protist-meta",
            "tool_version": "fixture",
            "command": "write empty routing evidence",
            "outputs": {"evidence": "evidence.tsv"},
        },
    )
    for name in ("eukaryotic", "prokaryotic", "viral", "unresolved"):
        (routing / name).mkdir()
    _write(routing / "evidence.tsv", "bin_id\tcandidate_class\n")

    bundle_path = collect_records(
        sample,
        [*stages[:3], screen, routing],
        tmp_path / "bundle.json",
    )
    bundle = validate_bundle(bundle_path)

    binning_stage = next(
        stage for stage in bundle.stages if stage.stage_id == "binning"
    )
    assert binning_stage.status == "completed"
    assert binning_stage.expected_result_keys == []
    assert bundle.bins == []
    assert bundle.memberships == []
    assert bundle.taxonomy == []
    routing_stage = next(
        stage for stage in bundle.stages if stage.stage_id == "routing"
    )
    assert routing_stage.status == "completed"
    for name in ("checkeuk", "checkm2", "gvclass"):
        stage = next(stage for stage in bundle.stages if stage.stage_id == name)
        assert stage.status == "skipped"
        assert stage.status_reason == "no_bins"


def test_skipped_checkm2_retains_no_annotation_evidence(tmp_path: Path) -> None:
    """An eligible CheckM2 no-result run remains queryable without QC rows."""
    sample, stages = _full_fixture(tmp_path)
    screen = _checkm2_no_annotation_fixture(tmp_path / "stages")

    bundle_path = collect_records(
        sample,
        [*stages[:3], screen],
        tmp_path / "bundle.json",
    )
    bundle = validate_bundle(bundle_path)
    database = build_catalog(bundle_path, tmp_path / "catalog.sqlite")

    stages_by_id = {stage.stage_id: stage for stage in bundle.stages}
    checkm2 = stages_by_id["checkm2"]
    assert checkm2.status == "skipped"
    assert checkm2.status_reason == "no_diamond_annotations"
    assert checkm2.expected_result_keys == []
    assert checkm2.notes == (
        "CheckM2 ran gene calling and DIAMOND, but DIAMOND returned no "
        "annotations; CheckM2 quality prediction is unavailable."
    )
    assert stages_by_id["checkeuk"].status == "completed"
    assert stages_by_id["gvclass"].status == "completed"
    assert not [row for row in bundle.qc if row.stage_id == "checkm2"]

    diagnostics = {
        Path(artifact.path).relative_to(screen.resolve()).as_posix(): artifact
        for artifact in bundle.artifacts
        if artifact.producer_stage_id == "checkm2"
        and artifact.artifact_id != "checkm2.manifest"
    }
    assert set(diagnostics) == {
        "checkm2/command.log",
        "checkm2/exit_code.txt",
        "checkm2/output/checkm2.log",
        "checkm2/output/diamond_output/DIAMOND_RESULTS.tsv",
        "checkm2/output/protein_files/euk1.faa",
        "checkm2/output/protein_files/euk2.faa",
        "checkm2/output/protein_files/prok.faa",
    }
    assert all(
        artifact.producer_stage_id == "checkm2" for artifact in diagnostics.values()
    )
    assert (
        diagnostics["checkm2/output/diamond_output/DIAMOND_RESULTS.tsv"].size_bytes == 0
    )
    assert not any("checkeuk" in path or "gvclass" in path for path in diagnostics)
    with sqlite3.connect(database) as connection:
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute(
            "SELECT COUNT(*) FROM qc WHERE stage_id = 'checkm2'"
        ).fetchone() == (0,)
        assert connection.execute(
            "SELECT COUNT(*) FROM artifacts WHERE producer_stage_id = 'checkm2'"
        ).fetchone() == (8,)


def test_zero_gene_source_builds_valid_catalog_with_original_contigs(
    tmp_path: Path,
) -> None:
    sample, stages = _full_fixture(tmp_path)
    prodigal = tmp_path / "stages/07_prodigal_gv"
    annotation = tmp_path / "stages/08_functional_annotation"

    sources = (prodigal / "inputs/sources.tsv").read_text(encoding="utf-8")
    _write(
        prodigal / "inputs/sources.tsv",
        "\n".join(
            line
            for line in sources.splitlines()
            if line.startswith(("source_id\t", "genomad-unbinned\t"))
        )
        + "\n",
    )
    origins = (prodigal / "inputs/origins.tsv").read_text(encoding="utf-8")
    _write(
        prodigal / "inputs/origins.tsv",
        "\n".join(
            line
            for line in origins.splitlines()
            if line.startswith(("contig_id\t", "viral_fragment\t"))
        )
        + "\n",
    )
    for relative in (
        "inputs/sources/prok-source.fna",
        "inputs/origins/prok-source.tsv",
        "proteins/prok-source.faa",
        "genes/prok-source.fna",
        "gff/prok-source.gff",
    ):
        (prodigal / relative).unlink()
    _write(
        prodigal / "gff/genomad-unbinned.gff",
        "##gff-version  3\n"
        '# Sequence Data: seqnum=1;seqlen=300;seqhdr="viral_fragment"\n'
        "# Model Data: version=Prodigal.v2.11.0-gv;"
        "run_type=Metagenomic;model=fixture;gc_cont=29.30;"
        "transl_table=4;uses_sd=1\n",
    )
    _write(prodigal / "proteins/genomad-unbinned.faa", "")
    _write(prodigal / "genes/genomad-unbinned.fna", "")

    manifest_path = annotation / "stage.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["status"] = "skipped"
    manifest["reason"] = "no_called_proteins"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _write(annotation / "proteins.faa", "")
    _write(
        annotation / "protein-map.tsv",
        "normalized_id\tnative_id\tnative_gene_id\tcaller\tsource_id\tbin_id\n",
    )
    (annotation / "eggnog/annotations.emapper.annotations").unlink()
    (annotation / "interpro/proteins.tsv").unlink()

    stages = [stage for stage in stages if stage.name != "06_braker3"]
    bundle_path = collect_records(sample, stages, tmp_path / "bundle.json")
    bundle = validate_bundle(bundle_path)
    database = build_catalog(bundle_path, tmp_path / "catalog.sqlite")

    stage = next(
        stage for stage in bundle.stages if stage.stage_id == "functional_annotation"
    )
    assert stage.status == "skipped"
    assert stage.status_reason == "no_called_proteins"
    assert stage.expected_result_keys == []
    assert bundle.genes == []
    assert bundle.annotations == []
    assert len(bundle.contigs) == 4
    artifacts = {artifact.artifact_id: artifact for artifact in bundle.artifacts}
    assert artifacts["functional_annotation.proteins"].size_bytes == 0
    assert artifacts["functional_annotation.protein_map"].size_bytes > 0
    with sqlite3.connect(database) as connection:
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute("SELECT COUNT(*) FROM contigs").fetchone() == (4,)
        assert connection.execute("SELECT COUNT(*) FROM genes").fetchone() == (0,)
        assert connection.execute("SELECT COUNT(*) FROM annotations").fetchone() == (0,)


def _full_fixture(tmp_path: Path) -> tuple[Path, list[Path]]:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    reads = _write(inputs / "reads.fastq.gz", "fixture reads\n")
    rna = _write(inputs / "rna.fastq.gz", "fixture RNA\n")
    proteins = _write(inputs / "reference.faa", ">reference\nMPEPTIDE\n")
    stages_root = tmp_path / "stages"
    stages_root.mkdir()

    read_qc = _stage(
        stages_root,
        "01_read_qc",
        {
            "stage": "read_qc",
            "status": "completed",
            "tool": "seqkit",
            "tool_version": "2.11.0",
            "command": "seqkit stats",
            "outputs": {},
        },
    )
    stats_header = "file\tformat\ttype\tnum_seqs\tsum_len\n"
    _write(read_qc / "raw_stats.tsv", stats_header + "raw\tFASTQ\tDNA\t10\t4000\n")
    _write(
        read_qc / "filtered_stats.tsv",
        stats_header + "filtered\tFASTQ\tDNA\t8\t3200\n",
    )
    _write(read_qc / "filtered.fastq.gz", "filtered reads\n")

    assembly = _stage(
        stages_root,
        "02_assembly",
        {
            "stage": "assembly",
            "status": "completed",
            "tool": "Flye",
            "tool_version": "2.9.6",
            "command": "flye --meta",
            "outputs": {"assembly": "assembly_primary.fa"},
        },
    )
    sequences = {
        "euk_contig_1": "A" * 1000,
        "euk_contig_2": "C" * 1000,
        "prok_contig": "G" * 1000,
        "virus_parent": "T" * 1000,
    }
    assembly_fasta = _write(
        assembly / "assembly_primary.fa",
        "".join(f">{name}\n{sequence}\n" for name, sequence in sequences.items()),
    )

    binning = _stage(
        stages_root,
        "03_binning",
        {
            "stage": "binning",
            "status": "completed",
            "tool": "QuickBin",
            "tool_version": "BBTools 40.02",
            "database_name": "quickclade_ref",
            "database_version": "fixture",
            "command": "quickbin.sh; quickclade.sh",
            "outputs": {
                "coverage": "coverage.tsv",
                "bins": "bins",
                "quickclade": "quickclade.tsv",
            },
        },
    )
    bins = binning / "bins"
    bins.mkdir()
    _write(bins / "euk1.fa", f">euk_contig_1\n{sequences['euk_contig_1']}\n")
    _write(bins / "euk2.fa", f">euk_contig_2\n{sequences['euk_contig_2']}\n")
    _write(bins / "prok.fa", f">prok_contig\n{sequences['prok_contig']}\n")
    coverage_header = (
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\n"
    )
    _write(
        binning / "coverage.tsv",
        coverage_header
        + "".join(f"{name}\t1\t1000\t2\t1000\t100\t2.0\n" for name in sequences),
    )
    _write(
        binning / "quickclade.tsv",
        "#QueryName\tlineage\tR_TaxID\n"
        "euk1.fa\tEukaryota\t2759\n"
        "euk2.fa\tEukaryota\t2759\n"
        "prok.fa\tBacteria\t2\n",
    )

    routing = _stage(
        stages_root,
        "04_routing",
        {
            "stage": "routing",
            "status": "completed",
            "tool": "protist-meta",
            "tool_version": "fixture",
            "command": "route bins",
            "outputs": {"evidence": "evidence.tsv"},
        },
    )
    _write(
        routing / "evidence.tsv",
        "bin_id\tcandidate_class\n"
        "euk1\teukaryotic_candidate\n"
        "euk2\teukaryotic_candidate\n"
        "prok\tbacterial_archaeal_candidate\n",
    )

    repeat = _stage(
        stages_root,
        "05_repeat_masking",
        {
            "stage": "repeat_masking",
            "status": "completed",
            "tool": "RepeatModeler;RepeatMasker",
            "tool_version": "2.0.9;4.2.4",
            "database_name": "de_novo_repeat_model",
            "database_version": "per-bin",
            "command": "RepeatModeler; RepeatMasker -xsmall",
            "outputs": {
                "masked_bins": "masked",
                "repeat_models": "models",
                "logs": "logs",
            },
        },
    )
    for name in ("masked", "models", "logs"):
        (repeat / name).mkdir()
    _write(
        repeat / "masked/euk1.softmasked.fna",
        f">euk_contig_1\n{sequences['euk_contig_1'].lower()}\n",
    )
    _write(
        repeat / "masked/euk2.softmasked.fna",
        f">euk_contig_2\n{sequences['euk_contig_2'].lower()}\n",
    )

    rna_stage = _stage(
        stages_root,
        "05b_rna_alignment",
        {
            "stage": "rna_alignment",
            "status": "completed",
            "tool": "minimap2;samtools",
            "tool_version": "2.30;1.22",
            "command": "minimap2 -ax splice; samtools sort",
            "outputs": {
                "whole_assembly_bam": "sample.splice.primary.unique.bam",
                "bin_bams": "bins",
                "stats": "stats",
            },
        },
    )
    (rna_stage / "bins").mkdir()
    (rna_stage / "stats").mkdir()
    _write(rna_stage / "sample.splice.primary.unique.bam", "bam\n")
    _write(rna_stage / "sample.splice.primary.unique.bam.bai", "index\n")
    _write(rna_stage / "bins/euk1.bam", "bam\n")
    _write(rna_stage / "bins/euk1.bam.bai", "index\n")
    _write(rna_stage / "bins/euk2.bam", "empty bam\n")
    _write(rna_stage / "stats/retained.flagstat.txt", "1 + 0 mapped\n")

    braker = _braker_fixture(stages_root)
    prodigal = _prodigal_fixture(stages_root)
    annotation = _annotation_fixture(stages_root)

    sample = tmp_path / "sample.json"
    sample.write_text(
        json.dumps(
            {
                "sample_id": "sample1",
                "run_id": "run1",
                "platform": "ONT",
                "reads": str(reads),
                "assembly": str(assembly_fasta),
                "rna_reads": str(rna),
                "rna_platform": "ont_cdna",
                "protein_reference": str(proteins),
                "genetic_code": 11,
                "workflow_version": "0.1.0",
                "source_revision": "0123456789abcdef",
            }
        ),
        encoding="utf-8",
    )
    return sample, [
        read_qc,
        assembly,
        binning,
        routing,
        repeat,
        rna_stage,
        braker,
        prodigal,
        annotation,
    ]


def _braker_fixture(root: Path) -> Path:
    stage = _stage(
        root,
        "06_braker3",
        {
            "stage": "braker3",
            "status": "pending",
            "tool": "BRAKER3",
            "tool_version": "3.0.8-longrna-source-lock",
            "database_name": "protein_reference",
            "database_version": "fixture lineage",
            "command": "braker.pl --gff3 --softmasking",
            "reason": "one_or_more_bins_lack_compatible_rna_or_protein_evidence",
            "outputs": {"annotations": "annotations", "proteins": "proteins"},
        },
    )
    (stage / "annotations/euk1").mkdir(parents=True)
    (stage / "proteins").mkdir()
    _write(
        stage / "annotations/euk1/braker.gff3",
        "##gff-version 3\n"
        "euk_contig_1\tAUGUSTUS\tgene\t100\t399\t.\t+\t.\tID=g1;\n"
        "euk_contig_1\tAUGUSTUS\tmRNA\t100\t399\t.\t+\t.\tID=g1.t1;Parent=g1;\n"
        "euk_contig_1\tAUGUSTUS\tCDS\t100\t399\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1;\n",
    )
    _write(stage / "proteins/euk1.faa", ">g1.t1\nMPEPTIDE\n")
    _write(stage / "pending.tsv", "euk2\tno_compatible_rna_or_protein_evidence\n")
    return stage


def _checkm2_no_annotation_fixture(root: Path) -> Path:
    stage = _stage(
        root,
        "03b_all_bin_screen",
        {
            "stage": "all_bin_screen",
            "status": "completed",
            "tool": "CheckEUK;CheckM2;GVClass",
            "tool_version": "0.12;1.1.0;2.1",
            "database_name": "checkeuk_db;checkm2_db;gvclass_db",
            "database_version": "fixture;fixture;fixture",
            "command": "run all-bin screens and retain native CheckM2 evidence",
            "tools": {
                "checkeuk": {
                    "version": "0.12",
                    "database_name": "checkeuk_db",
                    "database_version": "fixture",
                    "output": "checkeuk/result_checkeuk/checkeuk_report.tsv",
                },
                "checkm2": {
                    "version": "1.1.0",
                    "status": "skipped",
                    "reason": "no_diamond_annotations",
                    "database_name": "checkm2_db",
                    "database_version": "fixture",
                    "output": "checkm2",
                },
                "gvclass": {
                    "version": "2.1",
                    "database_name": "gvclass_db",
                    "database_version": "fixture",
                    "output": "gvclass/output/gvclass_summary.tsv",
                },
            },
            "outputs": {
                "checkeuk": "checkeuk/result_checkeuk/checkeuk_report.tsv",
                "checkm2": "checkm2",
                "gvclass": "gvclass/output/gvclass_summary.tsv",
            },
        },
    )
    checkeuk = stage / "checkeuk/result_checkeuk"
    checkeuk.mkdir(parents=True)
    checkm2 = stage / "checkm2/output"
    (checkm2 / "diamond_output").mkdir(parents=True)
    (checkm2 / "protein_files").mkdir()
    gvclass = stage / "gvclass/output"
    gvclass.mkdir(parents=True)

    checkeuk_header = (
        "genome\tstatus\tcompleteness\tcompleteness_basis\tcompleteness_bom\t"
        "contamination\tlineage\tcompleteness_note\tpfam_status\t"
        "pfam_union_families\tpfam_observed_families\tpfam_completeness\n"
    )
    checkeuk_row = (
        "\tok\t0\tpfam_union1085_ridge\t\t0\tk__Eukaryota\t"
        "markers 0/74\tscored\t1085\t0\t0\n"
    )
    _write(
        checkeuk / "checkeuk_report.tsv",
        checkeuk_header
        + "euk1.fa"
        + checkeuk_row
        + "euk2.fa"
        + checkeuk_row
        + "prok.fa"
        + checkeuk_row,
    )
    _write(
        checkeuk / "checkeuk_report_detail.tsv",
        "genome\tn_present_core\tn_present_expected\tn_expected\n"
        "euk1.fa\t0\t0\t74\n"
        "euk2.fa\t0\t0\t74\n"
        "prok.fa\t0\t0\t74\n",
    )

    no_annotations = "No DIAMOND annotation was generated. Exiting\n"
    _write(stage / "checkm2/command.log", no_annotations)
    _write(stage / "checkm2/exit_code.txt", "1\n")
    _write(checkm2 / "checkm2.log", no_annotations)
    _write(checkm2 / "diamond_output/DIAMOND_RESULTS.tsv", "")
    _write(checkm2 / "protein_files/euk1.faa", ">euk1_1\nMPEPTIDE\n")
    _write(checkm2 / "protein_files/euk2.faa", ">euk2_1\nMPEPTIDE\n")
    _write(checkm2 / "protein_files/prok.faa", ">prok_1\nMPEPTIDE\n")

    gvclass_header = (
        "query\ttaxonomy_majority\ttaxonomy_confidence\t"
        "estimated_completeness\testimated_completeness_strategy\t"
        "completeness_query_support\tcompleteness_model_group\t"
        "completeness_model_reliability\testimated_contamination\t"
        "estimated_contamination_strategy\tcontamination_query_support\t"
        "contamination_model_reliability\tcontamination_type\n"
    )
    gvclass_row = (
        "\tk__Eukaryota\thigh\t50\tgvog8 regression\t1\tNCLDV\thigh\t0\t"
        "marker duplication\t1\thigh\tnone\n"
    )
    _write(
        gvclass / "gvclass_summary.tsv",
        gvclass_header
        + "euk1.fa"
        + gvclass_row
        + "euk2.fa"
        + gvclass_row
        + "prok.fa"
        + gvclass_row,
    )
    return stage


def _prodigal_fixture(root: Path) -> Path:
    stage = _stage(
        root,
        "07_prodigal_gv",
        {
            "stage": "prodigal_gv",
            "status": "completed",
            "tool": "Prodigal-gv",
            "tool_version": "2.11.0",
            "database_name": "prodigal_gv_models",
            "database_version": "bundled viral models",
            "command": "prodigal-gv -p meta",
            "outputs": {
                "sources": "inputs/sources.tsv",
                "origins": "inputs/origins.tsv",
                "proteins": "proteins",
                "genes": "genes",
                "gff": "gff",
            },
        },
    )
    for name in ("inputs/sources", "inputs/origins", "proteins", "genes", "gff"):
        (stage / name).mkdir(parents=True, exist_ok=True)
    header = "source_id\tcategory\tbin_id\tfasta\torigin_map\tstatus\treason\n"
    _write(
        stage / "inputs/sources.tsv",
        header + "prok-source\tprokaryotic_bin\tprok\tsources/prok-source.fna\t"
        "origins/prok-source.tsv\tincluded\t\n"
        + "genomad-unbinned\tgenomad_unbinned\t\tsources/genomad-unbinned.fna\t"
        "origins/genomad-unbinned.tsv\tincluded\t\n",
    )
    origin_header = "contig_id\tparent_contig\tstart\tend\tsource_id\t"
    origin_header += "bin_id\tcategory\tstatus\treason\n"
    prok_origin = (
        "prok_contig\tprok_contig\t1\t1000\tprok-source\tprok\t"
        "prokaryotic_bin\tincluded\t\n"
    )
    virus_origin = (
        "viral_fragment\tvirus_parent\t101\t400\tgenomad-unbinned\t\t"
        "genomad_unbinned\tincluded\t\n"
    )
    _write(stage / "inputs/origins.tsv", origin_header + prok_origin + virus_origin)
    _write(stage / "inputs/origins/prok-source.tsv", origin_header + prok_origin)
    _write(
        stage / "inputs/origins/genomad-unbinned.tsv",
        origin_header + virus_origin,
    )
    _write(
        stage / "inputs/sources/prok-source.fna",
        ">prok_contig\n" + "G" * 1000 + "\n",
    )
    _write(
        stage / "inputs/sources/genomad-unbinned.fna",
        ">viral_fragment\n" + "T" * 300 + "\n",
    )
    _write_prodigal_call(stage, "prok-source", "prok_contig", 1, 300, 11)
    _write_prodigal_call(stage, "genomad-unbinned", "viral_fragment", 4, 300, 15)
    return stage


def _write_prodigal_call(
    stage: Path,
    source: str,
    contig: str,
    start: int,
    end: int,
    code: int,
) -> None:
    _write(
        stage / f"gff/{source}.gff",
        "##gff-version  3\n"
        f"# Model Data: version=Prodigal.v2.11.0-gv;transl_table={code};uses_sd=1\n"
        f"{contig}\tProdigal_v2.11.0-gv\tCDS\t{start}\t{end}\t.\t+\t0\t"
        f"ID=1_1;partial=00;start_type=ATG;genetic_code={code}\n",
    )
    _write(
        stage / f"proteins/{source}.faa",
        f">{contig}_1 # {start} # {end} # 1 # "
        f"ID=1_1;partial=00;start_type=ATG;genetic_code={code}\nMPEPTIDE\n",
    )
    _write(stage / f"genes/{source}.fna", f">{contig}_1\n" + "A" * 300 + "\n")


def _annotation_fixture(root: Path) -> Path:
    stage = _stage(
        root,
        "08_functional_annotation",
        {
            "stage": "functional_annotation",
            "status": "completed",
            "tool": "eggNOG-mapper;InterProScan",
            "tool_version": "2.1.15;5.76-107.0",
            "database_name": "eggnog_db;interproscan_db",
            "database_version": "eggNOG 5.0.2;InterPro 107.0",
            "command": "emapper.py; interproscan.sh",
            "tools": {
                "eggnog_mapper": {
                    "version": "2.1.15",
                    "database_name": "eggnog_db",
                    "database_version": "5.0.2",
                    "output": "eggnog/annotations.emapper.annotations",
                },
                "interproscan": {
                    "version": "5.76-107.0",
                    "database_name": "interproscan_db",
                    "database_version": "107.0",
                    "output": "interpro",
                },
            },
            "outputs": {
                "proteins": "proteins.faa",
                "protein_map": "protein-map.tsv",
                "eggnog": "eggnog/annotations.emapper.annotations",
                "interpro": "interpro",
            },
        },
    )
    (stage / "eggnog").mkdir()
    (stage / "interpro").mkdir()
    braker = "braker3:braker3:euk1:g1.t1"
    prok = "prodigal-gv:prok-source:prok:prok_contig_1"
    virus = "prodigal-gv:genomad-unbinned:viral_fragment_1"
    _write(
        stage / "proteins.faa",
        f">{braker}\nMPEPTIDE\n>{prok}\nMPEPTIDE\n>{virus}\nMPEPTIDE\n",
    )
    map_header = "normalized_id\tnative_id\tnative_gene_id\tcaller\tsource_id\tbin_id\n"
    _write(
        stage / "protein-map.tsv",
        map_header
        + f"{braker}\tg1.t1\tg1\tbraker3\tbraker3\teuk1\n"
        + f"{prok}\tprok_contig_1\t1_1\tprodigal-gv\tprok-source\tprok\n"
        + f"{virus}\tviral_fragment_1\t1_1\tprodigal-gv\tgenomad-unbinned\t\n",
    )
    values = dict.fromkeys(EGGNOG_HEADER, "-")
    values.update(
        {
            "#query": braker,
            "evalue": "1e-40",
            "score": "200",
            "eggNOG_OGs": "KOG0001@2759|Eukaryota",
            "Description": "fixture protein",
        }
    )
    _write(
        stage / "eggnog/annotations.emapper.annotations",
        "## emapper-2.1.15\n"
        + "\t".join(EGGNOG_HEADER)
        + "\n"
        + "\t".join(values[field] for field in EGGNOG_HEADER)
        + "\n",
    )
    _write(
        stage / "interpro/proteins.tsv",
        f"{prok}\t0123456789abcdef0123456789abcdef\t100\tPfam\tPF00001\t"
        "Fixture domain\t3\t45\t2.4E-7\tT\t19-09-2026\n",
    )
    return stage


def _stage(root: Path, name: str, manifest: dict[str, object]) -> Path:
    directory = root / name
    directory.mkdir()
    (directory / "stage.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    return directory


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path
