"""Input-boundary checks for sequencing technology and sample identity."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from protist_meta.cli import main
from protist_meta.inputs import (
    _source_digest,
    prepare_inputs,
    read_registry,
    read_samples,
)


def test_source_digest_changes_when_python_runtime_helper_changes(tmp_path):
    (tmp_path / "conf").mkdir()
    (tmp_path / "conf/databases.yaml").write_text("entries: {}\n", encoding="utf-8")
    (tmp_path / "scripts").mkdir()
    runtime_helper = tmp_path / "scripts/check-gvclass-reference.py"
    runtime_helper.write_text("raise SystemExit(0)\n", encoding="utf-8")
    original = _source_digest(tmp_path)

    runtime_helper.write_text("raise SystemExit(1)\n", encoding="utf-8")
    changed = _source_digest(tmp_path)

    assert changed != original


def test_manifest_resolves_relative_reads(tmp_path):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\nS1\tONT\treads.fq\n", encoding="utf-8"
    )

    samples = read_samples(sheet)

    assert samples[0].reads == tmp_path / "reads.fq"
    assert samples[0].genetic_code == "auto"


@pytest.mark.parametrize("platform", ["PACBIO_CLR", "ILLUMINA"])
def test_manifest_rejects_unsupported_dna_platform(tmp_path, platform):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"sample_id\tplatform\treads\nS1\t{platform}\treads.fq\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match="platform"):
        read_samples(sheet)


@pytest.mark.parametrize(
    ("platform", "basecaller", "model", "message"),
    [
        ("PACBIO_HIFI", "Dorado", "", "non-ONT sample"),
        ("PACBIO_HIFI", "Dorado", "sup_model", "non-ONT sample"),
    ],
)
def test_manifest_rejects_invalid_basecaller_provenance(
    tmp_path, platform, basecaller, model, message
):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\tbasecaller\tbasecaller_model\n"
        f"S1\t{platform}\treads.fq\t{basecaller}\t{model}\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=message):
        read_samples(sheet)


@pytest.mark.parametrize("basecaller", ["Dorado", ""])
def test_manifest_preserves_ont_basecaller_provenance(tmp_path, basecaller):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\tbasecaller\tbasecaller_model\n"
        f"S1\tONT\treads.fq\t{basecaller}\tsup_model\n",
        encoding="utf-8",
    )

    sample = read_samples(sheet)[0]

    assert sample.basecaller == (basecaller or None)
    assert sample.basecaller_model == "sup_model"


def test_manifest_rejects_dna_as_rna_evidence(tmp_path):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\trna_reads\trna_platform\n"
        "S1\tONT\treads.fq\treads.fq\tONT_CDNA\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="RNA evidence"):
        read_samples(sheet)


@pytest.mark.parametrize("code", ["29", "30", "106", "129"])
def test_manifest_rejects_unsupported_explicit_gene_codes(tmp_path, code):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"sample_id\tplatform\treads\tgenetic_code\nS1\tONT\treads.fq\t{code}\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="genetic_code"):
        read_samples(sheet)


@pytest.mark.parametrize(
    ("dependency", "skip_gene_calling", "skip_annotation"),
    [
        ("quickclade_ref", False, False),
        ("quickclade_ref", True, True),
        ("gene_manifest", False, False),
        ("gene_manifest", False, True),
        ("eggnog_db", False, False),
    ],
)
def test_registry_fails_without_local_reference(
    tmp_path: Path,
    dependency: str,
    *,
    skip_gene_calling: bool,
    skip_annotation: bool,
) -> None:
    """Enabled dependencies and QC references still require local files."""
    registry = tmp_path / "db.yaml"
    registry.write_text(
        f"schema_version: 1.0.0\nentries:\n  {dependency}:\n"
        "    path: missing.spectra.gz\n    kind: file\n"
        "    version: fixture\n    mode: full\n",
        encoding="utf-8",
    )

    with pytest.raises(FileNotFoundError, match=r"missing\.spectra"):
        read_registry(
            registry,
            mode="full",
            skip_gene_calling=skip_gene_calling,
            skip_annotation=skip_annotation,
        )


@pytest.mark.parametrize(
    ("disabled", "skip_gene_calling", "skip_annotation"),
    [
        (
            (
                "gene_manifest",
                "braker_manifest",
                "dfam_db",
                "annotation_manifest",
                "interpro_manifest",
                "eggnog_db",
                "interproscan_db",
            ),
            True,
            False,
        ),
        (
            (
                "annotation_manifest",
                "interpro_manifest",
                "eggnog_db",
                "interproscan_db",
            ),
            False,
            True,
        ),
    ],
)
def test_registry_skips_only_disabled_stage_dependencies(
    tmp_path: Path,
    disabled: tuple[str, ...],
    *,
    skip_gene_calling: bool,
    skip_annotation: bool,
) -> None:
    """Disabled references may be absent while shared dependencies remain resolved."""
    reference = tmp_path / "present.spectra.gz"
    reference.write_text("reference\n", encoding="utf-8")
    registry = tmp_path / "db.yaml"
    entries = {
        name: {"path": "missing", "kind": "file", "version": "fixture", "mode": "full"}
        for name in disabled
    }
    entries["quickclade_ref"] = {
        "path": reference.name,
        "kind": "file",
        "version": "fixture",
        "mode": "core",
    }
    registry.write_text(
        yaml.safe_dump({"schema_version": "1.0.0", "entries": entries}),
        encoding="utf-8",
    )

    _, resolved = read_registry(
        registry,
        mode="full",
        skip_gene_calling=skip_gene_calling,
        skip_annotation=skip_annotation,
    )

    assert resolved == {"quickclade_ref": reference}


@pytest.mark.parametrize(
    ("flags", "skip_gene_calling", "skip_annotation", "skip_repeat_masking"),
    [
        ((), False, False, True),
        (("--skip-gene-calling",), True, False, True),
        (("--skip-annotation",), False, True, True),
        (("--skip-gene-calling", "--skip-annotation"), True, True, True),
        (("--enable-repeat-masking",), False, False, False),
        (("--skip-repeat-masking",), False, False, True),
    ],
)
def test_prepare_records_requested_stage_choices(
    tmp_path: Path,
    flags: tuple[str, ...],
    *,
    skip_gene_calling: bool,
    skip_annotation: bool,
    skip_repeat_masking: bool,
) -> None:
    """Preparation records requested flags, including the enabled legacy defaults."""
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\nS1\tONT\treads.fq\n", encoding="utf-8"
    )
    registry = tmp_path / "db.yaml"
    registry.write_text("schema_version: 1.0.0\nentries: {}\n", encoding="utf-8")
    output = tmp_path / "prepared"

    status = main(
        [
            "prepare",
            "--input",
            str(sheet),
            "--registry",
            str(registry),
            "--output-dir",
            str(output),
            *flags,
        ]
    )

    assert status == 0
    provenance = json.loads((output / "inputs.json").read_text(encoding="utf-8"))
    assert provenance["skip_gene_calling"] is skip_gene_calling
    assert provenance["skip_annotation"] is skip_annotation
    assert provenance["skip_repeat_masking"] is skip_repeat_masking
    config = (output / "dependencies.config").read_text(encoding="utf-8")
    assert f"params.skip_gene_calling = {str(skip_gene_calling).lower()}\n" in config
    assert f"params.skip_annotation = {str(skip_annotation).lower()}\n" in config
    assert (
        f"params.skip_repeat_masking = {str(skip_repeat_masking).lower()}\n" in config
    )
    assert (
        f"params.prepared_skip_repeat_masking = {str(skip_repeat_masking).lower()}\n"
        in config
    )
    assert (
        f"params.prepared_skip_gene_calling = {str(skip_gene_calling).lower()}\n"
        in config
    )
    assert (
        f"params.prepared_skip_annotation = {str(skip_annotation).lower()}\n" in config
    )


@pytest.mark.parametrize(
    ("mode", "skip_genes"),
    [("full", False), ("core", False), ("full", True)],
)
def test_softmasked_input_requires_repeat_masking(
    tmp_path: Path, mode: str, skip_genes: bool
):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    (tmp_path / "masked.fa").write_text(">c\nACgt\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\tassembly\tsoftmasked\n"
        "S1\tONT\treads.fq\tmasked.fa\tmasked.fa\n",
        encoding="utf-8",
    )
    registry = tmp_path / "db.yaml"
    registry.write_text("schema_version: 1.0.0\nentries: {}\n", encoding="utf-8")

    if mode == "full" and not skip_genes:
        with pytest.raises(ValueError, match=r"softmasked.*enable-repeat-masking"):
            prepare_inputs(sheet, registry, tmp_path / "disabled", mode=mode)
    else:
        assert prepare_inputs(
            sheet,
            registry,
            tmp_path / "disabled",
            mode=mode,
            skip_gene_calling=skip_genes,
        ).is_file()
    manifest = prepare_inputs(
        sheet,
        registry,
        tmp_path / "enabled",
        mode="full",
        skip_repeat_masking=False,
    )
    assert manifest.is_file()


@pytest.mark.parametrize(
    ("saved", "mode", "option", "requested", "expected"),
    [
        ({}, "resume", "", "true", "false"),
        ({}, "resume-reviewed", "", "true", "false"),
        ({}, "resume", "--skip-repeat-masking", "true", None),
        ({}, "resume-reviewed", "--skip-repeat-masking", "true", "true"),
        ({"skip_repeat_masking": True}, "resume", "", "true", "true"),
        (
            {"skip_repeat_masking": False},
            "resume",
            "--enable-repeat-masking",
            "false",
            "false",
        ),
        (
            {"skip_repeat_masking": True},
            "resume-reviewed",
            "--enable-repeat-masking",
            "false",
            None,
        ),
        (
            {"mode": "core", "samples": [{"softmasked": "supplied.fa"}]},
            "resume",
            "--skip-repeat-masking",
            "true",
            "true",
        ),
        (
            {"skip_gene_calling": True, "samples": [{"softmasked": "supplied.fa"}]},
            "resume",
            "--skip-repeat-masking",
            "true",
            "true",
        ),
    ],
)
def test_launcher_preserves_saved_repeat_choice_unless_reviewed_disable(
    tmp_path: Path,
    saved: dict,
    mode: str,
    option: str,
    requested: str,
    expected: str | None,
) -> None:
    """Execute the launcher's actual preparation validator with legacy metadata."""
    script = (Path(__file__).parents[1] / "scripts/run-allocation.slurm").read_text(
        encoding="utf-8"
    )
    validator = script.split("<<'PYTHON'\n")[1].split("\nPYTHON", 1)[0]
    registry = tmp_path / "db.yaml"
    registry.write_text("schema_version: 1.0.0\nentries: {}\n", encoding="utf-8")
    inputs = tmp_path / "inputs.json"
    inputs.write_text(json.dumps({"samples": [], **saved}), encoding="utf-8")
    before = inputs.read_bytes()

    result = subprocess.run(
        [
            sys.executable,
            "-",
            str(registry),
            saved.get("mode", "full"),
            str(saved.get("skip_gene_calling", False)).lower(),
            "false",
            str(inputs),
            mode,
            requested,
            option,
        ],
        input=validator,
        text=True,
        capture_output=True,
        check=False,
    )

    assert inputs.read_bytes() == before
    if expected is None:
        assert result.returncode != 0
        assert "Repeat choice differs" in result.stderr
    else:
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == expected


def test_reviewed_receipt_records_repeat_override_without_rewriting_preparation(
    tmp_path: Path,
) -> None:
    script = (Path(__file__).parents[1] / "scripts/run-allocation.slurm").read_text(
        encoding="utf-8"
    )
    receipt_code = script.split("<<'PYTHON'\n")[2].split("\nPYTHON", 1)[0]
    identity = tmp_path / "resume-identity.tsv"
    identity.write_text("source_revision\told-source\n", encoding="utf-8")
    inputs = tmp_path / "inputs.json"
    inputs.write_text('{"samples":[{"sample_id":"S1"}]}\n', encoding="utf-8")
    dependencies = tmp_path / "dependencies.config"
    dependencies.write_text("params.run_id = 'saved-run'\n", encoding="utf-8")
    before = (identity.read_bytes(), inputs.read_bytes(), dependencies.read_bytes())
    execution = tmp_path / "execution"
    execution.mkdir()
    config = execution / "allocation.config"
    config.write_text(
        "params.skip_repeat_masking = true\n"
        "params.prepared_skip_repeat_masking = true\n",
        encoding="utf-8",
    )

    subprocess.run(
        [
            sys.executable,
            "-",
            str(identity),
            str(execution),
            "old-source",
            "new-source",
            "true",
        ],
        input=receipt_code,
        text=True,
        capture_output=True,
        check=True,
    )

    receipt = json.loads((execution / "recovery.json").read_text(encoding="utf-8"))
    assert receipt["prepared_skip_repeat_masking"] is False
    assert receipt["effective_skip_repeat_masking"] is True
    assert receipt["original_source_revision"] == "old-source"
    assert receipt["execution_source_revision"] == "new-source"
    assert (
        receipt["allocation_config_sha256"]
        == hashlib.sha256(config.read_bytes()).hexdigest()
    )
    assert (
        identity.read_bytes(),
        inputs.read_bytes(),
        dependencies.read_bytes(),
    ) == before


def test_manifest_rejects_shell_expansion_in_protein_lineage(tmp_path):
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    (tmp_path / "reference.faa").write_text(">protein\nMPEPTIDE\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\tprotein_reference\tprotein_lineage\n"
        "S1\tONT\treads.fq\treference.faa\t$(printf injected)\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="protein_lineage"):
        read_samples(sheet)


@pytest.mark.parametrize("kind", ["file", "directory"])
def test_registry_checks_tools_inside_each_pixi_workspace(tmp_path, kind):
    workspace = tmp_path / "tool"
    commands = workspace / ".pixi/envs/default/bin"
    commands.mkdir(parents=True)
    command = commands / "tool-command"
    command.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    command.chmod(0o755)
    manifest = workspace / "pixi.toml"
    manifest.write_text("# fixture manifest\n", encoding="utf-8")
    target = "tool/pixi.toml" if kind == "file" else "tool"
    registry = tmp_path / "db.yaml"
    registry.write_text(
        "schema_version: 1.0.0\nentries:\n  test_tool:\n"
        f"    path: {target}\n    kind: {kind}\n"
        "    version: fixture\n    mode: core\n"
        "    executables:\n      default: [tool-command]\n",
        encoding="utf-8",
    )

    read_registry(registry, mode="core")
    command.chmod(0o644)

    with pytest.raises(ValueError, match="not executable"):
        read_registry(registry, mode="core")
