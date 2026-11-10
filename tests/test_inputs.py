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
    Registry,
    _source_digest,
    prepare_inputs,
    read_registry,
    read_samples,
)

REGISTRY = Path(__file__).parents[1] / "conf/databases.yaml"


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
        ((), True, False, True),
        (("--skip-gene-calling",), True, False, True),
        (("--enable-gene-calling",), False, False, True),
        (("--skip-annotation",), True, True, True),
        (("--enable-gene-calling", "--skip-annotation"), False, True, True),
        (("--skip-gene-calling", "--skip-annotation"), True, True, True),
        (("--enable-repeat-masking",), True, False, False),
        (("--enable-gene-calling", "--enable-repeat-masking"), False, False, False),
        (("--skip-repeat-masking",), True, False, True),
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
    """Preparation records requested flags; gene calling is off unless enabled."""
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


def test_prepare_rejects_conflicting_gene_calling_options(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as error:
        main(
            [
                "prepare",
                "--input",
                str(tmp_path / "samples.tsv"),
                "--registry",
                str(tmp_path / "db.yaml"),
                "--output-dir",
                str(tmp_path / "prepared"),
                "--enable-gene-calling",
                "--skip-gene-calling",
            ]
        )

    assert error.value.code == 2


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
            prepare_inputs(
                sheet,
                registry,
                tmp_path / "disabled",
                mode=mode,
                skip_gene_calling=False,
            )
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
        skip_gene_calling=False,
        skip_repeat_masking=False,
    )
    assert manifest.is_file()


def test_prepare_accepts_softmasked_input_while_gene_calling_is_off_by_default(
    tmp_path: Path,
) -> None:
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

    manifest = prepare_inputs(sheet, registry, tmp_path / "prepared", mode="full")

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
            "",
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


@pytest.mark.parametrize(
    ("saved_skip", "requested", "option", "message"),
    [
        (False, "true", "", "required nonempty file is missing"),
        (True, "true", "", None),
        (True, "true", "--skip-gene-calling", None),
        (False, "false", "--enable-gene-calling", "required nonempty file is missing"),
        (False, "true", "--skip-gene-calling", "Resume stage choices differ"),
        (True, "false", "--enable-gene-calling", "Resume stage choices differ"),
    ],
    ids=[
        "adopts-enabled",
        "adopts-skipped",
        "repeats-skipped",
        "repeats-enabled",
        "explicit-skip-differs",
        "explicit-enable-differs",
    ],
)
def test_launcher_adopts_the_saved_gene_calling_choice_unless_explicit(
    tmp_path: Path,
    saved_skip: bool,
    requested: str,
    option: str,
    message: str | None,
) -> None:
    """Resume keeps a run's gene-calling choice after the default flipped to off.

    The registry lists an absent gene manifest, so only a validator that treats
    gene calling as enabled can notice it.
    """
    script = (Path(__file__).parents[1] / "scripts/run-allocation.slurm").read_text(
        encoding="utf-8"
    )
    validator = script.split("<<'PYTHON'\n")[1].split("\nPYTHON", 1)[0]
    registry = tmp_path / "db.yaml"
    registry.write_text(
        "schema_version: 1.0.0\nentries:\n  gene_manifest:\n"
        "    path: missing.toml\n    kind: file\n"
        "    version: fixture\n    mode: full\n",
        encoding="utf-8",
    )
    inputs = tmp_path / "inputs.json"
    saved = {"skip_gene_calling": saved_skip, "skip_repeat_masking": True}
    inputs.write_text(json.dumps({**saved, "samples": []}), encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            "-",
            str(registry),
            "full",
            requested,
            "false",
            str(inputs),
            "resume",
            "true",
            "",
            option,
        ],
        input=validator,
        text=True,
        capture_output=True,
        check=False,
    )

    if message is None:
        assert result.returncode == 0, result.stderr
    else:
        assert result.returncode != 0
        assert message in result.stderr


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


def _gene_registry(tmp_path: Path) -> Path:
    """A gene workspace that lacks the RepeatMasker environment and Dfam."""
    commands = tmp_path / "gene/.pixi/envs/prodigal-gv/bin"
    commands.mkdir(parents=True)
    command = commands / "prodigal-gv"
    command.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    command.chmod(0o755)
    (tmp_path / "gene/pixi.toml").write_text("# fixture manifest\n", encoding="utf-8")
    registry = tmp_path / "db.yaml"
    registry.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0.0",
                "entries": {
                    "gene_manifest": {
                        "path": "gene/pixi.toml",
                        "kind": "file",
                        "version": "fixture",
                        "mode": "full",
                        "executables": {
                            "repeatmasker": ["RepeatMasker"],
                            "prodigal-gv": ["prodigal-gv"],
                        },
                    },
                    "dfam_db": {
                        "path": "dfam",
                        "kind": "directory",
                        "version": "fixture",
                        "mode": "full",
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    return registry


def test_registry_skips_repeat_masking_dependencies_when_masking_is_off(
    tmp_path: Path,
) -> None:
    """Gene calling without repeat masking needs neither Dfam nor RepeatMasker."""
    registry = _gene_registry(tmp_path)

    _, resolved = read_registry(
        registry, mode="full", skip_gene_calling=False, skip_repeat_masking=True
    )

    assert set(resolved) == {"gene_manifest"}


def test_registry_requires_repeat_masking_dependencies_when_masking_is_on(
    tmp_path: Path,
) -> None:
    registry = _gene_registry(tmp_path)

    with pytest.raises(FileNotFoundError, match=r"dfam|RepeatMasker"):
        read_registry(
            registry, mode="full", skip_gene_calling=False, skip_repeat_masking=False
        )


def test_registry_ignores_repeat_masking_dependencies_by_default(
    tmp_path: Path,
) -> None:
    registry = _gene_registry(tmp_path)

    _, resolved = read_registry(registry, mode="full")

    assert resolved == {}


def _git(directory: Path, *arguments: str) -> str:
    return subprocess.run(
        [
            "git",
            "-C",
            str(directory),
            "-c",
            "user.name=fixture",
            "-c",
            "user.email=fixture@example.invalid",
            "-c",
            "commit.gpgsign=false",
            *arguments,
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _git_app(tmp_path: Path) -> tuple[Path, str]:
    """Create a one-commit application clone and return it with its commit."""
    app = tmp_path / "app"
    app.mkdir()
    (app / "pixi.toml").write_text("# fixture manifest\n", encoding="utf-8")
    _git(app, "init", "-q")
    _git(app, "add", "pixi.toml")
    _git(app, "commit", "-q", "-m", "fixture")
    return app, _git(app, "rev-parse", "HEAD")


def _app_registry(tmp_path: Path, app: Path, commit: str) -> Path:
    registry = tmp_path / "db.yaml"
    registry.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0.0",
                "entries": {
                    "fixture_app": {
                        "path": str(app),
                        "kind": "directory",
                        "version": "1.0",
                        "mode": "full",
                        "git_commit": commit,
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    return registry


def test_registry_accepts_an_app_at_its_pinned_clean_commit(tmp_path: Path) -> None:
    app, commit = _git_app(tmp_path)

    _, resolved = read_registry(_app_registry(tmp_path, app, commit), mode="full")

    assert resolved == {"fixture_app": app}


def test_registry_rejects_an_app_at_another_commit(tmp_path: Path) -> None:
    app, commit = _git_app(tmp_path)
    (app / "extra.txt").write_text("more\n", encoding="utf-8")
    _git(app, "add", "extra.txt")
    _git(app, "commit", "-q", "-m", "moved")

    with pytest.raises(ValueError, match=f"fixture_app is at commit .*pins {commit}"):
        read_registry(_app_registry(tmp_path, app, commit), mode="full")


@pytest.mark.parametrize("change", ["pixi.toml", "untracked.txt"])
def test_registry_rejects_an_app_with_uncommitted_changes(
    tmp_path: Path, change: str
) -> None:
    app, commit = _git_app(tmp_path)
    (app / change).write_text("changed\n", encoding="utf-8")

    with pytest.raises(ValueError, match="fixture_app has uncommitted changes"):
        read_registry(_app_registry(tmp_path, app, commit), mode="full")


@pytest.mark.parametrize("commit", ["86ab78a", "G" * 40, "$(id)"])
def test_registry_rejects_a_malformed_git_commit(tmp_path: Path, commit: str) -> None:
    app, _ = _git_app(tmp_path)

    with pytest.raises(ValueError, match="git_commit"):
        read_registry(_app_registry(tmp_path, app, commit), mode="full")


def test_registry_rejects_a_pinned_directory_that_is_not_a_clone(
    tmp_path: Path,
) -> None:
    app = tmp_path / "app"
    app.mkdir()

    with pytest.raises(ValueError, match="cannot read git state"):
        read_registry(_app_registry(tmp_path, app, "a" * 40), mode="full")


def test_prepare_writes_the_pinned_app_commit(tmp_path: Path) -> None:
    app, commit = _git_app(tmp_path)
    (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tplatform\treads\nS1\tONT\treads.fq\n", encoding="utf-8"
    )

    prepare_inputs(
        sheet,
        _app_registry(tmp_path, app, commit),
        tmp_path / "prepared",
        mode="full",
    )

    config = (tmp_path / "prepared/dependencies.config").read_text(encoding="utf-8")
    assert f"params.fixture_app_commit = '{commit}'\n" in config
    assert "params.fixture_app_version = '1.0'\n" in config


def test_source_digest_changes_when_an_app_moves_to_another_commit(
    tmp_path: Path,
) -> None:
    app, commit = _git_app(tmp_path)
    (tmp_path / "conf").mkdir()
    registry = _app_registry(tmp_path, app, commit)
    (tmp_path / "conf/databases.yaml").write_text(
        registry.read_text(encoding="utf-8"), encoding="utf-8"
    )
    original = _source_digest(tmp_path)
    (app / "extra.txt").write_text("more\n", encoding="utf-8")
    _git(app, "add", "extra.txt")
    _git(app, "commit", "-q", "-m", "moved")

    assert _source_digest(tmp_path) != original


def test_shipped_registry_pins_each_application_clone_to_a_commit() -> None:
    registry = Registry.model_validate(yaml.safe_load(REGISTRY.read_text("utf-8")))

    commits = {
        name: entry.git_commit
        for name, entry in registry.entries.items()
        if name.endswith("_app")
    }

    assert commits == {
        "checkeuk_app": "86ab78adab6b93031ef092a2dc47bcdd86243825",
        "gvclass_app": "1517cc057a484cace264e5853a94b0038b4d32b3",
        "ssuextract_app": "4dcf85d3796ef4ade3881594ad828428557d43f2",
    }


def test_shipped_registry_registers_the_tiara_and_busco_dependencies() -> None:
    registry = Registry.model_validate(yaml.safe_load(REGISTRY.read_text("utf-8")))

    assert {"contig_manifest", "eukqc_manifest", "busco_db"} <= set(registry.entries)
    assert registry.entries["busco_db"].sentinels == ["dataset.cfg"]
    assert registry.entries["eukqc_manifest"].executables == {
        "busco": ["busco", "python"],
        "tidk": ["tidk"],
    }
    assert registry.entries["viral_manifest"].executables == {
        "genomad": ["genomad", "python"]
    }
    assert registry.entries["checkv_manifest"].executables == {
        "checkv": ["checkv", "python"]
    }
