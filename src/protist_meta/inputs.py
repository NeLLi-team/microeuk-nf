"""Validate sample manifests and resolve registered dependencies before compute."""

import csv
import hashlib
import json
import os
import tomllib
from pathlib import Path
from typing import Literal, Self

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

__all__ = ["Registry", "SampleInput", "prepare_inputs", "read_registry", "read_samples"]

SAMPLE_COLUMNS = (
    "sample_id",
    "platform",
    "reads",
    "assembly",
    "rna_reads",
    "rna_platform",
    "protein_reference",
    "protein_lineage",
    "genetic_code",
    "softmasked",
    "basecaller",
    "basecaller_model",
    "library_prep",
)
PATH_COLUMNS = ("reads", "assembly", "rna_reads", "protein_reference", "softmasked")
UNSAFE_PATH_CHARACTERS = frozenset("\r\n\0'\"`$")


class SampleInput(BaseModel):
    """One DNA metagenome and its separately identified RNA evidence."""

    model_config = ConfigDict(extra="forbid", strict=True)

    sample_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$")
    platform: Literal["ONT", "PACBIO_HIFI"]
    reads: Path
    assembly: Path | None = None
    rna_reads: Path | None = None
    rna_platform: (
        Literal["ONT_CDNA", "ONT_DIRECT", "PACBIO_ISOSEQ", "ILLUMINA"] | None
    ) = None
    protein_reference: Path | None = None
    protein_lineage: str | None = Field(
        default=None, pattern=r"^[A-Za-z0-9][A-Za-z0-9_.:; /-]*$"
    )
    genetic_code: str = Field(
        default="auto",
        pattern=r"^(auto|[1-6]|9|1[0-6]|2[1-5])$",
        description="Prodigal-GV single-mode code, or automatic metagenomic selection.",
    )
    softmasked: Path | None = None
    basecaller: str | None = None
    basecaller_model: str | None = None
    library_prep: str | None = None

    @model_validator(mode="after")
    def check_evidence(self) -> Self:
        """Require evidence labels and prevent DNA/RNA path reuse."""
        if (self.rna_reads is None) != (self.rna_platform is None):
            raise ValueError("RNA reads and RNA platform must be supplied together")
        if self.rna_reads == self.reads:
            raise ValueError("RNA evidence cannot be the DNA read file")
        if (self.protein_reference is None) != (self.protein_lineage is None):
            raise ValueError("protein reference and lineage must be supplied together")
        if self.softmasked is not None and self.assembly is None:
            raise ValueError(
                "a supplied softmasked assembly needs the matching assembly"
            )
        if self.platform != "ONT" and (
            self.basecaller is not None or self.basecaller_model is not None
        ):
            raise ValueError("non-ONT sample contains ONT basecaller provenance")
        return self


class Dependency(BaseModel):
    """One explicit local tool or reference; no remote fallback is permitted."""

    model_config = ConfigDict(extra="forbid", strict=True)

    path: str
    kind: Literal["file", "directory"]
    version: str
    mode: Literal["core", "full"]
    sentinels: list[str] = Field(default_factory=list)
    executables: dict[str, list[str]] = Field(default_factory=dict)


class Registry(BaseModel):
    """Paths resolved relative to the registry, with recorded reference versions."""

    model_config = ConfigDict(extra="forbid", strict=True)

    schema_version: Literal["1.0.0"]
    entries: dict[str, Dependency]


def read_samples(path: Path) -> list[SampleInput]:
    """Read TSV and reject duplicate identifiers or unsupported input columns."""
    path = path.resolve()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        if not {"sample_id", "platform", "reads"}.issubset(columns):
            raise ValueError("sample sheet needs sample_id, platform and reads columns")
        if columns - set(SAMPLE_COLUMNS):
            raise ValueError("sample sheet contains unsupported columns")
        samples = [_parse_sample(row, path.parent) for row in reader]
    if not samples:
        raise ValueError("sample sheet contains no samples")
    identifiers = [sample.sample_id for sample in samples]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("sample identifiers must be unique")
    return samples


def read_registry(
    path: Path,
    *,
    mode: str,
    skip_gene_calling: bool = False,
    skip_annotation: bool = False,
) -> tuple[Registry, dict[str, Path]]:
    """Validate enabled local references without copying or downloading them."""
    path = path.resolve()
    registry = Registry.model_validate(yaml.safe_load(path.read_text(encoding="utf-8")))
    resolved: dict[str, Path] = {}
    excluded: set[str] = set()
    if skip_gene_calling:
        excluded.update(("gene_manifest", "braker_manifest", "dfam_db"))
    if skip_gene_calling or skip_annotation:
        excluded.update(
            ("annotation_manifest", "interpro_manifest", "eggnog_db", "interproscan_db")
        )
    for name, dependency in registry.entries.items():
        if name in excluded or (mode == "core" and dependency.mode == "full"):
            continue
        if not name.replace("_", "").isalnum():
            raise ValueError(f"invalid dependency parameter: {name}")
        target = _resolve_path(dependency.path, path.parent)
        if dependency.kind == "file":
            _require_file(target)
        elif not target.is_dir():
            raise FileNotFoundError(f"required directory is missing: {target}")
        for sentinel in dependency.sentinels:
            _require_file(target / sentinel)
        _check_executables(target, dependency)
        resolved[name] = target
    return registry, resolved


def _check_executables(target: Path, dependency: Dependency) -> None:
    root = target.parent if dependency.kind == "file" else target
    for environment, executables in dependency.executables.items():
        for executable in executables:
            command = root / ".pixi/envs" / environment / "bin" / executable
            _require_file(command)
            if not os.access(command, os.X_OK):
                raise ValueError(f"Pixi dependency is not executable: {command}")


def prepare_inputs(
    samples_path: Path,
    registry_path: Path,
    output_dir: Path,
    *,
    mode: str,
    skip_gene_calling: bool = False,
    skip_annotation: bool = False,
) -> Path:
    """Write validated TSV and a literal Nextflow dependency configuration."""
    samples = read_samples(samples_path)
    registry, dependencies = read_registry(
        registry_path,
        mode=mode,
        skip_gene_calling=skip_gene_calling,
        skip_annotation=skip_annotation,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = output_dir / "samples.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample.model_dump(mode="json") for sample in samples)
    config_lines: list[str] = [
        f"params.skip_gene_calling = {str(skip_gene_calling).lower()}",
        f"params.skip_annotation = {str(skip_annotation).lower()}",
        f"params.prepared_skip_gene_calling = {str(skip_gene_calling).lower()}",
        f"params.prepared_skip_annotation = {str(skip_annotation).lower()}",
    ]
    for name, path in dependencies.items():
        version = registry.entries[name].version
        if any(character in version for character in UNSAFE_PATH_CHARACTERS):
            raise ValueError(f"invalid characters in dependency version: {name}")
        config_lines.extend(
            (f"params.{name} = '{path}'", f"params.{name}_version = '{version}'")
        )
    config = "\n".join(config_lines)
    project = Path(__file__).resolve().parents[2]
    with (project / "pyproject.toml").open("rb") as handle:
        version = tomllib.load(handle)["project"]["version"]
    config += f"\nparams.workflow_version = '{version}'"
    config += f"\nparams.source_revision = '{_source_digest(project)}'"
    run_digest = hashlib.sha256(str(output_dir.resolve()).encode()).hexdigest()[:16]
    config += f"\nparams.run_id = 'run_{run_digest}'"
    (output_dir / "dependencies.config").write_text(config + "\n", encoding="utf-8")
    provenance = {
        "mode": mode,
        "skip_gene_calling": skip_gene_calling,
        "skip_annotation": skip_annotation,
        "registry": registry.model_dump(mode="json"),
        "resolved": {name: str(path) for name, path in dependencies.items()},
        "samples": [sample.model_dump(mode="json") for sample in samples],
    }
    (output_dir / "inputs.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def _parse_sample(row: dict[str, str], base: Path) -> SampleInput:
    values: dict[str, object] = {key: value for key, value in row.items() if value}
    for key in PATH_COLUMNS:
        if key in values:
            value = values[key]
            if not isinstance(value, str):
                raise ValueError(f"invalid sample path field: {key}")
            target = _resolve_path(value, base)
            _require_file(target)
            values[key] = target
    return SampleInput.model_validate(values)


def _resolve_path(value: str, base: Path) -> Path:
    if any(character in value for character in UNSAFE_PATH_CHARACTERS):
        raise ValueError("paths cannot contain quotes, shell expansions or line breaks")
    target = Path(value).expanduser()
    return (target if target.is_absolute() else base / target).resolve()


def _require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"required nonempty file is missing: {path}")


def _source_digest(root: Path) -> str:
    digest = hashlib.sha256()
    patterns = (
        "src/**/*.py",
        "modules/*.nf",
        "schema/*.yaml",
        "schema/templates/*.jinja",
        "*.nf",
        "*.config",
        "conf/*.config",
        "conf/*.yaml",
        "pixi.toml",
        "pixi.lock",
        "pyproject.toml",
        "bin/*",
        "scripts/*.py",
        "scripts/*.sh",
        "scripts/*.slurm",
        "envs/*/pixi.*",
        "envs/*/*/pixi.*",
        "envs/*/*/package/pixi.toml",
        "envs/*/*/package/recipe/*",
    )
    for pattern in patterns:
        for path in sorted(root.glob(pattern)):
            digest.update(str(path.relative_to(root)).encode())
            digest.update(path.read_bytes())
    registry_path = root / "conf/databases.yaml"
    entries = yaml.safe_load(registry_path.read_text(encoding="utf-8"))["entries"]
    for name, entry in sorted(entries.items()):
        if not name.endswith(("_manifest", "_app")):
            continue
        target = _resolve_path(entry["path"], registry_path.parent)
        manifest = target / "pixi.toml" if name.endswith("_app") else target
        for path in (manifest, manifest.with_name("pixi.lock")):
            digest.update(f"registered/{name}/{path.name}".encode())
            digest.update(path.read_bytes() if path.is_file() else b"MISSING")
    return "sha256:" + digest.hexdigest()
