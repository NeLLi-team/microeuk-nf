#!/usr/bin/env python3
"""Generate strict Pydantic models from the catalog LinkML schema."""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import tempfile
from contextlib import chdir
from pathlib import Path
from types import ModuleType

from linkml.generators.pydanticgen import PydanticGenerator

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCHEMA = PROJECT_ROOT / "schema" / "catalog.yaml"
DEFAULT_OUTPUT = PROJECT_ROOT / "src" / "protist_meta" / "models.py"
EXPECTED_CLASSES = (
    "MetadataBundle",
    "Sample",
    "Run",
    "Artifact",
    "Stage",
    "ReadStatistics",
    "Assembly",
    "Contig",
    "Bin",
    "Membership",
    "QualityAssessment",
    "TaxonomyEvidence",
    "SsuLocus",
    "VirusCall",
    "Phenotype",
    "Gene",
    "FunctionalAnnotation",
)


def generate_source(schema_path: Path) -> str:
    """Render importable Pydantic source with strict boundary validation."""
    resolved_schema = schema_path.resolve()
    with chdir(resolved_schema.parent):
        source = PydanticGenerator(
            resolved_schema.name,
            extra_fields="forbid",
            metadata_mode=None,
            template_dir=PROJECT_ROOT / "schema" / "templates",
        ).serialize()
    compile(source, str(resolved_schema), "exec")
    return source


def format_source(source: str, output_path: Path) -> str:
    """Format generated source with the project's pinned Ruff configuration."""
    for operation in (("format",), ("check", "--fix"), ("format",)):
        completed = subprocess.run(
            [
                "ruff",
                *operation,
                "--config",
                str(PROJECT_ROOT / "ruff.toml"),
                "--stdin-filename",
                str(output_path),
                "-",
            ],
            input=source,
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"generated model violates house style: {completed.stderr}"
            )
        source = completed.stdout
    compile(source, str(output_path), "exec")
    return source


def load_module(path: Path) -> ModuleType:
    """Import a generated module and return it after executing class definitions."""
    module_name = f"generated_{path.stem}"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load generated model: {path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)
    return module


def publish_source(output_path: Path, source: str) -> str:
    """Atomically publish generated source, or leave an identical file unchanged."""
    if output_path.exists() and output_path.read_text(encoding="utf-8") == source:
        return "unchanged"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        handle.write(source)
        temporary_path = Path(handle.name)
    temporary_path.replace(output_path)
    return "generated"


def check_source(output_path: Path, source: str) -> str:
    """Reject a missing or stale generated module without changing it."""
    if not output_path.is_file():
        raise FileNotFoundError(f"generated model is missing: {output_path}")
    if output_path.read_text(encoding="utf-8") != source:
        raise ValueError(f"generated model differs from schema: {output_path}")
    return "verified"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse model-generation paths and verification mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--schema",
        type=Path,
        default=DEFAULT_SCHEMA,
        help="LinkML schema to generate from",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Generated Python module to write or verify",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify that the generated module matches the schema without writing",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Generate or verify the catalog models and import the resulting module."""
    args = parse_args(argv)
    schema_path = args.schema.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    if not schema_path.is_file():
        raise FileNotFoundError(f"schema does not exist: {schema_path}")

    source = format_source(generate_source(schema_path), output_path)
    status = (
        check_source(output_path, source)
        if args.check
        else publish_source(output_path, source)
    )
    module = load_module(output_path)
    missing_classes = [
        class_name for class_name in EXPECTED_CLASSES if not hasattr(module, class_name)
    ]
    if missing_classes:
        names = ", ".join(missing_classes)
        raise RuntimeError(f"generated model is missing expected classes: {names}")

    print(f"Catalog model {status}: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
