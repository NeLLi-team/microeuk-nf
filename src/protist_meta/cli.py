"""Command-line boundaries for manifest preparation and result publication."""

import argparse
import json
import logging
from collections.abc import Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

from protist_meta.collect import collect_records
from protist_meta.gene_commands import merge_proteins
from protist_meta.gene_inputs import prepare_gene_inputs
from protist_meta.inputs import prepare_inputs
from protist_meta.visualization import export_visualization

__all__ = ["main"]


def main(argv: Sequence[str] | None = None) -> int:
    """Execute one explicit operation, returning nonzero on invalid input."""
    parser = _parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    match args.command:
        case "prepare":
            prepare_inputs(
                args.input,
                args.registry,
                args.output_dir,
                mode=args.mode,
                skip_gene_calling=args.skip_gene_calling,
                skip_annotation=args.skip_annotation,
            )
        case "collect":
            collect_records(args.sample_json, args.stages, args.output)
        case "prepare-gene-inputs":
            prepare_gene_inputs(args.routing_dir, args.viral_dir, args.output_dir)
        case "merge-proteins":
            merge_proteins(args.prodigal_dir, args.braker_dir, args.output_dir)
        case "route":
            from protist_meta.routing import (  # noqa: PLC0415  # Routing is independent of catalog model generation.
                route_bins,
            )

            route_bins(
                args.sample_json,
                args.bins_dir,
                {
                    "quickclade": args.quickclade_dir,
                    "checkeuk": args.checkeuk_dir,
                    "gvclass": args.gvclass_dir,
                    "ssuextract": args.ssuextract_dir,
                    "viral": args.viral_dir,
                },
                args.output_dir,
            )
        case "catalog":
            _build_catalog(args.records, args.output)
        case "visualization":
            export_visualization(args.catalog, args.output_dir)
        case "report":
            from protist_meta.report import (  # noqa: PLC0415  # Avoid notebook startup cost for input validation.
                build_report,
            )

            build_report(
                args.catalog,
                args.output_dir,
                execution_provenance=args.execution_provenance,
            )
    return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser(
        "prepare", help="Validate inputs and local Pixi dependencies"
    )
    prepare.add_argument("--input", type=Path, required=True, help="Sample TSV")
    prepare.add_argument(
        "--registry", type=Path, required=True, help="Local database/tool YAML"
    )
    prepare.add_argument(
        "--output-dir", type=Path, required=True, help="Prepared run directory"
    )
    prepare.add_argument(
        "--mode", choices=("core", "full"), default="full", help="Analysis scope"
    )
    prepare.add_argument(
        "--skip-gene-calling",
        action="store_true",
        help="Skip dedicated gene calling and functional annotation",
    )
    prepare.add_argument(
        "--skip-annotation",
        action="store_true",
        help="Skip functional annotation while retaining gene calling",
    )
    collect = commands.add_parser("collect", help="Normalize exact stage outputs")
    collect.add_argument(
        "--sample-json", type=Path, required=True, help="Sample/run metadata"
    )
    collect.add_argument(
        "--stages", type=Path, nargs="+", required=True, help="Stage directories"
    )
    collect.add_argument(
        "--output", type=Path, required=True, help="Normalized bundle JSON"
    )
    gene_inputs = commands.add_parser(
        "prepare-gene-inputs", help="Stage exclusive bin and geNomad gene inputs"
    )
    gene_inputs.add_argument("--routing-dir", type=Path, required=True)
    gene_inputs.add_argument("--viral-dir", type=Path, required=True)
    gene_inputs.add_argument("--output-dir", type=Path, required=True)
    proteins = commands.add_parser(
        "merge-proteins", help="Create uniquely named proteins and their source map"
    )
    proteins.add_argument("--prodigal-dir", type=Path, required=True)
    proteins.add_argument("--braker-dir", type=Path, required=True)
    proteins.add_argument("--output-dir", type=Path, required=True)
    route = commands.add_parser("route", help="Choose gene callers from bin evidence")
    route.add_argument(
        "--sample-json", type=Path, required=True, help="Sample/run metadata"
    )
    route.add_argument(
        "--bins-dir", type=Path, required=True, help="QuickBin FASTA directory"
    )
    for tool in ("quickclade", "checkeuk", "gvclass", "ssuextract", "viral"):
        route.add_argument(
            f"--{tool}-dir", type=Path, required=True, help=f"Native {tool} results"
        )
    route.add_argument(
        "--output-dir", type=Path, required=True, help="Bin copies and routing evidence"
    )
    catalog = commands.add_parser(
        "catalog", help="Validate and atomically publish SQLite"
    )
    catalog.add_argument(
        "--records",
        type=Path,
        nargs="+",
        required=True,
        help="Normalized sample bundle JSON files",
    )
    catalog.add_argument("--output", type=Path, required=True, help="SQLite output")
    report = commands.add_parser(
        "report", help="Execute database-backed notebook and HTML"
    )
    report.add_argument(
        "--catalog", type=Path, required=True, help="Validated SQLite catalog"
    )
    report.add_argument(
        "--output-dir", type=Path, required=True, help="Report directory"
    )
    report.add_argument(
        "--execution-provenance",
        type=Path,
        help="Reviewed-source recovery receipt to validate and embed",
    )
    visualization = commands.add_parser(
        "visualization", help="Portable native evidence"
    )
    exports = visualization.add_subparsers(dest="operation", required=True)
    export = exports.add_parser("export", help="Export native visualization bundles")
    export.add_argument(
        "--catalog", type=Path, required=True, help="Validated SQLite catalog"
    )
    export.add_argument(
        "--output-dir", type=Path, required=True, help="Bundle directory"
    )
    return parser


def _build_catalog(paths: list[Path], output: Path) -> None:
    from protist_meta.catalog import (  # noqa: PLC0415  # Input checks do not need generated models.
        build_catalog,
        validate_bundle,
    )

    if len(paths) == 1:
        build_catalog(paths[0], output)
        return
    bundles = [validate_bundle(path).model_dump(mode="json") for path in paths]
    metadata_fields = (
        "schema_version",
        "workflow_name",
        "workflow_version",
        "source_revision",
    )
    metadata = {field: bundles[0][field] for field in metadata_fields}
    if any(
        bundle[field] != metadata[field]
        for bundle in bundles
        for field in metadata_fields
    ):
        raise ValueError(
            "sample bundles have different schema or workflow source identities"
        )
    merged = dict(metadata)
    for field, value in bundles[0].items():
        if isinstance(value, list):
            merged[field] = [record for bundle in bundles for record in bundle[field]]
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="catalog-merge-", dir=output.parent) as temporary:
        path = Path(temporary) / "records.json"
        path.write_text(json.dumps(merged), encoding="utf-8")
        build_catalog(path, output)


if __name__ == "__main__":
    raise SystemExit(main())
