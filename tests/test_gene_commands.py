"""Tests for preparing functional-annotation protein inputs."""

from __future__ import annotations

from pathlib import Path

import pytest

from protist_meta.gene_commands import merge_proteins


def _write_source_ledger(path: Path) -> None:
    path.parent.mkdir(parents=True)
    path.write_text(
        "source_id\tcategory\tbin_id\tfasta\torigin_map\tstatus\treason\n"
        "short-virus\tviral_candidate\t\tshort-virus.fna\t"
        "short-virus.origins.tsv\tincluded\t\n",
        encoding="utf-8",
    )


def test_merge_proteins_preserves_zero_call_source(tmp_path: Path) -> None:
    prodigal = tmp_path / "prodigal"
    braker = tmp_path / "braker"
    output = tmp_path / "output"
    _write_source_ledger(prodigal / "inputs" / "sources.tsv")
    (prodigal / "proteins").mkdir()
    (prodigal / "proteins" / "short-virus.faa").touch()
    (braker / "proteins").mkdir(parents=True)

    result = merge_proteins(prodigal, braker, output)

    assert result == output / "proteins.faa"
    assert result.read_text(encoding="utf-8") == ""
    assert (output / "protein-map.tsv").read_text(encoding="utf-8") == (
        "normalized_id\tnative_id\tnative_gene_id\tcaller\tsource_id\tbin_id\n"
    )


def test_merge_proteins_rejects_missing_zero_call_fasta(tmp_path: Path) -> None:
    prodigal = tmp_path / "prodigal"
    braker = tmp_path / "braker"
    _write_source_ledger(prodigal / "inputs" / "sources.tsv")
    (prodigal / "proteins").mkdir()
    (braker / "proteins").mkdir(parents=True)

    with pytest.raises(ValueError, match=r"missing=\['short-virus'\]"):
        merge_proteins(prodigal, braker, tmp_path / "output")
