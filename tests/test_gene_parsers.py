"""Tests for native gene-call parsing and annotation protein naming."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from protist_meta.gene_parsers import (
    ProteinSource,
    merge_annotation_proteins,
    parse_braker3,
    parse_prodigal_gv,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_parse_prodigal_gv_preserves_native_codes_and_coordinates(
    tmp_path: Path,
) -> None:
    gff = _write(
        tmp_path / "genes.gff",
        (
            "##gff-version  3\n"
            '# Sequence Data: seqnum=1;seqlen=1200;seqhdr="viral_contig"\n'
            "# Model Data: version=Prodigal.v2.11.0-gv;"
            'run_type=Metagenomic;model="A";gc_cont=0.4;'
            "transl_table=11;uses_sd=1\n"
            "viral_contig\tProdigal_v2.11.0-gv\tCDS\t4\t303\t90.0\t+\t0\t"
            "ID=1_1;partial=00;start_type=ATG;genetic_code=11;gc_cont=0.4\n"
            '# Sequence Data: seqnum=2;seqlen=900;seqhdr="alt_code_contig"\n'
            "# Model Data: version=Prodigal.v2.11.0-gv;"
            'run_type=Metagenomic;model="B";gc_cont=0.3;'
            "transl_table=15;uses_sd=0\n"
            "alt_code_contig\tProdigal_v2.11.0-gv\tCDS\t100\t699\t80.0\t-\t0\t"
            "ID=2_1;partial=01;start_type=Edge;genetic_code=15;gc_cont=0.3\n"
        ),
    )
    proteins = _write(
        tmp_path / "genes.faa",
        (
            ">viral_contig_1 # 4 # 303 # 1 # "
            "ID=1_1;partial=00;start_type=ATG;genetic_code=11;gc_cont=0.4\n"
            "MPEPTIDE\n"
            ">alt_code_contig_1 # 100 # 699 # -1 # "
            "ID=2_1;partial=01;start_type=Edge;genetic_code=15;gc_cont=0.3\n"
            "MALTSEQ\n"
        ),
    )

    rows = parse_prodigal_gv(
        gff,
        proteins,
        source_id="assembly1",
        bin_id=None,
    )

    assert rows == [
        {
            "gene_id": "prodigal-gv:assembly1:viral_contig_1",
            "contig_id": "viral_contig",
            "bin_id": None,
            "start": 4,
            "end": 303,
            "strand": "plus",
            "phase": 0,
            "genetic_code": 11,
        },
        {
            "gene_id": "prodigal-gv:assembly1:alt_code_contig_1",
            "contig_id": "alt_code_contig",
            "bin_id": None,
            "start": 100,
            "end": 699,
            "strand": "minus",
            "phase": 0,
            "genetic_code": 15,
        },
    ]


def test_parse_prodigal_gv_rejects_override_and_join_mismatches(
    tmp_path: Path,
) -> None:
    gff = _write(
        tmp_path / "genes.gff",
        """##gff-version  3
# Model Data: version=Prodigal.v2.11.0-gv;transl_table=15;uses_sd=0
contig_1\tProdigal_v2.11.0-gv\tCDS\t1\t300\t.\t+\t0\tID=1_1;genetic_code=15
""",
    )
    proteins = _write(
        tmp_path / "genes.faa",
        ">contig_1_1 # 1 # 300 # 1 # ID=1_1;genetic_code=15\nMPEPTIDE\n",
    )

    with pytest.raises(ValueError, match="uses genetic code 15, expected 11"):
        parse_prodigal_gv(
            gff,
            proteins,
            source_id="assembly1",
            bin_id=None,
            genetic_code=11,
        )

    mismatched = _write(
        tmp_path / "mismatched.faa",
        ">contig_1_1 # 1 # 300 # 1 # ID=1_2;genetic_code=15\nMPEPTIDE\n",
    )
    with pytest.raises(ValueError, match="GFF CDS/protein ID mismatch"):
        parse_prodigal_gv(
            gff,
            mismatched,
            source_id="assembly1",
            bin_id=None,
        )


def test_parse_prodigal_gv_accepts_native_zero_call_outputs(tmp_path: Path) -> None:
    gff = _write(
        tmp_path / "zero.gff",
        (
            "##gff-version  3\n"
            '# Sequence Data: seqnum=1;seqlen=9;seqhdr="short"\n'
            "# Model Data: version=Prodigal.v2.11.0-gv;"
            "run_type=Metagenomic;"
            'model="0|Mycoplasma_bovis_PG45|B|29.3|4|1";'
            "gc_cont=29.30;transl_table=4;uses_sd=1\n"
        ),
    )
    proteins = _write(tmp_path / "zero.faa", "")

    assert (
        parse_prodigal_gv(
            gff,
            proteins,
            source_id="short-virus",
            bin_id=None,
        )
        == []
    )
    with pytest.raises(ValueError, match="uses genetic code 4, expected 11"):
        parse_prodigal_gv(
            gff,
            proteins,
            source_id="short-virus",
            bin_id=None,
            genetic_code=11,
        )
    inconsistent = _write(
        tmp_path / "unexpected.faa",
        ">short_1 # 1 # 9 # 1 # ID=1_1;genetic_code=4\nMK\n",
    )
    with pytest.raises(ValueError, match="proteins exist without GFF CDS"):
        parse_prodigal_gv(
            gff,
            inconsistent,
            source_id="short-virus",
            bin_id=None,
        )


def test_parse_braker3_preserves_isoforms_and_minus_strand_phase(
    tmp_path: Path,
) -> None:
    gff = _write(
        tmp_path / "braker.gff3",
        """##gff-version 3
contig_2\tAUGUSTUS\tgene\t2529\t5201\t.\t-\t.\tID=g1713;
contig_2\tAUGUSTUS\tmRNA\t2529\t5117\t0.64\t-\t.\tID=g1713.t1;Parent=g1713;
contig_2\tAUGUSTUS\tCDS\t2529\t5117\t0.64\t-\t0\tID=g1713.t1.CDS1;Parent=g1713.t1;
contig_2\tgmst\tmRNA\t2529\t5201\t.\t-\t.\tID=g1713.t2;Parent=g1713;
contig_2\tgmst\tCDS\t2529\t5201\t.\t-\t0\tID=g1713.t2.CDS1;Parent=g1713.t2;
contig_11\tAUGUSTUS\tgene\t3406\t4598\t.\t-\t.\tID=g4;
contig_11\tAUGUSTUS\tmRNA\t3406\t4598\t0.55\t-\t.\tID=g4.t1;Parent=g4;
contig_11\tAUGUSTUS\tCDS\t3406\t3782\t0.71\t-\t2\tID=g4.t1.CDS1;Parent=g4.t1;
contig_11\tAUGUSTUS\tCDS\t4388\t4598\t0.77\t-\t0\tID=g4.t1.CDS4;Parent=g4.t1;
""",
    )
    proteins = _write(
        tmp_path / "braker.aa",
        ">g1713.t1\nMISOFORMONE\n>g1713.t2\nMISOFORMTWO\n>g4.t1\nMMINUS\n",
    )

    rows = parse_braker3(
        gff,
        proteins,
        source_id="bin42",
        bin_id="bin42",
    )

    assert len(rows) == 3
    assert rows[0] == {
        "gene_id": "braker3:bin42:bin42:g1713.t1",
        "contig_id": "contig_2",
        "bin_id": "bin42",
        "start": 2529,
        "end": 5117,
        "strand": "minus",
        "phase": 0,
        "genetic_code": 1,
    }
    assert rows[1]["gene_id"] == "braker3:bin42:bin42:g1713.t2"
    assert rows[1]["end"] == 5201
    assert rows[2]["phase"] == 0


def test_parse_braker3_requires_one_protein_per_transcript(tmp_path: Path) -> None:
    gff = _write(
        tmp_path / "braker.gff3",
        """contig_1\tAUGUSTUS\tgene\t1\t300\t.\t+\t.\tID=g1;
contig_1\tAUGUSTUS\tmRNA\t1\t300\t.\t+\t.\tID=g1.t1;Parent=g1;
contig_1\tAUGUSTUS\tCDS\t1\t300\t.\t+\t0\tID=g1.t1.CDS1;Parent=g1.t1;
""",
    )
    proteins = _write(tmp_path / "braker.aa", ">g1.t2\nMPEPTIDE\n")

    with pytest.raises(ValueError, match="transcripts/protein ID mismatch"):
        parse_braker3(
            gff,
            proteins,
            source_id="bin1",
            bin_id="bin1",
        )


def test_merge_annotation_proteins_prefixes_bins_and_writes_mapping(
    tmp_path: Path,
) -> None:
    first_braker = _write(tmp_path / "bin1.faa", ">g1.t1\nMONE*\n")
    second_braker = _write(tmp_path / "bin2.faa", ">g1.t1\nMTWO\n")
    prodigal = _write(
        tmp_path / "viral.faa",
        ">virus_1 # 10 # 210 # 1 # ID=1_1;genetic_code=11\nMVIRUS*\n",
    )
    merged = tmp_path / "annotation.faa"
    mapping = tmp_path / "protein-map.tsv"

    rows = merge_annotation_proteins(
        [
            ProteinSource(first_braker, "braker3", "assembly", "bin1"),
            ProteinSource(second_braker, "braker3", "assembly", "bin2"),
            ProteinSource(prodigal, "prodigal-gv", "assembly", None),
        ],
        merged,
        mapping,
    )

    assert [row.normalized_id for row in rows] == [
        "braker3:assembly:bin1:g1.t1",
        "braker3:assembly:bin2:g1.t1",
        "prodigal-gv:assembly:virus_1",
    ]
    assert [row.native_gene_id for row in rows] == ["g1", "g1", "1_1"]
    assert merged.read_text(encoding="utf-8") == (
        ">braker3:assembly:bin1:g1.t1\nMONE\n"
        ">braker3:assembly:bin2:g1.t1\nMTWO\n"
        ">prodigal-gv:assembly:virus_1\nMVIRUS\n"
    )
    with mapping.open(encoding="utf-8", newline="") as handle:
        mapping_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert mapping_rows[0] == {
        "normalized_id": "braker3:assembly:bin1:g1.t1",
        "native_id": "g1.t1",
        "native_gene_id": "g1",
        "caller": "braker3",
        "source_id": "assembly",
        "bin_id": "bin1",
    }
    assert mapping_rows[2]["bin_id"] == ""
    assert first_braker.read_text(encoding="utf-8") == ">g1.t1\nMONE*\n"
    assert prodigal.read_text(encoding="utf-8").endswith("\nMVIRUS*\n")


@pytest.mark.parametrize("sequence", ["MPEP*TIDE", "MPEPTIDE**", "*"])
def test_merge_annotation_proteins_rejects_invalid_stop_markers(
    tmp_path: Path,
    sequence: str,
) -> None:
    source = _write(tmp_path / "invalid.faa", f">bad.t1\n{sequence}\n")
    output = tmp_path / "merged.faa"
    mapping = tmp_path / "mapping.tsv"

    with pytest.raises(
        ValueError, match="expected at most one terminal stop marker"
    ) as error:
        merge_annotation_proteins(
            [ProteinSource(source, "braker3", "assembly", "bin1")],
            output,
            mapping,
        )

    assert str(source) in str(error.value)
    assert "bad.t1" in str(error.value)
    assert not output.exists()
    assert not mapping.exists()


def test_merge_annotation_proteins_writes_valid_empty_outputs(tmp_path: Path) -> None:
    empty = _write(tmp_path / "empty.faa", "")
    output = tmp_path / "merged.faa"
    mapping = tmp_path / "mapping.tsv"

    rows = merge_annotation_proteins(
        [ProteinSource(empty, "prodigal-gv", "short-virus", None)],
        output,
        mapping,
    )

    assert rows == []
    assert output.read_text(encoding="utf-8") == ""
    assert mapping.read_text(encoding="utf-8") == (
        "normalized_id\tnative_id\tnative_gene_id\tcaller\tsource_id\tbin_id\n"
    )
