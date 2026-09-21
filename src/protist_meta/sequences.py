"""Read FASTA records without loading a whole metagenome into memory."""

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

__all__ = ["Sequence", "read_fasta"]


@dataclass(frozen=True, slots=True)
class Sequence:
    """One DNA or protein record, identified by the first header token."""

    identifier: str
    sequence: str


def read_fasta(path: Path) -> Iterator[Sequence]:
    """Yield nonempty records and reject sequence before the first header."""
    identifier: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if identifier is not None:
                    yield _record(identifier, chunks)
                header = line[1:].split()
                if not header:
                    raise ValueError(f"FASTA header lacks an identifier: {path}")
                identifier = header[0]
                chunks = []
            elif line.strip():
                if identifier is None:
                    raise ValueError(f"FASTA sequence precedes header: {path}")
                chunks.append(line.strip())
    if identifier is not None:
        yield _record(identifier, chunks)


def _record(identifier: str, chunks: list[str]) -> Sequence:
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(f"empty FASTA sequence: {identifier}")
    return Sequence(identifier=identifier, sequence=sequence)
