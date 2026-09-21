"""Install the pinned Dfam 4.0 FamDB files with checksum receipts."""

import argparse
import gzip
import hashlib
import logging
import re
import tempfile
import urllib.request
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

LOGGER = logging.getLogger(__name__)
BASE_URL = "https://www.dfam.org/releases/Dfam_4.0/families/FamDB"
BUFFER_BYTES = 1024 * 1024
MD5_PATTERN = re.compile(r"^([0-9a-fA-F]{32})(?:\s+\*?\S+)?$")


@dataclass(frozen=True, slots=True)
class ReleaseFile:
    """One compressed release file and its published size."""

    compressed_name: str
    compressed_bytes: int
    publisher_md5: str

    @property
    def name(self) -> str:
        """Return the decompressed FamDB filename."""
        return self.compressed_name.removesuffix(".gz")

    @property
    def url(self) -> str:
        """Return the release-specific download URL."""
        return f"{BASE_URL}/{self.compressed_name}"


RELEASE_FILES = (
    ReleaseFile(
        "dfam40.0.h5.gz",
        60_708_386,
        "234d177775f1bf3445b1fe146bc6e65e",
    ),
    ReleaseFile(
        "dfam40.curated.consensus.0.h5.gz",
        28_260_067,
        "7892e18016fc820264e625cbb9ec607b",
    ),
)


def install_release(destination: Path) -> None:
    """Download, verify, and atomically install the pinned release subset."""
    destination = destination.resolve()
    if destination.exists():
        raise FileExistsError(f"Dfam destination already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".dfam-4.0-", dir=destination.parent
    ) as temporary:
        work_dir = Path(temporary)
        install_dir = work_dir / "famdb"
        downloads_dir = work_dir / "downloads"
        install_dir.mkdir()
        downloads_dir.mkdir()
        receipts = [
            _install_file(release_file, downloads_dir, install_dir)
            for release_file in RELEASE_FILES
        ]
        _write_receipts(install_dir, receipts)
        install_dir.replace(destination)
    LOGGER.info("installed Dfam 4.0 FamDB files at %s", destination)


def _install_file(
    release_file: ReleaseFile, downloads_dir: Path, install_dir: Path
) -> tuple[ReleaseFile, str, str, str, int]:
    publisher_md5 = _read_publisher_md5(f"{release_file.url}.md5")
    if publisher_md5 != release_file.publisher_md5:
        raise ValueError(
            f"publisher MD5 changed for {release_file.compressed_name}: "
            f"{publisher_md5} != {release_file.publisher_md5}"
        )
    compressed_path = downloads_dir / release_file.compressed_name
    compressed_bytes, observed_md5, compressed_sha256 = _download(
        release_file.url, compressed_path
    )
    if compressed_bytes != release_file.compressed_bytes:
        raise ValueError(
            f"unexpected size for {release_file.compressed_name}: "
            f"{compressed_bytes} != {release_file.compressed_bytes}"
        )
    if observed_md5 != release_file.publisher_md5:
        raise ValueError(
            f"publisher MD5 mismatch for {release_file.compressed_name}: "
            f"{observed_md5} != {release_file.publisher_md5}"
        )
    decompressed_sha256, decompressed_bytes = _decompress(
        compressed_path, install_dir / release_file.name
    )
    return (
        release_file,
        publisher_md5,
        compressed_sha256,
        decompressed_sha256,
        decompressed_bytes,
    )


def _read_publisher_md5(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "protist-meta-nf"})
    with urllib.request.urlopen(request, timeout=60) as response:
        line = response.read(1024).decode("ascii").strip()
    match = MD5_PATTERN.fullmatch(line)
    if match is None:
        raise ValueError(f"invalid publisher MD5 response from {url}")
    return match.group(1).lower()


def _download(url: str, destination: Path) -> tuple[int, str, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "protist-meta-nf"})
    md5 = hashlib.md5()  # Dfam publishes MD5 for transport checks.
    sha256 = hashlib.sha256()
    size = 0
    with (
        urllib.request.urlopen(request, timeout=60) as response,
        destination.open("xb") as handle,
    ):
        while chunk := response.read(BUFFER_BYTES):
            handle.write(chunk)
            md5.update(chunk)
            sha256.update(chunk)
            size += len(chunk)
    return size, md5.hexdigest(), sha256.hexdigest()


def _decompress(source: Path, destination: Path) -> tuple[str, int]:
    sha256 = hashlib.sha256()
    size = 0
    with gzip.open(source, "rb") as compressed, destination.open("xb") as handle:
        while chunk := compressed.read(BUFFER_BYTES):
            handle.write(chunk)
            sha256.update(chunk)
            size += len(chunk)
    if size == 0:
        raise ValueError(f"decompressed Dfam file is empty: {source}")
    return sha256.hexdigest(), size


def _write_receipts(
    destination: Path,
    receipts: Sequence[tuple[ReleaseFile, str, str, str, int]],
) -> None:
    publisher_lines: list[str] = []
    sha256_lines: list[str] = []
    source_header = (
        "file\turl\tpublisher_md5_url\tpublisher_md5\tcompressed_bytes\t"
        "compressed_sha256\tbytes\tsha256"
    )
    source_lines = [source_header]
    for release_file, publisher_md5, compressed_sha256, sha256, size in receipts:
        publisher_lines.append(f"{publisher_md5}  {release_file.compressed_name}")
        sha256_lines.append(f"{sha256}  {release_file.name}")
        source_lines.append(
            "\t".join(
                (
                    release_file.name,
                    release_file.url,
                    f"{release_file.url}.md5",
                    publisher_md5,
                    str(release_file.compressed_bytes),
                    compressed_sha256,
                    str(size),
                    sha256,
                )
            )
        )
    (destination / "PUBLISHER_MD5SUMS").write_text(
        "\n".join(publisher_lines) + "\n", encoding="utf-8"
    )
    (destination / "SHA256SUMS").write_text(
        "\n".join(sha256_lines) + "\n", encoding="utf-8"
    )
    (destination / "sources.tsv").write_text(
        "\n".join(source_lines) + "\n", encoding="utf-8"
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse the fixed-release installer arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "destination",
        type=Path,
        help="new directory that will contain the Dfam 4.0 FamDB files",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Install the pinned release subset and return the process exit code."""
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    install_release(args.destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
