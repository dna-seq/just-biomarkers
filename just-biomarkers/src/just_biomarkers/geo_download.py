"""Download GEO methylation example datasets without a biolearn dependency.

Protocol
--------
GEO series matrix files are gzip-compressed tab-separated text files hosted
on the NCBI FTP server.  Their structure is:

  <header lines starting with "!">
  !series_matrix_table_begin
  "ID_REF"   "GSM..."  "GSM..."  ...     <- column header
  "cg..."    0.812     0.345     ...     <- beta-value rows
  ...
  !series_matrix_table_end

We skip everything above ``!series_matrix_table_begin``, parse the tab table
with Polars, rename ``ID_REF`` → ``CpGmarker``, and write a CSV.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import shutil
import urllib.request
from pathlib import Path
from typing import Optional

import platformdirs
import polars as pl

# ---------------------------------------------------------------------------
# Curated catalogue of well-tested GEO datasets usable as example inputs.
# Each entry: dataset_id → (ftp_url, human description, matrix_start_keyword)
# The FTP URLs are the same ones used by biolearn's library.yaml.
# ---------------------------------------------------------------------------

GEO_EXAMPLE_DATASETS: dict[str, dict[str, str]] = {
    "GSE112618": {
        "description": "FACS validation – IlluminaEPIC, 6 blood samples",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112618/matrix/GSE112618_series_matrix.txt.gz",
        "format": "IlluminaEPIC",
        "samples": "6",
        "tissue": "Blood",
    },
    "GSE110554": {
        "description": "FlowSorted.Blood.EPIC – 49 blood cell-type reference samples",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110554/matrix/GSE110554_series_matrix.txt.gz",
        "format": "IlluminaEPIC",
        "samples": "49",
        "tissue": "Blood",
    },
    "GSE40279": {
        "description": "Genome-wide methylation profiles across age – 656 blood samples (450k)",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz",
        "format": "Illumina450k",
        "samples": "656",
        "tissue": "Blood",
    },
    "GSE41169": {
        "description": "Blood DNA methylation – 95 samples, schizophrenia study (450k)",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41169/matrix/GSE41169_series_matrix.txt.gz",
        "format": "Illumina450k",
        "samples": "95",
        "tissue": "Blood",
    },
    "GSE164056": {
        "description": "Social Anxiety Disorder – 143 blood samples (EPIC)",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164056/matrix/GSE164056_series_matrix.txt.gz",
        "format": "IlluminaEPIC",
        "samples": "143",
        "tissue": "Blood",
    },
}

_TABLE_BEGIN = b"!series_matrix_table_begin"
_TABLE_END = "!series_matrix_table_end"
_CACHE_APP_NAME = "just-biomarkers"


def _cache_dir() -> Path:
    """Return the platform-appropriate cache directory for raw downloads."""
    base = Path(platformdirs.user_cache_dir(_CACHE_APP_NAME))
    geo_dir = base / "geo_downloads"
    geo_dir.mkdir(parents=True, exist_ok=True)
    return geo_dir


def _cached_gz_path(url: str) -> Path:
    """Deterministic local cache path for a URL (sha256 of url + original ext)."""
    digest = hashlib.sha256(url.encode()).hexdigest()
    ext = Path(url).suffix  # e.g. ".gz"
    return _cache_dir() / f"{digest}{ext}"


def _download_gz(url: str, dest: Path) -> None:
    """Stream-download *url* to *dest*, showing a simple progress indicator."""
    with urllib.request.urlopen(url) as response:  # noqa: S310
        with dest.open("wb") as out_file:
            shutil.copyfileobj(response, out_file)


def _read_matrix_from_gz(gz_path: Path, max_samples: Optional[int]) -> pl.DataFrame:
    """Parse a GEO series matrix .gz and return a CpGmarker-first Polars DataFrame.

    The function locates the ``!series_matrix_table_begin`` sentinel, then
    reads the tab-separated beta-value table that follows it.  The trailing
    ``!series_matrix_table_end`` sentinel row is dropped.
    """
    with gzip.open(gz_path, "rb") as fh:
        # Scan for the table-begin sentinel, accumulate lines after it.
        in_table = False
        table_lines: list[bytes] = []
        for raw_line in fh:
            if not in_table:
                if raw_line.strip().startswith(_TABLE_BEGIN):
                    in_table = True
                continue
            # We are inside the table
            stripped = raw_line.strip().decode("utf-8", errors="replace")
            if stripped.startswith(_TABLE_END):
                break
            table_lines.append(raw_line)

    if not table_lines:
        raise ValueError(
            f"No matrix table found in {gz_path.name}. "
            "The series matrix may use an unsupported format."
        )

    raw_bytes = b"".join(table_lines)
    df = pl.read_csv(
        io.BytesIO(raw_bytes),
        separator="\t",
        infer_schema_length=500,
        null_values=["", "NA", "NaN", "null"],
    )

    # Rename the CpG-ID column to our standard name
    first_col = df.columns[0]
    if first_col != "CpGmarker":
        df = df.rename({first_col: "CpGmarker"})

    # Strip surrounding quotes that GEO sometimes leaves on the CpGmarker values
    df = df.with_columns(
        pl.col("CpGmarker").str.strip_chars('"').str.strip_chars("'")
    )

    sample_cols = [c for c in df.columns if c != "CpGmarker"]
    if max_samples is not None and max_samples < len(sample_cols):
        sample_cols = sample_cols[:max_samples]
        df = df.select(["CpGmarker"] + sample_cols)

    # Ensure all sample columns are Float64
    df = df.with_columns(
        [pl.col(c).cast(pl.Float64, strict=False) for c in sample_cols]
    )
    return df


def download_geo_example(
    dataset_id: str,
    output_dir: Path | str,
    *,
    max_samples: Optional[int] = None,
    force: bool = False,
) -> Path:
    """Download a GEO methylation dataset and save it as a CpGmarker CSV.

    The raw ``.txt.gz`` file is cached under the OS user-cache directory so
    repeated calls (or CLI re-runs) are instant after the first download.

    Parameters
    ----------
    dataset_id:
        A key from :data:`GEO_EXAMPLE_DATASETS`, e.g. ``"GSE112618"``.
    output_dir:
        Directory to write the output CSV into.
    max_samples:
        If given, keep only the first *N* sample columns (useful for the very
        large 450k datasets like GSE40279 with 656 samples).
    force:
        Re-download and re-parse even if the output CSV already exists.

    Returns
    -------
    Path
        The saved CSV file path, ready for :func:`just_biomarkers.io.read_methylation_csv`.
    """
    if dataset_id not in GEO_EXAMPLE_DATASETS:
        known = ", ".join(GEO_EXAMPLE_DATASETS)
        raise ValueError(
            f"Unknown dataset '{dataset_id}'. Known datasets: {known}"
        )

    info = GEO_EXAMPLE_DATASETS[dataset_id]
    url: str = info["url"]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"_max{max_samples}" if max_samples is not None else ""
    dest_csv = output_dir / f"{dataset_id}{suffix}_methylation.csv"

    if dest_csv.exists() and not force:
        return dest_csv

    gz_path = _cached_gz_path(url)
    if not gz_path.exists() or force:
        _download_gz(url, gz_path)

    df = _read_matrix_from_gz(gz_path, max_samples=max_samples)
    df.write_csv(str(dest_csv))
    return dest_csv
