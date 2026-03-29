"""Download ComputAgeBench datasets from Hugging Face for clock benchmarking.

ComputAgeBench is an epigenetic aging clocks benchmark dataset hosted at
``computage/computage_bench`` on Hugging Face.  It contains 65 GEO studies
with ~10k samples and ~900k CpG features, along with metadata (age, condition,
tissue, etc.).

The parquet files follow the same orientation as our standard: rows = CpG sites,
columns = sample IDs.  Metadata is a TSV with per-sample annotations.

Reference
---------
Kriukov et al., "ComputAgeBench: Epigenetic Aging Clocks Benchmark", bioRxiv 2024.
https://doi.org/10.1101/2024.06.06.597715
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import platformdirs
import polars as pl

_HF_REPO_ID = "computage/computage_bench"
_CACHE_APP_NAME = "just-biomarkers"


def _cache_dir() -> Path:
    """Return the platform-appropriate cache directory for ComputAgeBench."""
    base = Path(platformdirs.user_cache_dir(_CACHE_APP_NAME))
    cab_dir = base / "computage_bench"
    cab_dir.mkdir(parents=True, exist_ok=True)
    return cab_dir


def _ensure_snapshot(force: bool = False) -> Path:
    """Download (or locate cached) ComputAgeBench HF snapshot.

    Uses ``huggingface_hub.snapshot_download`` which handles caching,
    resumable downloads, and LFS file resolution automatically.

    Returns the local directory path of the downloaded snapshot.
    """
    from huggingface_hub import snapshot_download

    local_dir = _cache_dir() / "repo"
    return Path(
        snapshot_download(
            repo_id=_HF_REPO_ID,
            repo_type="dataset",
            local_dir=str(local_dir),
            force_download=force,
        )
    )


def list_computage_datasets(split: str = "benchmark") -> list[str]:
    """List available ComputAgeBench dataset IDs (GEO accessions).

    Parameters
    ----------
    split:
        Either ``"benchmark"`` or ``"train"``.

    Returns
    -------
    list[str]
        Sorted list of GEO dataset IDs available in the split.
    """
    snap = _ensure_snapshot()
    data_dir = snap / "data" / split
    if not data_dir.exists():
        raise FileNotFoundError(
            f"ComputAgeBench split '{split}' not found at {data_dir}. "
            "Try downloading with force=True."
        )
    prefix = f"computage_{split}_data_"
    ids = sorted(
        p.stem.replace(prefix, "")
        for p in data_dir.glob("*.parquet")
    )
    return ids


def load_computage_meta(split: str = "benchmark") -> pl.DataFrame:
    """Load the ComputAgeBench metadata table.

    Parameters
    ----------
    split:
        Either ``"benchmark"`` or ``"train"``.

    Returns
    -------
    pl.DataFrame
        Metadata with columns: sample_id, DatasetID, PlatformID, Tissue,
        CellType, Gender, Age, Condition, Class.
    """
    snap = _ensure_snapshot()
    meta_file = snap / f"computage_{split}_meta.tsv"
    if not meta_file.exists():
        raise FileNotFoundError(
            f"ComputAgeBench metadata file not found: {meta_file}"
        )
    df = pl.read_csv(str(meta_file), separator="\t")
    first_col = df.columns[0]
    if first_col != "sample_id":
        df = df.rename({first_col: "sample_id"})
    return df


def load_computage_dataset(
    dataset_id: str,
    *,
    split: str = "benchmark",
    max_samples: Optional[int] = None,
) -> pl.DataFrame:
    """Load a single ComputAgeBench parquet file as a CpGmarker DataFrame.

    The returned DataFrame follows the just-biomarkers standard:
    first column is ``CpGmarker`` (CpG site IDs), remaining columns are
    sample IDs with Float64 beta values.

    Parameters
    ----------
    dataset_id:
        GEO accession, e.g. ``"GSE100264"``.
    split:
        Either ``"benchmark"`` or ``"train"``.
    max_samples:
        If given, keep only the first *N* sample columns.

    Returns
    -------
    pl.DataFrame
        Methylation matrix ready for ``score_clocks()``.
    """
    snap = _ensure_snapshot()
    parquet_path = (
        snap / "data" / split / f"computage_{split}_data_{dataset_id}.parquet"
    )
    if not parquet_path.exists():
        available = list_computage_datasets(split=split)
        raise FileNotFoundError(
            f"Dataset '{dataset_id}' not found in ComputAgeBench {split} split. "
            f"Available: {', '.join(available[:10])}{'...' if len(available) > 10 else ''}"
        )

    df = pl.read_parquet(str(parquet_path))

    first_col = df.columns[0]
    if first_col != "CpGmarker":
        df = df.rename({first_col: "CpGmarker"})

    sample_cols = [c for c in df.columns if c != "CpGmarker"]
    if max_samples is not None and max_samples < len(sample_cols):
        sample_cols = sample_cols[:max_samples]
        df = df.select(["CpGmarker"] + sample_cols)

    df = df.with_columns(
        [pl.col(c).cast(pl.Float64, strict=False) for c in sample_cols]
    )
    return df


def download_computage_dataset(
    dataset_id: str,
    output_dir: Path | str,
    *,
    split: str = "benchmark",
    max_samples: Optional[int] = None,
    force: bool = False,
) -> Path:
    """Download a ComputAgeBench dataset and save as CpGmarker CSV.

    Parameters
    ----------
    dataset_id:
        GEO accession, e.g. ``"GSE100264"``.
    output_dir:
        Directory to write the output CSV into.
    split:
        Either ``"benchmark"`` or ``"train"``.
    max_samples:
        If given, keep only the first *N* sample columns.
    force:
        Re-download and re-convert even if the output CSV already exists.

    Returns
    -------
    Path
        The saved CSV file path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"_max{max_samples}" if max_samples is not None else ""
    dest_csv = output_dir / f"computage_{split}_{dataset_id}{suffix}_methylation.csv"

    if dest_csv.exists() and not force:
        return dest_csv

    df = load_computage_dataset(
        dataset_id, split=split, max_samples=max_samples
    )
    df.write_csv(str(dest_csv))
    return dest_csv


def download_computage_meta(
    output_dir: Path | str,
    *,
    split: str = "benchmark",
    force: bool = False,
) -> Path:
    """Download ComputAgeBench metadata and save as TSV.

    Parameters
    ----------
    output_dir:
        Directory to write the metadata TSV.
    split:
        Either ``"benchmark"`` or ``"train"``.
    force:
        Re-download even if the output file already exists.

    Returns
    -------
    Path
        The saved TSV file path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    dest_tsv = output_dir / f"computage_{split}_meta.tsv"
    if dest_tsv.exists() and not force:
        return dest_tsv

    meta = load_computage_meta(split=split)
    meta.write_csv(str(dest_tsv), separator="\t")
    return dest_tsv
