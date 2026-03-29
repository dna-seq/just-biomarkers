"""IO helpers for loading and validating methylation matrices.

Methylation standard (following biolearn):
  - Rows  = CpG site IDs (e.g. cg00075967)
  - Columns = sample IDs
  - Values  = beta values in [0, 1] or null/NaN
"""
from __future__ import annotations

from pathlib import Path
from typing import Literal

import polars as pl

IdatPreprocess = Literal["raw", "noob"]


def read_methylation_csv(
    path: Path | str,
    *,
    separator: str = ",",
) -> pl.DataFrame:
    """Read a methylation beta-value matrix from CSV.

    The first column is assumed to contain CpG site IDs and is renamed to
    ``CpGmarker``.  All other columns are sample IDs with float values.
    """
    path = Path(path)
    df = pl.read_csv(str(path), separator=separator, infer_schema_length=10_000)

    first_col = df.columns[0]
    if first_col != "CpGmarker":
        df = df.rename({first_col: "CpGmarker"})

    sample_cols = [c for c in df.columns if c != "CpGmarker"]
    df = df.with_columns(
        [pl.col(c).cast(pl.Float64, strict=False) for c in sample_cols]
    )
    return df


def read_methylation_matrix(
    path: Path | str,
    *,
    idat_preprocess: IdatPreprocess = "raw",
) -> pl.DataFrame:
    """Read methylation matrix from IDAT/CSV/TSV/TXT/Parquet with auto-detection.

    Notes
    -----
    Raw Illumina ``*.idat`` inputs require ``pylluminator-modern``.
    """
    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Input path does not exist: {resolved}")

    if resolved.is_dir():
        idat_files = sorted(resolved.glob("*.idat"))
        if idat_files:
            return read_methylation_idat(resolved, preprocess=idat_preprocess)
        raise ValueError(
            "Directory input is not supported unless it contains a pre-merged "
            "tabular matrix file. Provide a CSV/TSV/TXT/Parquet matrix path."
        )

    suffixes = [s.lower() for s in resolved.suffixes]
    name_lower = resolved.name.lower()
    if name_lower.endswith(".idat"):
        return read_methylation_idat(resolved.parent, preprocess=idat_preprocess)

    if suffixes and suffixes[-1] == ".parquet":
        return read_methylation_parquet(resolved)

    if suffixes and suffixes[-1] in {".csv", ".tsv", ".txt"}:
        separator = "\t" if suffixes[-1] == ".tsv" else ","
        return read_methylation_csv(resolved, separator=separator)

    if (
        len(suffixes) >= 2
        and suffixes[-1] == ".gz"
        and suffixes[-2] in {".csv", ".tsv", ".txt"}
    ):
        separator = "\t" if suffixes[-2] == ".tsv" else ","
        return read_methylation_csv(resolved, separator=separator)

    raise ValueError(
        "Unsupported methylation input format. "
        "Supported: .idat (directory or file pair), .csv, .tsv, .txt, "
        ".csv.gz, .tsv.gz, .txt.gz, .parquet. "
        f"Got: {resolved.name}"
    )


def read_methylation_idat(
    path: Path | str,
    *,
    preprocess: IdatPreprocess = "raw",
    include_out_of_band: bool = True,
    apply_mask: bool = True,
    strip_probe_suffix: bool = True,
) -> pl.DataFrame:
    """Read Illumina IDAT files and return a CpG x sample beta matrix.

    Requires the optional dependency ``pylluminator-modern`` (module name:
    ``pylluminator``).
    """
    try:
        from pylluminator.samples import read_samples  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "IDAT support requires 'pylluminator-modern'. "
            "Install with: uv add \"pylluminator-modern; python_version >= '3.12'\""
        ) from exc

    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"IDAT input path does not exist: {resolved}")

    samples = read_samples(resolved, annotation=None, keep_idat=False)
    if samples is None:
        raise ValueError(f"Could not parse IDAT input at: {resolved}")

    if preprocess == "noob":
        try:
            samples.infer_type1_channel()
            samples.dye_bias_correction_nl()
            samples.poobah()
            samples.noob_background_correction()
        except Exception as exc:
            raise ValueError(
                "IDAT 'noob' preprocessing failed. "
                "Try preprocess='raw' for this input. "
                f"Original error: {exc}"
            ) from exc
    elif preprocess != "raw":
        raise ValueError(
            f"Unsupported IDAT preprocess mode: {preprocess}. "
            "Supported: 'raw', 'noob'."
        )

    samples.calculate_betas(include_out_of_band=include_out_of_band)
    betas = samples.get_betas(apply_mask=apply_mask)
    if betas is None:
        raise ValueError("Failed to compute beta values from IDAT input.")

    index = betas.index
    probe_ids = (
        index.get_level_values("probe_id")
        if getattr(index, "names", None) and "probe_id" in index.names
        else index.get_level_values(0)
        if getattr(index, "nlevels", 1) > 1
        else index
    )

    betas = betas.copy()
    betas.insert(0, "CpGmarker", probe_ids.astype(str))
    if strip_probe_suffix:
        betas["CpGmarker"] = betas["CpGmarker"].str.split("_", n=1).str[0]

    # Multiple probe variants can map to the same cg ID after suffix stripping.
    betas = betas.groupby("CpGmarker", as_index=False).mean(numeric_only=True)

    df = pl.from_pandas(betas)
    sample_cols = [c for c in df.columns if c != "CpGmarker"]
    return df.with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in sample_cols])


def read_methylation_parquet(path: Path | str) -> pl.DataFrame:
    """Read a methylation matrix stored as parquet."""
    path = Path(path)
    df = pl.read_parquet(str(path))
    first_col = df.columns[0]
    if first_col != "CpGmarker":
        df = df.rename({first_col: "CpGmarker"})
    return df


def validate_methylation_matrix(df: pl.DataFrame) -> list[str]:
    """Validate a methylation DataFrame, return list of warnings (empty if OK)."""
    warnings: list[str] = []
    if "CpGmarker" not in df.columns:
        warnings.append("Missing 'CpGmarker' column")
        return warnings

    sample_cols = [c for c in df.columns if c != "CpGmarker"]
    if not sample_cols:
        warnings.append("No sample columns found")
        return warnings

    for col in sample_cols:
        if df[col].dtype != pl.Float64:
            warnings.append(f"Column '{col}' is not Float64 (got {df[col].dtype})")

    n_rows = df.height
    if n_rows == 0:
        warnings.append("Empty methylation matrix")
    return warnings
