"""Utilities for loading coefficient CSV data files shipped with the package."""
from __future__ import annotations

from importlib import resources as importlib_resources
from pathlib import Path

import polars as pl


def get_data_dir() -> Path:
    """Return the path to the package ``data/`` directory."""
    ref = importlib_resources.files("just_biomarkers") / "data"
    return Path(str(ref))


def load_coefficients(filename: str) -> pl.DataFrame:
    """Load a coefficient CSV from the package data directory.

    The CSV is expected to have an index column (CpG site names) as the first
    column, followed by ``CoefficientTraining`` or ``Weight``.
    """
    path = get_data_dir() / filename
    if not path.exists():
        raise FileNotFoundError(
            f"Coefficient file not found: {path}. "
            "Make sure the data files are included in the package."
        )
    df = pl.read_csv(str(path))
    first_col = df.columns[0]
    if first_col != "CpGmarker":
        df = df.rename({first_col: "CpGmarker"})
    return df
