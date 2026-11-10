"""Utilities for loading clock coefficient CSV files.

Files come from the package ``data/`` directory when bundled, else from the
verified download cache (see :mod:`just_biomarkers.coefficients`).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from just_biomarkers.coefficients import coefficient_path, package_data_dir


def get_data_dir() -> Path:
    """Return the path to the package ``data/`` directory (may not exist)."""
    return package_data_dir()


def load_coefficients(filename: str) -> pl.DataFrame:
    """Load a coefficient CSV, downloading it from biolearn if not available locally.

    The CSV is expected to have an index column (CpG site names) as the first
    column, followed by ``CoefficientTraining`` or ``Weight``.
    """
    df = pl.read_csv(str(coefficient_path(filename)))
    first_col = df.columns[0]
    if first_col != "CpGmarker":
        df = df.rename({first_col: "CpGmarker"})
    return df
