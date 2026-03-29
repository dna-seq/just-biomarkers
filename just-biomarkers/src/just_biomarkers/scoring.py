"""Core scoring engine for linear methylation clocks (Polars implementation).

This is the primary compute path.  Given a methylation matrix and a clock
definition, it performs:
  1. Load coefficient table
  2. Join coefficients with the methylation matrix on CpG site
  3. Add intercept row
  4. Weighted sum per sample
  5. Apply optional transform
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import polars as pl

from just_biomarkers.data_utils import load_coefficients
from just_biomarkers.models import BatchClockResult, ClockDefinition, ClockResult
from just_biomarkers.registry import CLOCK_DEFINITIONS, get_clock
from just_biomarkers.transforms import TRANSFORM_REGISTRY, identity


def _resolve_transform(clock: ClockDefinition) -> Callable:
    if clock.transform_name is None:
        return identity
    if clock.transform_name not in TRANSFORM_REGISTRY:
        raise ValueError(
            f"Unknown transform '{clock.transform_name}' for clock '{clock.name}'"
        )
    return TRANSFORM_REGISTRY[clock.transform_name]


def _resolve_preprocess(clock: ClockDefinition) -> Callable[[pl.DataFrame], pl.DataFrame] | None:
    if clock.preprocess_name is None:
        return None
    from just_biomarkers.preprocessing import PREPROCESS_REGISTRY
    if clock.preprocess_name not in PREPROCESS_REGISTRY:
        raise ValueError(
            f"Unknown preprocess '{clock.preprocess_name}' for clock '{clock.name}'"
        )
    return PREPROCESS_REGISTRY[clock.preprocess_name]


def score_clock(
    dnam: pl.DataFrame,
    clock: ClockDefinition | str,
) -> pl.DataFrame:
    """Score a single clock against a methylation matrix.

    Parameters
    ----------
    dnam
        Methylation DataFrame with ``CpGmarker`` column and one float column
        per sample.
    clock
        Either a :class:`ClockDefinition` instance or the string name of a
        registered clock.

    Returns
    -------
    pl.DataFrame
        Columns: ``sample_id``, ``clock_name``, ``score``,
        ``cpgs_matched``, ``cpgs_required``.
    """
    if isinstance(clock, str):
        clock = get_clock(clock)

    preprocess = _resolve_preprocess(clock)
    dnam_to_score = preprocess(dnam) if preprocess is not None else dnam

    coeff_df = load_coefficients(clock.coefficient_file)

    coeff_col = clock.coefficient_column
    if coeff_col not in coeff_df.columns:
        alt = "Weight"
        if alt in coeff_df.columns:
            coeff_col = alt
        else:
            raise ValueError(
                f"Coefficient column '{clock.coefficient_column}' not found "
                f"in {clock.coefficient_file}. Columns: {coeff_df.columns}"
            )

    coeff_only = coeff_df.select(["CpGmarker", coeff_col])

    sample_cols = [c for c in dnam_to_score.columns if c != "CpGmarker"]
    required_cpgs = set(coeff_only["CpGmarker"].to_list()) - {clock.intercept_name}
    n_required = len(required_cpgs)

    joined = dnam_to_score.join(coeff_only, on="CpGmarker", how="inner")

    matched_cpgs = set(joined["CpGmarker"].to_list()) - {clock.intercept_name}
    n_matched = len(matched_cpgs)

    intercept_rows = joined.filter(pl.col("CpGmarker") == clock.intercept_name)
    if intercept_rows.height == 0:
        intercept_row_in_coeff = coeff_only.filter(
            pl.col("CpGmarker") == clock.intercept_name
        )
        if intercept_row_in_coeff.height > 0:
            intercept_val = intercept_row_in_coeff[coeff_col][0]
            new_row = pl.DataFrame(
                {
                    "CpGmarker": [clock.intercept_name],
                    coeff_col: [intercept_val],
                    **{c: [1.0] for c in sample_cols},
                }
            )
            joined = pl.concat(
                [joined, new_row.select(joined.columns)], how="diagonal_relaxed"
            )

    transform = _resolve_transform(clock)

    results_rows: list[dict] = []
    for col in sample_cols:
        weighted = joined.select(
            (pl.col(col) * pl.col(coeff_col)).alias("weighted")
        )
        raw_sum = weighted["weighted"].sum()
        score = float(transform(raw_sum))
        results_rows.append(
            {
                "sample_id": col,
                "clock_name": clock.name,
                "score": score,
                "cpgs_matched": n_matched,
                "cpgs_required": n_required,
            }
        )

    return pl.DataFrame(results_rows)


def score_clocks(
    dnam: pl.DataFrame,
    clock_names: list[str] | None = None,
) -> BatchClockResult:
    """Score multiple clocks against a methylation matrix.

    Parameters
    ----------
    dnam
        Methylation DataFrame (``CpGmarker`` + sample columns).
    clock_names
        List of clock names to score.  If ``None``, all registered linear
        clocks are used.

    Returns
    -------
    BatchClockResult
    """
    if clock_names is None:
        clock_names = sorted(CLOCK_DEFINITIONS.keys())

    all_results: list[ClockResult] = []
    warnings: list[str] = []

    for name in clock_names:
        clock = get_clock(name)
        result_df = score_clock(dnam, clock)
        for row in result_df.iter_rows(named=True):
            match_rate = row["cpgs_matched"] / max(row["cpgs_required"], 1)
            all_results.append(
                ClockResult(
                    clock_name=row["clock_name"],
                    sample_id=row["sample_id"],
                    score=row["score"],
                    output=clock.output,
                    cpgs_matched=row["cpgs_matched"],
                    cpgs_required=row["cpgs_required"],
                    match_rate=match_rate,
                )
            )
            if match_rate < 0.8:
                warnings.append(
                    f"Low CpG match rate for {name} on sample {row['sample_id']}: "
                    f"{row['cpgs_matched']}/{row['cpgs_required']} ({match_rate:.1%})"
                )

    return BatchClockResult(results=all_results, warnings=warnings)
