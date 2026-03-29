"""Imputation strategies for missing CpG sites in methylation matrices.

All functions operate on Polars DataFrames where:
  - Column ``CpGmarker`` holds CpG site IDs
  - Remaining columns are sample-level beta values (Float64)
"""
from __future__ import annotations

import polars as pl


def impute_from_average(
    dnam: pl.DataFrame,
    required_cpgs: list[str] | None = None,
) -> pl.DataFrame:
    """Fill NaN values with the row mean across all samples.

    If *required_cpgs* is given, only rows matching those CpG sites are
    imputed; other rows are left untouched.
    """
    sample_cols = [c for c in dnam.columns if c != "CpGmarker"]
    if not sample_cols:
        return dnam

    row_mean = pl.mean_horizontal(*[pl.col(c) for c in sample_cols])

    if required_cpgs is not None:
        mask = pl.col("CpGmarker").is_in(required_cpgs)
        filled_exprs = [
            pl.when(mask)
            .then(pl.col(c).fill_null(row_mean))
            .otherwise(pl.col(c))
            .alias(c)
            for c in sample_cols
        ]
    else:
        filled_exprs = [
            pl.col(c).fill_null(row_mean).alias(c)
            for c in sample_cols
        ]

    return dnam.with_columns(filled_exprs)


def impute_from_reference(
    dnam: pl.DataFrame,
    reference: pl.DataFrame,
    required_cpgs: list[str] | None = None,
) -> pl.DataFrame:
    """Fill NaN values using a reference table of CpG site mean/median values.

    ``reference`` must have columns ``CpGmarker`` and ``value``.
    Missing CpG rows present in *required_cpgs* but absent from *dnam* are
    added with the reference value for all samples.
    """
    sample_cols = [c for c in dnam.columns if c != "CpGmarker"]
    if not sample_cols:
        return dnam

    ref_map = dict(
        zip(
            reference["CpGmarker"].to_list(),
            reference["value"].to_list(),
        )
    )

    if required_cpgs is not None:
        existing = set(dnam["CpGmarker"].to_list())
        missing = [cpg for cpg in required_cpgs if cpg not in existing]
        if missing:
            missing_rows = pl.DataFrame(
                {
                    "CpGmarker": missing,
                    **{
                        c: [ref_map.get(cpg, None) for cpg in missing]
                        for c in sample_cols
                    },
                }
            ).cast({c: pl.Float64 for c in sample_cols})
            dnam = pl.concat([dnam, missing_rows], how="diagonal_relaxed")

    ref_series = (
        dnam.select("CpGmarker")
        .with_columns(
            pl.col("CpGmarker")
            .replace_strict(ref_map, default=None)
            .cast(pl.Float64)
            .alias("_ref_val")
        )
    )["_ref_val"]

    fill_exprs = [
        pl.col(c).fill_null(ref_series).alias(c)
        for c in sample_cols
    ]
    return dnam.with_columns(fill_exprs)


def hybrid_impute(
    dnam: pl.DataFrame,
    reference: pl.DataFrame,
    required_cpgs: list[str],
    threshold: float = 0.8,
) -> pl.DataFrame:
    """Hybrid imputation: drop low-coverage rows, average-fill the rest,
    then add missing required CpGs from the reference.

    Mirrors biolearn's ``hybrid_impute`` logic.
    """
    sample_cols = [c for c in dnam.columns if c != "CpGmarker"]
    if not sample_cols:
        return dnam

    n_samples = len(sample_cols)
    non_null_count = pl.sum_horizontal(*[pl.col(c).is_not_null().cast(pl.UInt32) for c in sample_cols])
    dnam = dnam.filter(non_null_count >= int(threshold * n_samples))

    dnam = impute_from_average(dnam)

    dnam = impute_from_reference(dnam, reference, required_cpgs)

    return dnam
