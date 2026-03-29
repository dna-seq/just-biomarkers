"""Pre-scoring preprocessing functions for methylation clocks.

Some clocks require the beta-value matrix to be transformed *before* the
standard weighted-sum step.  Currently the only such clock is DunedinPACE,
which applies quantile normalization against 20,000 gold-standard probe means.

Preprocessing functions have the signature:
    preprocess(dnam: pl.DataFrame) -> pl.DataFrame
where *dnam* has a ``CpGmarker`` column + one float column per sample.
The returned DataFrame has the same structure but potentially different row
count (only probes covered by the preprocessing reference are kept).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import polars as pl

from just_biomarkers.data_utils import get_data_dir, load_coefficients


def _rankdata_average(a: np.ndarray) -> np.ndarray:
    """Average-rank of array elements — pure NumPy equivalent of
    ``scipy.stats.rankdata(a, method='average')``.

    Rank starts at 1.  Tied elements receive the average of the ranks they
    would occupy.
    """
    n = len(a)
    sorter = np.argsort(a, kind="stable")
    a_sorted = a[sorter]

    # Locate boundaries of tie groups
    obs = np.concatenate(([True], a_sorted[1:] != a_sorted[:-1]))
    # 1-based group index for each element (in sorted order)
    dense = obs.cumsum()
    # Positions where each new group starts, plus sentinel = n
    count = np.concatenate((np.where(obs)[0], [n]))
    # Average rank for each group: (first_rank + last_rank) / 2
    avg_ranks = 0.5 * (count[dense - 1] + count[dense]) + 0.5

    result = np.empty(n, dtype=np.float64)
    result[sorter] = avg_ranks
    return result


def _quantile_normalize_column(
    column: np.ndarray,
    sorted_target: np.ndarray,
) -> np.ndarray:
    """Quantile-normalize one sample column against a sorted target distribution.

    Replicates the per-column logic of biolearn's
    ``quantile_normalize_using_target``.

    Parameters
    ----------
    column:
        1-D float array of beta values for one sample (length = n_probes).
        The array is modified **in place** and returned.
    sorted_target:
        Sorted gold-standard mean values (ascending).

    Returns
    -------
    np.ndarray
        The normalized column (same object as *column*).
    """
    ranks = _rankdata_average(column)
    floor_ranks = np.floor(ranks).astype(int)
    has_decimal_above_0_4 = (ranks - floor_ranks) > 0.4

    # Ties / half-integer ranks: interpolate between adjacent target values
    column[has_decimal_above_0_4] = 0.5 * (
        sorted_target[floor_ranks[has_decimal_above_0_4] - 1]
        + sorted_target[floor_ranks[has_decimal_above_0_4]]
    )
    # Exact integer ranks: pick directly from sorted target
    column[~has_decimal_above_0_4] = sorted_target[
        floor_ranks[~has_decimal_above_0_4] - 1
    ]
    return column


def _dunedin_pace_preprocess(dnam: pl.DataFrame) -> pl.DataFrame:
    """Apply DunedinPACE quantile normalization to a methylation matrix.

    Algorithm (mirrors ``dunedin_pace_normalization`` from biolearn):

    1. Load the 20,000-probe gold-standard means shipped with the package.
    2. Restrict *dnam* to the intersection of its probes and the gold standard.
    3. For missing gold-standard probes, impute with their gold-standard mean.
    4. Quantile-normalize each sample column to match the gold-standard mean
       distribution.

    Parameters
    ----------
    dnam:
        Full methylation matrix (``CpGmarker`` + sample columns).

    Returns
    -------
    pl.DataFrame
        Quantile-normalized matrix containing only gold-standard probes
        (≤20,000 rows), same sample columns.
    """
    gold_path = get_data_dir() / "DunedinPACE_Gold_Means.csv"
    gold_df = pl.read_csv(str(gold_path))
    first_col = gold_df.columns[0]
    if first_col != "CpGmarker":
        gold_df = gold_df.rename({first_col: "CpGmarker"})

    gold_means: dict[str, float] = dict(
        zip(gold_df["CpGmarker"].to_list(), gold_df["mean"].to_list())
    )
    gold_probes = list(gold_means.keys())
    sorted_target = np.sort(list(gold_means.values()))

    sample_cols = [c for c in dnam.columns if c != "CpGmarker"]

    # --- Filter dnam to gold-standard probes present in dnam ------------------
    present = dnam.filter(pl.col("CpGmarker").is_in(gold_probes))

    # --- Impute missing gold-standard probes using their gold mean ------------
    present_probes = set(present["CpGmarker"].to_list())
    missing_probes = [p for p in gold_probes if p not in present_probes]

    if missing_probes:
        imputed_rows = pl.DataFrame(
            {
                "CpGmarker": missing_probes,
                **{
                    col: [gold_means[p] for p in missing_probes]
                    for col in sample_cols
                },
            }
        )
        present = pl.concat([present, imputed_rows], how="diagonal_relaxed")

    # Re-order to the canonical gold-probe order so ranking is consistent
    gold_order = pl.DataFrame({"CpGmarker": gold_probes})
    present = gold_order.join(present, on="CpGmarker", how="left")

    # --- Quantile-normalize each sample column in NumPy -----------------------
    # Convert to numpy for in-place column-wise operations
    matrix = present.select(sample_cols).to_numpy(allow_copy=True).astype(np.float64)

    for col_idx in range(matrix.shape[1]):
        _quantile_normalize_column(matrix[:, col_idx], sorted_target)

    # Rebuild as Polars DataFrame
    normalized = pl.DataFrame(
        {"CpGmarker": present["CpGmarker"].to_list()}
        | {sample_cols[i]: matrix[:, i] for i in range(len(sample_cols))}
    )
    return normalized


PREPROCESS_REGISTRY: dict[str, Callable[[pl.DataFrame], pl.DataFrame]] = {
    "dunedin_pace": _dunedin_pace_preprocess,
}
