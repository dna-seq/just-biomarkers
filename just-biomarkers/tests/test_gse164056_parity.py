"""Parity tests using GSE164056 (Social Anxiety Disorder, IlluminaEPIC, 143 blood samples).

GSE164056 has ~690k CpGs (vs 866k for full EPIC), so ALL clocks are partial-match here —
biolearn cannot score any without imputation. We validate:
  - finite, non-NaN output for all 143 samples
  - regression stability (scores don't drift from pre-computed anchors)
  - age-predicting clocks stay within plausible biological range
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from pathlib import Path

GSE164056_CSV = Path("data/input/examples/GSE164056_methylation.csv")

ALL_CLOCKS = [
    "PhenoAge", "Horvathv2", "DunedinPoAm38", "YingDamAge", "YingAdaptAge",
    "Horvathv1", "Hannum", "Lin", "Zhang_10", "YingCausAge",
]

AGE_CLOCKS = ["PhenoAge", "Horvathv2", "Horvathv1", "Hannum", "Lin"]

# Regression anchors from first 3 samples (GSM4995676, GSM4995677, GSM4995678)
EXPECTED_SCORES: dict[str, dict[str, float]] = {
    "PhenoAge":     {"GSM4995676": 42.8678, "GSM4995677": 50.4090, "GSM4995678": 34.2430},
    "Horvathv2":    {"GSM4995676": 31.0515, "GSM4995677": 38.2961, "GSM4995678": 20.3871},
    "DunedinPoAm38":{"GSM4995676":  0.6823, "GSM4995677":  0.7699, "GSM4995678":  0.6922},
    "YingDamAge":   {"GSM4995676": 46.1073, "GSM4995677": 57.5274, "GSM4995678": 30.7908},
    "YingAdaptAge": {"GSM4995676": 27.1000, "GSM4995677": 29.8899, "GSM4995678": 26.8303},
    "Horvathv1":    {"GSM4995676": 40.5345, "GSM4995677": 48.8599, "GSM4995678": 32.6033},
    "Hannum":       {"GSM4995676": 20.6702, "GSM4995677": 28.8392, "GSM4995678": 14.0798},
    "Lin":          {"GSM4995676": 26.1162, "GSM4995677": 29.0862, "GSM4995678": 13.3098},
    "Zhang_10":     {"GSM4995676": -1.3992, "GSM4995677": -1.4350, "GSM4995678": -1.5290},
    "YingCausAge":  {"GSM4995676": -5.7213, "GSM4995677": -3.9984, "GSM4995678":-20.6908},
}


@pytest.fixture(scope="module")
def gse164056_csv() -> Path:
    if not GSE164056_CSV.exists():
        pytest.skip(
            f"GSE164056 not found at {GSE164056_CSV}. "
            "Run: uv run biomarkers download-example GSE164056"
        )
    return GSE164056_CSV


@pytest.fixture(scope="module")
def our_dnam(gse164056_csv: Path) -> pl.DataFrame:
    from just_biomarkers.io import read_methylation_csv
    return read_methylation_csv(gse164056_csv)


def _our_scores(dnam: pl.DataFrame, clock_name: str) -> dict[str, float]:
    from just_biomarkers.scoring import score_clock
    result = score_clock(dnam, clock_name)
    return {row["sample_id"]: row["score"] for row in result.iter_rows(named=True)}


@pytest.mark.parametrize("clock_name", ALL_CLOCKS)
def test_gse164056_finite(clock_name: str, our_dnam: pl.DataFrame) -> None:
    """All 143 samples must yield finite scores."""
    ours = _our_scores(our_dnam, clock_name)
    assert len(ours) == 143, f"{clock_name}: expected 143 samples, got {len(ours)}"
    bad = [s for s, v in ours.items() if not np.isfinite(v)]
    assert not bad, f"{clock_name}: non-finite scores for {bad[:5]}"


@pytest.mark.parametrize("clock_name", ALL_CLOCKS)
def test_gse164056_regression(clock_name: str, our_dnam: pl.DataFrame) -> None:
    """Scores must not drift from pre-computed anchor values."""
    ours = _our_scores(our_dnam, clock_name)
    for sample, expected in EXPECTED_SCORES[clock_name].items():
        diff = abs(ours[sample] - expected)
        assert diff < 1e-3, (
            f"{clock_name}/{sample}: expected={expected:.4f} got={ours[sample]:.4f} diff={diff:.2e}"
        )


@pytest.mark.parametrize("clock_name", AGE_CLOCKS)
def test_gse164056_age_plausible(clock_name: str, our_dnam: pl.DataFrame) -> None:
    """Age-predicting clocks must yield mostly plausible values (0-120 years)."""
    ours = _our_scores(our_dnam, clock_name)
    scores = list(ours.values())
    in_range = sum(0 <= s <= 120 for s in scores)
    assert in_range >= len(scores) * 0.8, (
        f"{clock_name}: only {in_range}/{len(scores)} scores in [0, 120]. "
        f"Min={min(scores):.2f} Max={max(scores):.2f}"
    )


@pytest.mark.parametrize("clock_name", ALL_CLOCKS)
def test_gse164056_match_rate(clock_name: str, our_dnam: pl.DataFrame) -> None:
    """Match rates must all be above 50% (enough CpGs to produce meaningful scores)."""
    from just_biomarkers.scoring import score_clock
    result = score_clock(our_dnam, clock_name)
    low = [
        (row["sample_id"], row["cpgs_matched"] / max(row["cpgs_required"], 1))
        for row in result.iter_rows(named=True)
        if row["cpgs_matched"] / max(row["cpgs_required"], 1) < 0.5
    ]
    assert not low, (
        f"{clock_name}: {len(low)} samples below 50% match rate"
    )
