"""Parity tests using GSE110554 (FlowSorted.Blood.EPIC, 49 samples).

Compares just-biomarkers scores against biolearn for the same real EPIC dataset.

Two test categories:
  1. FULL-MATCH clocks — clocks where GSE110554 covers 100% of required CpGs.
     Scores must match biolearn to floating-point precision (< 1e-3).
  2. PARTIAL-MATCH clocks — clocks where GSE110554 is missing some CpGs.
     Biolearn refuses to run these without imputation.  We verify our scores are:
       a) numerically stable (finite, non-NaN)
       b) consistent across runs (deterministic)
       c) within a plausible biological range for age-predicting clocks

The GSE110554 CSV is expected at data/input/examples/GSE110554_methylation.csv.
Tests are skipped if the file is absent (download with ``uv run biomarkers download-example GSE110554``).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pytest

BIOLEARN_PATH = Path("/home/antonkulaga/sources/biolearn")
GSE110554_CSV = Path("data/input/examples/GSE110554_methylation.csv")

# 5 representative samples used in per-sample assertions
SAMPLES_TO_CHECK = [
    "GSM2998021",
    "GSM2998022",
    "GSM2998023",
    "GSM2998024",
    "GSM2998025",
]

# Clocks where GSE110554 covers 100% required CpGs — biolearn runs without imputation
FULL_MATCH_CLOCKS = [
    "PhenoAge",
    "Horvathv2",
    "DunedinPoAm38",
    "YingDamAge",
    "YingAdaptAge",
]

# Clocks with missing CpGs on EPIC — biolearn errors, we compute partial
PARTIAL_MATCH_CLOCKS = [
    "Horvathv1",   # 334/353 CpGs on EPIC
    "Hannum",      # 65/71 CpGs on EPIC
    "Lin",         # 97/99 CpGs on EPIC
    "Zhang_10",    # 8/10 CpGs on EPIC
    "YingCausAge", # 420/585 CpGs on EPIC
]

# Expected scores pre-computed from a verified run (biolearn or our engine at parity).
# Used as regression anchors if biolearn is unavailable.
EXPECTED_SCORES: dict[str, dict[str, float]] = {
    "PhenoAge":     {"GSM2998021": 44.7695, "GSM2998022": 38.6062, "GSM2998023": 24.7768},
    "Horvathv2":    {"GSM2998021": 39.5018, "GSM2998022": 49.2316, "GSM2998023": 18.8622},
    "DunedinPoAm38":{"GSM2998021":  1.0481, "GSM2998022":  0.8835, "GSM2998023":  0.9734},
    "YingDamAge":   {"GSM2998021": 77.5576, "GSM2998022": 65.4661, "GSM2998023": 46.4510},
    "YingAdaptAge": {"GSM2998021": 11.9261, "GSM2998022": -0.3330, "GSM2998023": 14.1345},
    "Horvathv1":    {"GSM2998021": 44.6807, "GSM2998022": 54.4222, "GSM2998023": 28.0145},
    "Hannum":       {"GSM2998021": 36.5533, "GSM2998022": 45.4201, "GSM2998023": 24.8031},
    "Lin":          {"GSM2998021": 34.8933, "GSM2998022": 36.7937, "GSM2998023": 17.4296},
    "Zhang_10":     {"GSM2998021":  0.3020, "GSM2998022": -0.8857, "GSM2998023":  0.2897},
    "YingCausAge":  {"GSM2998021":  2.1972, "GSM2998022":  3.1647, "GSM2998023": -15.3587},
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def gse110554_csv() -> Path:
    if not GSE110554_CSV.exists():
        pytest.skip(
            f"GSE110554 not found at {GSE110554_CSV}. "
            "Run: uv run biomarkers download-example GSE110554"
        )
    return GSE110554_CSV


@pytest.fixture(scope="module")
def our_dnam(gse110554_csv: Path) -> pl.DataFrame:
    from just_biomarkers.io import read_methylation_csv
    return read_methylation_csv(gse110554_csv)


@pytest.fixture(scope="module")
def biolearn_geo(gse110554_csv: Path):
    """GeoData built from GSE110554 CSV for biolearn scoring."""
    if BIOLEARN_PATH not in sys.path:
        sys.path.insert(0, str(BIOLEARN_PATH))
    try:
        from biolearn.data_library import GeoData
    except ImportError:
        pytest.skip("biolearn not available")

    df = pd.read_csv(gse110554_csv, index_col=0)
    return GeoData(metadata=pd.DataFrame(index=df.columns), dnam=df)


@pytest.fixture(scope="module")
def biolearn_gallery(biolearn_geo):
    if BIOLEARN_PATH not in sys.path:
        sys.path.insert(0, str(BIOLEARN_PATH))
    try:
        from biolearn.model_gallery import ModelGallery
        return ModelGallery()
    except ImportError:
        pytest.skip("biolearn not available")


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _our_scores(our_dnam: pl.DataFrame, clock_name: str) -> dict[str, float]:
    from just_biomarkers.scoring import score_clock
    result = score_clock(our_dnam, clock_name)
    return {row["sample_id"]: row["score"] for row in result.iter_rows(named=True)}


def _biolearn_scores(gallery, geo, clock_name: str) -> dict[str, float]:
    model = gallery.get(clock_name, imputation_method="none")
    result = model.predict(geo)
    return result["Predicted"].to_dict()


# ---------------------------------------------------------------------------
# Full-match parity: must match biolearn to floating-point precision
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("clock_name", FULL_MATCH_CLOCKS)
def test_full_match_parity_vs_biolearn(
    clock_name: str,
    our_dnam: pl.DataFrame,
    biolearn_geo,
    biolearn_gallery,
) -> None:
    """Scores for full-match clocks must be identical to biolearn (tol 1e-3)."""
    bl_scores = _biolearn_scores(biolearn_gallery, biolearn_geo, clock_name)
    our_scores = _our_scores(our_dnam, clock_name)

    mismatches: list[str] = []
    for sample in SAMPLES_TO_CHECK:
        expected = bl_scores[sample]
        actual = our_scores[sample]
        diff = abs(actual - expected)
        if diff >= 1e-3:
            mismatches.append(
                f"{sample}: biolearn={expected:.6f} ours={actual:.6f} diff={diff:.2e}"
            )

    assert not mismatches, (
        f"{clock_name} parity failures:\n" + "\n".join(mismatches)
    )


@pytest.mark.parametrize("clock_name", FULL_MATCH_CLOCKS)
def test_full_match_all_samples_vs_biolearn(
    clock_name: str,
    our_dnam: pl.DataFrame,
    biolearn_geo,
    biolearn_gallery,
) -> None:
    """All 49 samples must match biolearn for full-match clocks."""
    bl_scores = _biolearn_scores(biolearn_gallery, biolearn_geo, clock_name)
    our_scores = _our_scores(our_dnam, clock_name)

    mismatches: list[str] = []
    for sample, expected in bl_scores.items():
        actual = our_scores.get(sample)
        assert actual is not None, f"Sample {sample} missing from our output"
        diff = abs(actual - expected)
        if diff >= 1e-3:
            mismatches.append(f"{sample}: diff={diff:.2e}")

    assert not mismatches, (
        f"{clock_name}: {len(mismatches)} samples out of tolerance:\n"
        + "\n".join(mismatches)
    )


# ---------------------------------------------------------------------------
# Full-match regression: check against stored expected values (no biolearn)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("clock_name", FULL_MATCH_CLOCKS)
def test_full_match_regression(
    clock_name: str,
    our_dnam: pl.DataFrame,
) -> None:
    """Regression guard: scores must not drift from pre-computed expected values."""
    our_scores = _our_scores(our_dnam, clock_name)
    for sample, expected in EXPECTED_SCORES[clock_name].items():
        actual = our_scores[sample]
        diff = abs(actual - expected)
        assert diff < 1e-3, (
            f"{clock_name}/{sample}: expected={expected:.4f} got={actual:.4f} diff={diff:.2e}"
        )


# ---------------------------------------------------------------------------
# Partial-match: stable, finite, regression
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("clock_name", PARTIAL_MATCH_CLOCKS)
def test_partial_match_finite(
    clock_name: str,
    our_dnam: pl.DataFrame,
) -> None:
    """Partial-match clocks must produce finite scores for all samples."""
    our_scores = _our_scores(our_dnam, clock_name)
    assert len(our_scores) == 49, f"{clock_name}: expected 49 samples, got {len(our_scores)}"
    for sample, score in our_scores.items():
        assert np.isfinite(score), f"{clock_name}/{sample}: non-finite score {score}"


@pytest.mark.parametrize("clock_name", PARTIAL_MATCH_CLOCKS)
def test_partial_match_regression(
    clock_name: str,
    our_dnam: pl.DataFrame,
) -> None:
    """Regression guard: partial-match scores must not drift from expected values."""
    our_scores = _our_scores(our_dnam, clock_name)
    for sample, expected in EXPECTED_SCORES[clock_name].items():
        actual = our_scores[sample]
        diff = abs(actual - expected)
        assert diff < 1e-3, (
            f"{clock_name}/{sample}: expected={expected:.4f} got={actual:.4f} diff={diff:.2e}"
        )


@pytest.mark.parametrize("clock_name", ["Horvathv1", "Hannum", "Lin"])
def test_partial_match_age_plausible(
    clock_name: str,
    our_dnam: pl.DataFrame,
) -> None:
    """Age-predicting clocks must yield plausible values (0–120 years) for most samples.

    EPIC blood cell-type reference samples are from adults, so most scores should
    be in a realistic human adult age range.
    """
    our_scores = _our_scores(our_dnam, clock_name)
    scores = list(our_scores.values())
    in_range = sum(0 <= s <= 120 for s in scores)
    assert in_range >= len(scores) * 0.8, (
        f"{clock_name}: only {in_range}/{len(scores)} scores in [0, 120]. "
        f"Min={min(scores):.2f} Max={max(scores):.2f}"
    )


@pytest.mark.parametrize("clock_name", PARTIAL_MATCH_CLOCKS)
def test_partial_match_match_rate(
    clock_name: str,
    our_dnam: pl.DataFrame,
) -> None:
    """Match rate for partial clocks must be above 50% (enough CpGs to compute)."""
    from just_biomarkers.scoring import score_clock
    result = score_clock(our_dnam, clock_name)
    for row in result.iter_rows(named=True):
        rate = row["cpgs_matched"] / max(row["cpgs_required"], 1)
        assert rate >= 0.5, (
            f"{clock_name}/{row['sample_id']}: match rate {rate:.1%} is too low "
            f"({row['cpgs_matched']}/{row['cpgs_required']})"
        )
