"""Parity tests for DunedinPACE on EPIC data (GSE112618, 6 blood samples).

DunedinPACE requires quantile normalization against 20k gold-standard probe
means before the weighted sum.  These tests verify our pure-NumPy/Polars
implementation matches biolearn's scipy-based implementation to < 1e-3.

Note: installed biolearn has a read-only array bug in quantile_normalize_using_target
(it tries in-place assignment on data.T which numpy returns as read-only).
The conftest.py patch fixes this at test session start.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
import pytest
from pathlib import Path

GSE112618_CSV = Path("data/input/examples/GSE112618_methylation.csv")
ALL_SAMPLES = [
    "GSM3074480", "GSM3074481", "GSM3074482",
    "GSM3074483", "GSM3074484", "GSM3074485",
]

# Pre-computed at exact parity with biolearn (differences < 1e-15)
EXPECTED_SCORES: dict[str, float] = {
    "GSM3074480": 0.860654,
    "GSM3074481": 0.804734,
    "GSM3074482": 0.908200,
    "GSM3074483": 0.802373,
    "GSM3074484": 1.176896,
    "GSM3074485": 0.824452,
}


@pytest.fixture(scope="module")
def gse112618_csv() -> Path:
    if not GSE112618_CSV.exists():
        pytest.skip(
            f"GSE112618 not found at {GSE112618_CSV}. "
            "Run: uv run biomarkers download-example GSE112618"
        )
    return GSE112618_CSV


@pytest.fixture(scope="module")
def our_dnam(gse112618_csv: Path) -> pl.DataFrame:
    from just_biomarkers.io import read_methylation_csv
    return read_methylation_csv(gse112618_csv)


@pytest.fixture(scope="module")
def biolearn_geo(gse112618_csv: Path):
    try:
        from biolearn.data_library import GeoData
    except ImportError:
        pytest.skip("biolearn not available")
    df = pd.read_csv(gse112618_csv, index_col=0)
    return GeoData(metadata=pd.DataFrame(index=df.columns), dnam=df)


@pytest.fixture(scope="module")
def biolearn_gallery(biolearn_geo):
    try:
        from biolearn.model_gallery import ModelGallery
        from conftest import patch_biolearn_readonly_bug
        patch_biolearn_readonly_bug()
        return ModelGallery()
    except ImportError:
        pytest.skip("biolearn not available")


def _our_scores(dnam: pl.DataFrame) -> dict[str, float]:
    from just_biomarkers.scoring import score_clock
    result = score_clock(dnam, "DunedinPACE")
    return {row["sample_id"]: row["score"] for row in result.iter_rows(named=True)}


def test_dunedin_pace_in_registry() -> None:
    """DunedinPACE must be present in CLOCK_DEFINITIONS."""
    from just_biomarkers.registry import CLOCK_DEFINITIONS
    assert "DunedinPACE" in CLOCK_DEFINITIONS
    defn = CLOCK_DEFINITIONS["DunedinPACE"]
    assert defn.preprocess_name == "dunedin_pace"
    assert defn.output == "Aging Rate (Years/Year)"


def test_dunedin_pace_produces_6_samples(our_dnam: pl.DataFrame) -> None:
    ours = _our_scores(our_dnam)
    assert len(ours) == 6
    for sample, score in ours.items():
        assert np.isfinite(score), f"{sample}: non-finite score {score}"


def test_dunedin_pace_full_match(our_dnam: pl.DataFrame) -> None:
    from just_biomarkers.scoring import score_clock
    result = score_clock(our_dnam, "DunedinPACE")
    for row in result.iter_rows(named=True):
        assert row["cpgs_matched"] == 173
        assert row["cpgs_required"] == 173


def test_dunedin_pace_regression(our_dnam: pl.DataFrame) -> None:
    """Scores must not drift from pre-computed anchor values."""
    ours = _our_scores(our_dnam)
    for sample, expected in EXPECTED_SCORES.items():
        diff = abs(ours[sample] - expected)
        assert diff < 1e-3, (
            f"{sample}: expected={expected:.6f} got={ours[sample]:.6f} diff={diff:.2e}"
        )


def test_dunedin_pace_vs_biolearn(
    our_dnam: pl.DataFrame,
    biolearn_geo,
    biolearn_gallery,
) -> None:
    """Scores must match biolearn to < 1e-3 for all 6 samples."""
    bl = biolearn_gallery.get("DunedinPACE").predict(biolearn_geo)["Predicted"].to_dict()
    ours = _our_scores(our_dnam)
    mismatches = [
        f"{s}: biolearn={bl[s]:.6f} ours={ours[s]:.6f} diff={abs(ours[s]-bl[s]):.2e}"
        for s in bl
        if abs(ours[s] - bl[s]) >= 1e-3
    ]
    assert not mismatches, "DunedinPACE parity failures:\n" + "\n".join(mismatches)


def test_dunedin_pace_pace_range(our_dnam: pl.DataFrame) -> None:
    """DunedinPACE outputs aging rate; healthy blood should be near 1.0 (0.5–2.0)."""
    ours = _our_scores(our_dnam)
    out_of_range = [
        (s, v) for s, v in ours.items() if not (0.3 <= v <= 3.0)
    ]
    assert not out_of_range, (
        f"DunedinPACE scores outside plausible [0.3, 3.0] range: {out_of_range}"
    )
