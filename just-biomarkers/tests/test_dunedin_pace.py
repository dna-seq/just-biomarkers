"""Parity tests for DunedinPACE on real EPIC data.

DunedinPACE uses quantile normalization against 20,000 gold-standard probe means
before computing the weighted sum.  These tests verify our pure-NumPy
implementation matches biolearn's scipy-based implementation to floating-point
precision (< 1e-3).

Note: biolearn ≥ some versions has a bug where `quantile_normalize_using_target`
modifies a read-only NumPy view in-place.  The fixture patches it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
import pytest
from pathlib import Path

GSE112618_CSV = Path("data/input/examples/GSE112618_methylation.csv")
GSE110554_CSV = Path("data/input/examples/GSE110554_methylation.csv")

# Pre-computed expected scores (verified at parity with biolearn)
EXPECTED_GSE112618: dict[str, float] = {
    "GSM3074480": 0.860654,
    "GSM3074481": 0.804734,
    "GSM3074482": 0.908200,
    "GSM3074483": 0.802373,
    "GSM3074484": 1.176896,
    "GSM3074485": 0.824452,
}


@pytest.fixture(scope="module", autouse=False)
def patch_biolearn_readonly():
    """Patch biolearn's quantile_normalize_using_target to work with read-only arrays."""
    try:
        import biolearn.dunedin_pace as dp
        _orig = dp.quantile_normalize_using_target

        def _patched(data, target_values):
            return _orig(np.array(data, copy=True), target_values)

        dp.quantile_normalize_using_target = _patched
        yield
        dp.quantile_normalize_using_target = _orig
    except ImportError:
        yield


@pytest.fixture(scope="module")
def gse112618_csv() -> Path:
    if not GSE112618_CSV.exists():
        pytest.skip(
            f"GSE112618 not found at {GSE112618_CSV}. "
            "Run: uv run biomarkers download-example GSE112618"
        )
    return GSE112618_CSV


@pytest.fixture(scope="module")
def gse110554_csv() -> Path:
    if not GSE110554_CSV.exists():
        pytest.skip(
            f"GSE110554 not found at {GSE110554_CSV}. "
            "Run: uv run biomarkers download-example GSE110554"
        )
    return GSE110554_CSV


@pytest.fixture(scope="module")
def biolearn_gallery(patch_biolearn_readonly):
    try:
        from biolearn.model_gallery import ModelGallery
        return ModelGallery()
    except ImportError:
        pytest.skip("biolearn not available")


def _our_pace(csv_path: Path) -> dict[str, float]:
    from just_biomarkers.io import read_methylation_csv
    from just_biomarkers.scoring import score_clock
    dnam = read_methylation_csv(csv_path)
    result = score_clock(dnam, "DunedinPACE")
    return {row["sample_id"]: row["score"] for row in result.iter_rows(named=True)}


def _biolearn_pace(csv_path: Path, gallery) -> dict[str, float]:
    from biolearn.data_library import GeoData
    df = pd.read_csv(csv_path, index_col=0)
    geo = GeoData(metadata=pd.DataFrame(index=df.columns), dnam=df)
    model = gallery.get("DunedinPACE", imputation_method="none")
    return model.predict(geo)["Predicted"].to_dict()


def test_dunedin_pace_in_registry() -> None:
    """DunedinPACE must be present in the clock registry."""
    from just_biomarkers.registry import CLOCK_DEFINITIONS
    assert "DunedinPACE" in CLOCK_DEFINITIONS
    clock = CLOCK_DEFINITIONS["DunedinPACE"]
    assert clock.preprocess_name == "dunedin_pace"


def test_dunedin_pace_regression_gse112618(gse112618_csv: Path) -> None:
    """Scores must not drift from pre-computed expected values."""
    ours = _our_pace(gse112618_csv)
    for sample, expected in EXPECTED_GSE112618.items():
        diff = abs(ours[sample] - expected)
        assert diff < 1e-3, (
            f"DunedinPACE/{sample}: expected={expected:.6f} got={ours[sample]:.6f} diff={diff:.2e}"
        )


def test_dunedin_pace_parity_gse112618(
    gse112618_csv: Path,
    biolearn_gallery,
    patch_biolearn_readonly,
) -> None:
    """All 6 samples must match biolearn to < 1e-3."""
    ours = _our_pace(gse112618_csv)
    bl = _biolearn_pace(gse112618_csv, biolearn_gallery)
    mismatches = [
        f"{s}: ours={ours[s]:.6f} biolearn={bl[s]:.6f} diff={abs(ours[s]-bl[s]):.2e}"
        for s in bl
        if abs(ours[s] - bl[s]) >= 1e-3
    ]
    assert not mismatches, "DunedinPACE parity failures:\n" + "\n".join(mismatches)


def test_dunedin_pace_parity_gse110554(
    gse110554_csv: Path,
    biolearn_gallery,
    patch_biolearn_readonly,
) -> None:
    """All 49 samples on FlowSorted.Blood.EPIC must match biolearn to < 1e-3."""
    ours = _our_pace(gse110554_csv)
    bl = _biolearn_pace(gse110554_csv, biolearn_gallery)
    mismatches = [
        f"{s}: ours={ours[s]:.6f} biolearn={bl[s]:.6f} diff={abs(ours[s]-bl[s]):.2e}"
        for s in bl
        if abs(ours[s] - bl[s]) >= 1e-3
    ]
    assert not mismatches, (
        f"DunedinPACE parity failures ({len(mismatches)} samples):\n"
        + "\n".join(mismatches)
    )


def test_dunedin_pace_finite_gse112618(gse112618_csv: Path) -> None:
    """All scores must be finite."""
    ours = _our_pace(gse112618_csv)
    bad = [s for s, v in ours.items() if not np.isfinite(v)]
    assert not bad, f"Non-finite DunedinPACE scores: {bad}"


def test_dunedin_pace_range_gse112618(gse112618_csv: Path) -> None:
    """DunedinPACE (aging rate) should be in a plausible range (0.5 – 2.0 yr/yr)."""
    ours = _our_pace(gse112618_csv)
    out_of_range = [(s, v) for s, v in ours.items() if not (0.5 <= v <= 2.0)]
    assert not out_of_range, (
        f"DunedinPACE scores outside [0.5, 2.0]: {out_of_range}"
    )
