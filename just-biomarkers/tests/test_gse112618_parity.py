"""Parity tests using GSE112618 (FACS validation, IlluminaEPIC, 6 blood samples).

All 5 full-match clocks must agree with biolearn to < 1e-3.
Partial-match clocks are validated for finite output and regression stability.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
import pytest
from pathlib import Path

GSE112618_CSV = Path("data/input/examples/GSE112618_methylation.csv")

FULL_MATCH_CLOCKS = ["PhenoAge", "Horvathv2", "DunedinPoAm38", "YingDamAge", "YingAdaptAge"]
PARTIAL_MATCH_CLOCKS = ["Horvathv1", "Hannum", "Lin", "Zhang_10", "YingCausAge"]

# Pre-computed from a verified run at exact parity with biolearn
EXPECTED_SCORES: dict[str, dict[str, float]] = {
    "PhenoAge":     {"GSM3074480": 34.7318, "GSM3074481": 17.2753, "GSM3074482": 14.7060},
    "Horvathv2":    {"GSM3074480": 39.8320, "GSM3074481": 24.3132, "GSM3074482": 20.2021},
    "DunedinPoAm38":{"GSM3074480":  0.9737, "GSM3074481":  0.9704, "GSM3074482":  1.0103},
    "YingDamAge":   {"GSM3074480":103.2990, "GSM3074481": 74.2227, "GSM3074482": 69.1253},
    "YingAdaptAge": {"GSM3074480":-25.3910, "GSM3074481":-21.4699, "GSM3074482":-15.8156},
    "Horvathv1":    {"GSM3074480": 53.6994, "GSM3074481": 34.4596, "GSM3074482": 33.4270},
    "Hannum":       {"GSM3074480": 38.0923, "GSM3074481": 24.0050, "GSM3074482": 19.7300},
    "Lin":          {"GSM3074480": 43.6067, "GSM3074481": 19.3869, "GSM3074482": 18.9857},
    "Zhang_10":     {"GSM3074480": -1.2439, "GSM3074481": -1.4839, "GSM3074482": -1.0306},
    "YingCausAge":  {"GSM3074480":  6.4772, "GSM3074481":-10.8438, "GSM3074482":-15.2830},
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
        return ModelGallery()
    except ImportError:
        pytest.skip("biolearn not available")


def _our_scores(dnam: pl.DataFrame, clock_name: str) -> dict[str, float]:
    from just_biomarkers.scoring import score_clock
    result = score_clock(dnam, clock_name)
    return {row["sample_id"]: row["score"] for row in result.iter_rows(named=True)}


def _biolearn_scores(gallery, geo, clock_name: str) -> dict[str, float]:
    model = gallery.get(clock_name, imputation_method="none")
    return model.predict(geo)["Predicted"].to_dict()


@pytest.mark.parametrize("clock_name", FULL_MATCH_CLOCKS)
def test_gse112618_full_match_vs_biolearn(
    clock_name: str,
    our_dnam: pl.DataFrame,
    biolearn_geo,
    biolearn_gallery,
) -> None:
    bl = _biolearn_scores(biolearn_gallery, biolearn_geo, clock_name)
    ours = _our_scores(our_dnam, clock_name)
    mismatches = [
        f"{s}: biolearn={bl[s]:.6f} ours={ours[s]:.6f} diff={abs(ours[s]-bl[s]):.2e}"
        for s in bl
        if abs(ours[s] - bl[s]) >= 1e-3
    ]
    assert not mismatches, f"{clock_name} parity failures:\n" + "\n".join(mismatches)


@pytest.mark.parametrize("clock_name", FULL_MATCH_CLOCKS + PARTIAL_MATCH_CLOCKS)
def test_gse112618_regression(clock_name: str, our_dnam: pl.DataFrame) -> None:
    ours = _our_scores(our_dnam, clock_name)
    for sample, expected in EXPECTED_SCORES[clock_name].items():
        diff = abs(ours[sample] - expected)
        assert diff < 1e-3, (
            f"{clock_name}/{sample}: expected={expected:.4f} got={ours[sample]:.4f} diff={diff:.2e}"
        )


@pytest.mark.parametrize("clock_name", PARTIAL_MATCH_CLOCKS)
def test_gse112618_partial_finite(clock_name: str, our_dnam: pl.DataFrame) -> None:
    ours = _our_scores(our_dnam, clock_name)
    assert len(ours) == 6
    for sample, score in ours.items():
        assert np.isfinite(score), f"{clock_name}/{sample}: non-finite score {score}"
