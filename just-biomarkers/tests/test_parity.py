"""Parity tests: compare just-biomarkers Polars engine vs biolearn pandas engine.

These tests run the same input methylation matrix through both libraries for
~10 representative linear clocks and assert that scores match within practical
tolerance (1e-4 to 1e-3 depending on preprocessing path).
"""
from __future__ import annotations

import os
from collections import defaultdict

import numpy as np
import pandas as pd
import polars as pl
import pytest

TOLERANCES: dict[str, float] = defaultdict(lambda: 1e-4)
TOLERANCES["Horvathv1"] = 1e-3
TOLERANCES["Horvathv2"] = 1e-3
TOLERANCES["PEDBE"] = 1e-3
TOLERANCES["DNAmClockCortical"] = 1e-3

LINEAR_CLOCKS_TO_TEST = [
    "Horvathv1",
    "Hannum",
    "PhenoAge",
    "Lin",
    "Horvathv2",
    "DunedinPoAm38",
    "Zhang_10",
    "YingCausAge",
    "YingDamAge",
    "YingAdaptAge",
]


@pytest.fixture(scope="module")
def biolearn_test_data(biolearn_testset_dir: str):
    """Load biolearn GeoData from the testset."""
    from biolearn.data_library import GeoData

    data = GeoData.load_csv(biolearn_testset_dir, "testset", validate=False)
    return data


@pytest.fixture(scope="module")
def polars_dnam(biolearn_test_data) -> pl.DataFrame:
    """Convert biolearn's pandas dnam matrix to our Polars format."""
    pandas_dnam: pd.DataFrame = biolearn_test_data.dnam
    pandas_dnam_reset = pandas_dnam.reset_index()
    pandas_dnam_reset.rename(columns={pandas_dnam_reset.columns[0]: "CpGmarker"}, inplace=True)
    return pl.from_pandas(pandas_dnam_reset)


@pytest.mark.parametrize("clock_name", LINEAR_CLOCKS_TO_TEST)
def test_linear_clock_parity(
    clock_name: str,
    biolearn_test_data,
    polars_dnam: pl.DataFrame,
) -> None:
    """Score one clock via both engines and compare per-sample results."""
    from biolearn.model_gallery import ModelGallery

    gallery = ModelGallery()

    biolearn_model = gallery.get(clock_name, imputation_method="none")
    required_cpgs = biolearn_model.methylation_sites()
    available_cpgs = set(biolearn_test_data.dnam.index)
    missing = set(required_cpgs) - available_cpgs
    if missing:
        pytest.skip(
            f"{clock_name}: {len(missing)}/{len(required_cpgs)} CpGs missing from testset"
        )

    biolearn_result = biolearn_model.predict(biolearn_test_data)
    biolearn_scores = biolearn_result["Predicted"].sort_index()

    from just_biomarkers.scoring import score_clock

    our_result = score_clock(polars_dnam, clock_name)

    our_scores_dict = {
        row["sample_id"]: row["score"]
        for row in our_result.iter_rows(named=True)
    }

    tol = TOLERANCES[clock_name]

    for sample_id in biolearn_scores.index:
        expected = float(biolearn_scores[sample_id])
        actual = our_scores_dict.get(sample_id)
        assert actual is not None, (
            f"Sample {sample_id} missing from just-biomarkers output"
        )
        assert abs(actual - expected) < tol, (
            f"{clock_name} / {sample_id}: "
            f"expected={expected:.6f}, got={actual:.6f}, "
            f"diff={abs(actual - expected):.2e}, tol={tol:.0e}"
        )


KNOWN_BIOLEARN_LINEAR_CLOCKS = {
    "Horvathv1", "Hannum", "PhenoAge", "Lin", "Horvathv2",
    "DunedinPoAm38", "DunedinPACE", "Zhang_10",
    "YingCausAge", "YingDamAge", "YingAdaptAge",
    "PEDBE", "DNAmClockCortical", "VidalBralo", "DNAmTL",
    "AlcoholMcCartney", "SmokingMcCartney",
    "BMI_McCartney", "EducationMcCartney", "BodyFatMcCartney",
    "HDLCholesterolMcCartney", "LDLCholesterolMcCartney",
    "TotalCholesterolMcCartney",
    "StocZ", "StocP", "StocH",
    "HRSInCHPhenoAge", "EpiTOC1",
    "Knight", "LeeControl", "LeeRefinedRobust", "LeeRobust",
    "BMI_Reed", "CVD_Westerman", "AD_Bahado-Singh",
    "DepressionBarbu", "Weidner", "Garagnani",
    "Mayne", "ProstateCancerKirby", "HepatoXu",
    "Bohlin", "Bocklandt", "DownSyndrome",
}


def test_clock_count_parity() -> None:
    """Verify our registry covers a reasonable subset of biolearn linear clocks."""
    from just_biomarkers.registry import CLOCK_DEFINITIONS

    our_clocks = set(CLOCK_DEFINITIONS.keys())
    covered = KNOWN_BIOLEARN_LINEAR_CLOCKS & our_clocks
    assert len(covered) >= 20, (
        f"Expected at least 20 biolearn linear clocks covered, got {len(covered)}. "
        f"Missing: {KNOWN_BIOLEARN_LINEAR_CLOCKS - our_clocks}"
    )
