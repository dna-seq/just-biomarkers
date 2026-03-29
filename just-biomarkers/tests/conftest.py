"""Shared fixtures for just-biomarkers tests."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest


def patch_biolearn_readonly_bug() -> None:
    """Patch biolearn's quantile_normalize_using_target to handle read-only arrays.

    Installed versions of biolearn try to modify numpy arrays in-place on
    ``data.T``, which fails with newer numpy/pandas that return read-only views.
    This affects DunedinPACE parity tests.
    """
    import numpy as np
    import biolearn.dunedin_pace as _dp  # type: ignore
    orig = _dp.quantile_normalize_using_target
    if getattr(orig, "_patched_writable", False):
        return
    def _writable(data, target_values):  # noqa: ANN001
        return orig(np.array(data, dtype=np.float64), target_values)
    _writable._patched_writable = True
    _dp.quantile_normalize_using_target = _writable


def _find_biolearn_testset() -> str | None:
    """Find biolearn's testset directory -- check common locations."""
    candidates = []

    import biolearn
    pkg_dir = Path(biolearn.__file__).parent
    candidates.append(pkg_dir / "test" / "data" / "testset")

    src_candidates = [
        Path.home() / "sources" / "biolearn" / "biolearn" / "test" / "data" / "testset",
        Path("/home/antonkulaga/sources/biolearn/biolearn/test/data/testset"),
    ]
    candidates.extend(src_candidates)

    for candidate in candidates:
        methylation_file = candidate / "testset_methylation_part0.csv"
        if methylation_file.exists():
            return str(candidate)
    return None


@pytest.fixture(scope="session")
def biolearn_testset_dir() -> str:
    """Return path to biolearn's testset directory."""
    testset = _find_biolearn_testset()
    if testset is None:
        pytest.skip(
            "biolearn testset not available -- run biolearn/test/generate.py first"
        )
    return testset
