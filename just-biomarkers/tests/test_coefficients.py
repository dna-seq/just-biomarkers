"""Offline tests for the coefficient downloader (served from a local file:// URL)."""
from __future__ import annotations

import re
from pathlib import Path

import polars as pl
import pytest

from just_biomarkers import coefficients
from just_biomarkers.coefficients import (
    COEFFICIENT_FILES,
    CoefficientChecksumError,
    coefficient_path,
    download_coefficient,
    git_blob_sha1,
)
from just_biomarkers.registry import CLOCK_DEFINITIONS

FAKE_CSV = b"CpGmarker,CoefficientTraining\nintercept,1.5\ncg00000001,2.0\n"


@pytest.fixture
def served(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> str:
    """Serve one fake coefficient file from a local directory; cache under tmp."""
    source = tmp_path / "source"
    source.mkdir()
    (source / "Fake.csv").write_bytes(FAKE_CSV)
    monkeypatch.setenv("BIOMARKERS_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setitem(COEFFICIENT_FILES, "Fake.csv", git_blob_sha1(FAKE_CSV))
    monkeypatch.setattr(coefficients, "BIOLEARN_DATA_URL", source.as_uri())
    monkeypatch.setattr(coefficients, "package_data_dir", lambda: tmp_path / "no-bundle")
    return source.as_uri()


def test_git_blob_sha1_matches_git() -> None:
    # `printf 'hello\n' | git hash-object --stdin`
    assert git_blob_sha1(b"hello\n") == "ce013625030ba8dba906f756967f9e9ca394464a"


def test_manifest_covers_every_registered_clock() -> None:
    needed = {c.coefficient_file for c in CLOCK_DEFINITIONS.values()}
    needed.add("DunedinPACE_Gold_Means.csv")
    assert needed <= set(COEFFICIENT_FILES)
    assert all(re.fullmatch(r"[0-9a-f]{40}", sha) for sha in COEFFICIENT_FILES.values())


def test_download_verifies_and_caches(served: str, tmp_path: Path) -> None:
    path = download_coefficient("Fake.csv")
    assert path.read_bytes() == FAKE_CSV
    assert path.is_relative_to(tmp_path / "cache")
    # Second call must reuse the cache: remove the source to prove no download.
    (tmp_path / "source" / "Fake.csv").unlink()
    assert download_coefficient("Fake.csv") == path


def test_checksum_mismatch_raises_and_saves_nothing(
    served: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(COEFFICIENT_FILES, "Fake.csv", "0" * 40)
    with pytest.raises(CoefficientChecksumError):
        download_coefficient("Fake.csv")
    assert not list((tmp_path / "cache").rglob("Fake.csv*"))


def test_corrupt_cache_is_replaced(served: str) -> None:
    path = download_coefficient("Fake.csv")
    path.write_bytes(b"truncated")
    assert download_coefficient("Fake.csv").read_bytes() == FAKE_CSV


def test_unknown_file_is_refused(served: str) -> None:
    with pytest.raises(KeyError, match="Unknown coefficient file"):
        download_coefficient("NotAClock.csv")


def test_bundled_file_wins(served: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    (bundle / "Fake.csv").write_bytes(b"CpGmarker,Weight\nintercept,9\n")
    monkeypatch.setattr(coefficients, "package_data_dir", lambda: bundle)
    assert coefficient_path("Fake.csv") == bundle / "Fake.csv"


def test_load_coefficients_downloads_on_demand(served: str) -> None:
    from just_biomarkers.data_utils import load_coefficients

    df = load_coefficients("Fake.csv")
    assert df.columns[0] == "CpGmarker"
    assert df.filter(pl.col("CpGmarker") == "cg00000001")["CoefficientTraining"][0] == 2.0
