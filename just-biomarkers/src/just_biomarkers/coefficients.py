"""Download clock coefficient files from biolearn without a biolearn dependency.

The coefficient CSVs (and DunedinPACE's gold-standard means) are not kept in
git.  They are fetched from biolearn's ``biolearn/data/`` directory at a
pinned commit, checked against the git blob SHA-1 recorded below, and cached
under ``<cache>/coefficients/<commit>/``, where ``<cache>`` is
``BIOMARKERS_CACHE_DIR`` or ``platformdirs.user_cache_dir("just-biomarkers")``.

A file found in the package's own ``data/`` directory (a wheel that bundles
them) is used first; otherwise the cached copy; otherwise it is downloaded.
A download whose checksum does not match is discarded and raises.
"""
from __future__ import annotations

import hashlib
import os
import urllib.request
from pathlib import Path
from typing import Iterable, Optional

import platformdirs

BIOLEARN_REPO = "bio-learn/biolearn"
BIOLEARN_COMMIT = "0d714f5a0c0ad5743b6315f3edab1f86112d394f"
BIOLEARN_DATA_URL = (
    f"https://raw.githubusercontent.com/{BIOLEARN_REPO}/{BIOLEARN_COMMIT}/biolearn/data"
)

#: File name -> git blob SHA-1 at :data:`BIOLEARN_COMMIT`.
COEFFICIENT_FILES: dict[str, str] = {
    "AD_Bahado-Singh.csv": "049b30610995c8b2b89c852581572c8f5aa7efae",
    "Alcohol.csv": "2f654cd350f949417c5e923cde27d7e250d1b6db",
    "BMI_McCartney.csv": "7c6c46b2aac0f2073265ffa9b6e39d2983e04164",
    "BMI_Reed.csv": "4b1003c993e1eb0bf9b011c6b1c9ea819d29a7c3",
    "Bocklandt.csv": "8d0dfeb895cb35b2b61fdb941792914077a1b4bf",
    "BodyFatMcCartney.csv": "aeb1ad0fa1621143f19e1b19540c2bf65ad69c59",
    "Bohlin.csv": "d28789529152c0eb9442ad4d7bcca56ef2dedab4",
    "CVD_Westermann.csv": "a66307c4cf6ce6263c8771ffdbade10041f44f4e",
    "DNAmClockCortical.csv": "fb63f3da68eaa47355f117ee4ebc6c86fe29fdb6",
    "DNAmTL.csv": "c8c45081b7c5bdf9f9409fc72216d6e1e258fcfe",
    "DepressionBarbu.csv": "9ae664566d5a08a20561ab69418b883316a73fda",
    "DunedinPACE.csv": "b65bd0933816d7f074b35907216c235630ce028b",
    "DunedinPACE_Gold_Means.csv": "9aabb8d382716c60595df80991d65f9ea8877db8",
    "DunedinPoAm38.csv": "81543a8bbc4e3dc22b114253a1b0bccfaddea5f8",
    "EducationMcCartney.csv": "035fdefa750fa2e57c6255e0653d10dc8250df12",
    "EpiTOC1.csv": "dbf53392a7180a2bbd1eaeefc7f00338189853bf",
    "Garagnani.csv": "c10f4a9bfbf19ebedf05a13c4a186252d3d82006",
    "HDLCholesterolMcCartney.csv": "62aeefc4f4e97efe2aa3f95dbf404ae3ec27c959",
    "HRSInCHPhenoAge.csv": "a6f790c08199f1ebf863d083d65d1a52e2376587",
    "Hannum.csv": "1a26380822523748d1092d8a9f1ccae4f4494e68",
    "HepatoXu.csv": "cd6e406c2580a8c46bd5c834cf633812cee57080",
    "Horvath1.csv": "2c889c852087374b6bd0c6cdd19545beb226e029",
    "Horvath2.csv": "fea1667ce4481c4969f50318d72e4b8048b50993",
    "Knight.csv": "b44e4129814f24b73cfc5d8bee7193d5ad98afba",
    "LDLCholesterolMcCartney.csv": "978c36325d13ff7eed29b073d4d6fee714fc0f01",
    "LeeControl.csv": "5e7e0c5da09522a1d2c5d3ac06e7f36c8c6166ae",
    "LeeRefinedRobust.csv": "e4e0dbb91b38a28aca20242992cb21fc947dc88c",
    "LeeRobust.csv": "face5ec055e89dcf29db9550a6f48a3aeb38b951",
    "Lin.csv": "f575affa18a0f54dca7bd2792f6387c9c124273b",
    "Mayne.csv": "eaceddba3ffa5cd8c03032543b12bf17c6db067f",
    "PEDBE.csv": "ef8ecef22a065fdf2308fabb3eb50cc8597b2924",
    "PhenoAge.csv": "2362abc9c736996056a95ff5cb3f02f08276ec1f",
    "ProstateCancerKirby.csv": "058dee1c88c14ea1d42b57626ae4fb78eab8a32e",
    "Smoking.csv": "8ea3fe03baa6d95bd2adc7e22b6b01dc24ea2a01",
    "StocH.csv": "e025e52cd84ade2113fa617f7ad4a54390e403a5",
    "StocP.csv": "3c3f565bbc27675dea83ddb568b9d06e8ee60e83",
    "StocZ.csv": "1e37c7940776ce1e2b6131eeb524e3b4b883eb9a",
    "TotalCholesterolMcCartney.csv": "1fa27c333cf715bfd2886ba245e95afca8f5b637",
    "VidalBralo.csv": "e61ccf426bd5f4ea32e56c9f4983c5095d3e633e",
    "Weidner.csv": "9da00acec18321a7ea5d5041941fe9f750fdd87c",
    "YingAdaptAge.csv": "87652423f0c8faa6e1f5450d9b98a2c2ab36375e",
    "YingCausAge.csv": "27fcc77f54110d51c129edea578710160355f6b8",
    "YingDamAge.csv": "42b0fa0086381212ff6651a7bbee1f50ea1164c4",
    "Zhang_10.csv": "7ac2526d41e609500fcdb3b85eb4ed45bd8e01f2",
    "down_syndrome.csv": "e32f413f0e374fd0d64bc03c13aa806c24af72a7",
}


class CoefficientChecksumError(RuntimeError):
    """A downloaded coefficient file does not match its recorded checksum."""


def git_blob_sha1(data: bytes) -> str:
    """Return the git blob SHA-1 of *data* (what ``git hash-object`` prints)."""
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()  # noqa: S324


def package_data_dir() -> Path:
    """Return the package's own ``data/`` directory (may not exist)."""
    from importlib import resources as importlib_resources

    return Path(str(importlib_resources.files("just_biomarkers") / "data"))


def coefficients_cache_dir() -> Path:
    """Return the cache directory for this pinned biolearn commit."""
    override = os.environ.get("BIOMARKERS_CACHE_DIR")
    base = Path(override) if override else Path(platformdirs.user_cache_dir("just-biomarkers"))
    return base / "coefficients" / BIOLEARN_COMMIT[:12]


def _is_valid(path: Path, expected_sha: str) -> bool:
    return path.is_file() and git_blob_sha1(path.read_bytes()) == expected_sha


def download_coefficient(
    filename: str,
    force: bool = False,
    base_url: Optional[str] = None,
) -> Path:
    """Download one coefficient file into the cache and return its path.

    A cached copy that matches its checksum is reused unless *force* is set.
    """
    if filename not in COEFFICIENT_FILES:
        raise KeyError(
            f"Unknown coefficient file '{filename}'. "
            f"Known files: {', '.join(sorted(COEFFICIENT_FILES))}"
        )
    expected = COEFFICIENT_FILES[filename]
    dest = coefficients_cache_dir() / filename
    if not force and _is_valid(dest, expected):
        return dest

    url = f"{base_url or BIOLEARN_DATA_URL}/{filename}"
    with urllib.request.urlopen(url) as response:  # noqa: S310
        data = response.read()
    actual = git_blob_sha1(data)
    if actual != expected:
        raise CoefficientChecksumError(
            f"{filename} from {url} has git blob SHA-1 {actual}, expected {expected}; "
            "the file was not saved."
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    tmp.write_bytes(data)
    tmp.replace(dest)
    return dest


def download_coefficients(
    filenames: Optional[Iterable[str]] = None,
    force: bool = False,
    base_url: Optional[str] = None,
) -> list[Path]:
    """Download every coefficient file (or the named ones); return their paths."""
    names = sorted(COEFFICIENT_FILES) if filenames is None else list(filenames)
    return [download_coefficient(name, force=force, base_url=base_url) for name in names]


def coefficient_path(filename: str) -> Path:
    """Return a local path to *filename*, downloading it if needed.

    Order: the package's ``data/`` directory, then the verified cache, then a
    download from biolearn at the pinned commit.
    """
    bundled = package_data_dir() / filename
    if bundled.is_file():
        return bundled
    return download_coefficient(filename)
