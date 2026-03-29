"""Dagster resources for the nanopore methylation pipeline."""

import os
from pathlib import Path

from dagster import ConfigurableResource
from dotenv import load_dotenv
from platformdirs import user_cache_dir


class CacheDirResource(ConfigurableResource):
    """Provides a configurable root cache directory for nanopore data.

    Resolution order:
    1. Explicit ``cache_dir`` config value
    2. ``BIOMARKERS_CACHE_DIR`` environment variable
    3. OS-appropriate user cache dir via platformdirs
    """

    cache_dir: str = ""

    def get_path(self) -> Path:
        """Return the resolved cache directory Path."""
        if self.cache_dir:
            return Path(self.cache_dir)
        raw = os.environ.get("BIOMARKERS_CACHE_DIR", "")
        if raw:
            return Path(raw)
        return Path(user_cache_dir("just-biomarkers"))

    def nanopore_dir(self) -> Path:
        """Return the nanopore-specific subdirectory."""
        d = self.get_path() / "nanopore"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def downloads_dir(self) -> Path:
        """Return the directory for raw downloaded BED files."""
        d = self.nanopore_dir() / "downloads"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def processed_dir(self) -> Path:
        """Return the directory for processed output files."""
        d = self.nanopore_dir() / "processed"
        d.mkdir(parents=True, exist_ok=True)
        return d


class SourceURLResource(ConfigurableResource):
    """Provides the configurable source URL for nanopore data.

    The Synology/GoFile URL serves a directory listing of .bed.gz files.
    """

    source_url: str = ""

    def get_url(self) -> str:
        """Return the resolved source URL."""
        if self.source_url:
            return self.source_url
        load_dotenv()
        return os.environ.get(
            "NANOPORE_SOURCE_URL",
            "http://gofile.me/76Mv5/irFIggw48",
        )
