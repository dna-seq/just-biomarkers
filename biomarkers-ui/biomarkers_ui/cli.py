"""CLI for launching the biomarkers Reflex UI."""
from __future__ import annotations

import os
from pathlib import Path


def launch_ui() -> None:
    """Start the Reflex dev server for biomarkers-ui."""
    biomarkers_ui_dir = Path(__file__).resolve().parent.parent
    print(f"Starting Biomarkers UI from {biomarkers_ui_dir} ...")
    os.chdir(biomarkers_ui_dir)
    from reflex.reflex import cli as reflex_cli

    reflex_cli(["run"])
