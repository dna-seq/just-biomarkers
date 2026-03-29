"""Compare just-biomarkers vs biolearn scores on downloaded EPIC datasets.

biolearn is a dev-only dependency — it is imported lazily and never required
at runtime.  This module is wired into the CLI as ``biomarkers compare-epic``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import polars as pl
from rich.console import Console
from rich import box
from rich.table import Table

console = Console()

EXAMPLES_DIR = Path("data/input/examples")

EPIC_DATASETS = ["GSE112618", "GSE110554", "GSE164056"]

FULL_MATCH_CLOCKS = [
    "PhenoAge",
    "Horvathv2",
    "DunedinPoAm38",
    "DunedinPACE",
    "YingDamAge",
    "YingAdaptAge",
]

PARTIAL_MATCH_CLOCKS = [
    "Horvathv1",
    "Hannum",
    "Lin",
    "Zhang_10",
    "YingCausAge",
]

ALL_CLOCKS = FULL_MATCH_CLOCKS + PARTIAL_MATCH_CLOCKS


def _load_our_dnam(csv_path: Path) -> pl.DataFrame:
    from just_biomarkers.io import read_methylation_csv
    return read_methylation_csv(csv_path)


def _load_biolearn_geo(csv_path: Path):  # type: ignore[return]
    """Load a methylation CSV as a biolearn GeoData object."""
    import pandas as pd
    from biolearn.data_library import GeoData
    df = pd.read_csv(csv_path, index_col=0)
    return GeoData(metadata=pd.DataFrame(index=df.columns), dnam=df)


def _patch_biolearn_readonly_bug() -> None:
    """Patch biolearn's quantile_normalize_using_target for read-only arrays.

    Newer numpy/pandas return read-only views from ``.T``; biolearn tries to
    modify them in-place, which raises ``ValueError``.
    """
    import numpy as np
    import biolearn.dunedin_pace as _dp
    orig = _dp.quantile_normalize_using_target
    if getattr(orig, "_patched_writable", False):
        return

    def _writable(data, target_values):  # type: ignore[no-untyped-def]
        return orig(np.array(data, dtype=np.float64), target_values)

    _writable._patched_writable = True  # type: ignore[attr-defined]
    _dp.quantile_normalize_using_target = _writable


def _score_ours(dnam: pl.DataFrame, clock_name: str) -> dict[str, float]:
    from just_biomarkers.scoring import score_clock
    result = score_clock(dnam, clock_name)
    return {row["sample_id"]: row["score"] for row in result.iter_rows(named=True)}


def _score_biolearn(gallery, geo, clock_name: str) -> dict[str, float] | None:  # type: ignore[no-untyped-def]
    """Return per-sample scores from biolearn, or None if CpGs are missing."""
    model = gallery.get(clock_name, imputation_method="none")
    available = set(geo.dnam.index)
    required = set(model.methylation_sites())
    if required - available:
        return None
    result = model.predict(geo)
    return result["Predicted"].to_dict()


def compare_dataset(
    dataset_id: str,
    clocks: list[str],
    max_samples: Optional[int],
) -> None:
    """Run a side-by-side comparison for one GEO dataset."""
    import pandas as pd

    csv_path = EXAMPLES_DIR / f"{dataset_id}_methylation.csv"
    if not csv_path.exists():
        console.print(
            f"[yellow]Skipping {dataset_id}: file not found at {csv_path}[/yellow]"
        )
        console.print("  Run: [cyan]uv run biomarkers download-epic[/cyan]")
        return

    console.rule(f"[bold cyan]{dataset_id}[/bold cyan]")

    our_dnam = _load_our_dnam(csv_path)
    geo = _load_biolearn_geo(csv_path)

    all_samples = [c for c in our_dnam.columns if c != "CpGmarker"]
    if max_samples is not None:
        samples_to_show = all_samples[:max_samples]
        our_dnam_sub = our_dnam.select(["CpGmarker"] + samples_to_show)
        from biolearn.data_library import GeoData
        geo_sub = GeoData(
            metadata=pd.DataFrame(index=samples_to_show),
            dnam=geo.dnam[samples_to_show],
        )
    else:
        our_dnam_sub = our_dnam
        geo_sub = geo
        samples_to_show = all_samples

    console.print(
        f"Samples: [bold]{len(samples_to_show)}[/bold]  "
        f"CpGs: [bold]{our_dnam.height:,}[/bold]  "
        f"Clocks: [bold]{len(clocks)}[/bold]"
    )

    from biolearn.model_gallery import ModelGallery
    _patch_biolearn_readonly_bug()
    gallery = ModelGallery()

    summary = Table(
        title=f"{dataset_id} — Score Comparison Summary",
        box=box.ROUNDED,
        show_lines=False,
    )
    summary.add_column("Clock", style="cyan", no_wrap=True)
    summary.add_column("CpG Match", justify="right")
    summary.add_column("Biolearn", justify="center")
    summary.add_column("Max |delta|", justify="right")
    summary.add_column("Mean |delta|", justify="right")
    summary.add_column("Status", justify="center")

    for clock_name in clocks:
        our_scores = _score_ours(our_dnam_sub, clock_name)

        from just_biomarkers.scoring import score_clock
        result_full = score_clock(our_dnam_sub, clock_name)
        first_row = result_full.row(0, named=True)
        match_pct = (
            f"{first_row['cpgs_matched']}/{first_row['cpgs_required']} "
            f"({first_row['cpgs_matched']/max(first_row['cpgs_required'],1):.0%})"
        )

        bl_scores = _score_biolearn(gallery, geo_sub, clock_name)

        if bl_scores is None:
            summary.add_row(
                clock_name, match_pct,
                "[yellow]skipped[/yellow]", "—", "—",
                "[yellow]partial match[/yellow]",
            )
            continue

        diffs = [
            abs(our_scores[s] - bl_scores[s])
            for s in samples_to_show
            if s in our_scores and s in bl_scores
        ]
        max_diff = max(diffs) if diffs else float("nan")
        mean_diff = sum(diffs) / len(diffs) if diffs else float("nan")

        status = "[green]OK[/green]" if max_diff < 1e-3 else "[red]DRIFT[/red]"
        summary.add_row(
            clock_name, match_pct,
            "[green]yes[/green]",
            f"{max_diff:.2e}", f"{mean_diff:.2e}",
            status,
        )

    console.print(summary)

    runnable = [
        c for c in clocks
        if _score_biolearn(gallery, geo_sub, c) is not None
    ]
    if not runnable:
        return

    detail_samples = samples_to_show[:5]
    detail = Table(
        title=f"{dataset_id} — Per-sample detail (first {len(detail_samples)} samples, full-match clocks)",
        box=box.SIMPLE_HEAVY,
        show_lines=True,
    )
    detail.add_column("Clock", style="cyan", no_wrap=True)
    detail.add_column("Sample", style="bold")
    detail.add_column("Ours", justify="right")
    detail.add_column("Biolearn", justify="right")
    detail.add_column("|delta|", justify="right")

    for clock_name in runnable:
        our_scores = _score_ours(our_dnam_sub, clock_name)
        bl_scores = _score_biolearn(gallery, geo_sub, clock_name)
        for sample in detail_samples:
            ours = our_scores.get(sample)
            bl = bl_scores.get(sample) if bl_scores else None
            if ours is None or bl is None:
                continue
            diff = abs(ours - bl)
            diff_color = "red" if diff >= 1e-3 else "green"
            detail.add_row(
                clock_name, sample,
                f"{ours:.4f}", f"{bl:.4f}",
                f"[{diff_color}]{diff:.2e}[/{diff_color}]",
            )

    console.print(detail)
    console.print()
