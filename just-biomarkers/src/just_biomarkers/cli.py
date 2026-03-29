"""Typer CLI for just-biomarkers methylation clock computation."""
from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(
    name="just-biomarkers",
    help="Polars-based methylation clock computation",
    no_args_is_help=True,
)
console = Console()


@app.command()
def list_clocks(
    tissue: Optional[str] = typer.Option(None, help="Filter by tissue"),
    output: Optional[str] = typer.Option(None, help="Filter by output type"),
) -> None:
    """List available methylation clocks."""
    from just_biomarkers.registry import CLOCK_DEFINITIONS, search_clocks

    if tissue or output:
        clocks = search_clocks(tissue=tissue, output=output)
    else:
        clocks = CLOCK_DEFINITIONS

    table = Table(title=f"Available Clocks ({len(clocks)})")
    table.add_column("Name", style="cyan")
    table.add_column("Year", justify="right")
    table.add_column("Tissue", style="green")
    table.add_column("Output", style="yellow")

    for name, defn in sorted(clocks.items()):
        table.add_row(name, str(defn.year), defn.tissue, defn.output)

    console.print(table)


@app.command()
def compute(
    methylation_input: Path = typer.Argument(
        ...,
        help=(
            "Path to methylation matrix file "
            "(CSV/TSV/TXT/Parquet or IDAT directory/file pair)"
        ),
    ),
    clocks: Optional[str] = typer.Option(
        None,
        help="Comma-separated clock names (default: all linear clocks)",
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output CSV path for results"
    ),
    idat_preprocess: Literal["raw", "noob"] = typer.Option(
        "raw",
        help="IDAT preprocessing mode (used only for IDAT inputs).",
    ),
) -> None:
    """Compute methylation clock scores from a methylation matrix file."""
    from just_biomarkers.io import read_methylation_matrix, validate_methylation_matrix
    from just_biomarkers.scoring import score_clocks

    console.print(f"Reading methylation matrix from [cyan]{methylation_input}[/cyan]")
    try:
        dnam = read_methylation_matrix(
            methylation_input, idat_preprocess=idat_preprocess
        )
    except (ImportError, ValueError) as exc:
        raise typer.BadParameter(str(exc), param_hint="methylation_input") from exc

    warnings = validate_methylation_matrix(dnam)
    for w in warnings:
        console.print(f"[yellow]Warning:[/yellow] {w}")

    clock_names = None
    if clocks:
        clock_names = [c.strip() for c in clocks.split(",")]

    console.print("Scoring clocks...")
    batch = score_clocks(dnam, clock_names=clock_names)

    for w in batch.warnings:
        console.print(f"[yellow]Warning:[/yellow] {w}")

    table = Table(title="Clock Scores")
    table.add_column("Sample", style="cyan")
    table.add_column("Clock", style="green")
    table.add_column("Score", justify="right", style="bold")
    table.add_column("Match Rate", justify="right")

    for r in batch.results:
        table.add_row(
            r.sample_id,
            r.clock_name,
            f"{r.score:.4f}",
            f"{r.match_rate:.1%}",
        )

    console.print(table)

    if output:
        import polars as pl

        results_df = pl.DataFrame(
            [
                {
                    "sample_id": r.sample_id,
                    "clock_name": r.clock_name,
                    "score": r.score,
                    "cpgs_matched": r.cpgs_matched,
                    "cpgs_required": r.cpgs_required,
                    "match_rate": r.match_rate,
                }
                for r in batch.results
            ]
        )
        results_df.write_csv(str(output))
        console.print(f"Results saved to [green]{output}[/green]")


def _print_dataset_table(datasets: dict[str, dict[str, str]], title: str) -> None:
    table = Table(title=title)
    table.add_column("ID", style="cyan")
    table.add_column("Format", style="green")
    table.add_column("Samples", justify="right")
    table.add_column("Tissue", style="yellow")
    table.add_column("Description")
    for ds_id, info in datasets.items():
        table.add_row(
            ds_id,
            info["format"],
            info["samples"],
            info["tissue"],
            info["description"],
        )
    console.print(table)


@app.command("download-example")
def download_example(
    dataset_id: str = typer.Argument(
        "GSE112618",
        help="GEO dataset ID to download (e.g. GSE112618, GSE110554, GSE40279).",
    ),
    output_dir: Path = typer.Option(
        Path("data/input/examples"),
        "--output-dir",
        "-o",
        help="Directory to write the output CSV file.",
    ),
    max_samples: Optional[int] = typer.Option(
        None,
        "--max-samples",
        help="Keep only the first N sample columns (handy for large datasets).",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Re-download and re-parse even if the output CSV already exists.",
    ),
    list_datasets: bool = typer.Option(
        False,
        "--list",
        "-l",
        help="List available datasets and exit.",
    ),
) -> None:
    """Download a single GEO methylation example dataset as a CpGmarker CSV."""
    from just_biomarkers.geo_download import GEO_EXAMPLE_DATASETS, download_geo_example

    if list_datasets:
        _print_dataset_table(GEO_EXAMPLE_DATASETS, "Available GEO Example Datasets")
        return

    console.print(f"Downloading [cyan]{dataset_id}[/cyan] …")
    dest = download_geo_example(
        dataset_id,
        output_dir,
        max_samples=max_samples,
        force=force,
    )
    console.print(f"[green]Saved to:[/green] {dest}")


@app.command("download-epic")
def download_epic(
    output_dir: Path = typer.Option(
        Path("data/input/examples"),
        "--output-dir",
        "-o",
        help="Directory to write CSV files into.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Re-download and re-parse even if output CSVs already exist.",
    ),
    list_datasets: bool = typer.Option(
        False,
        "--list",
        "-l",
        help="List EPIC datasets that would be downloaded, then exit.",
    ),
) -> None:
    """Download all IlluminaEPIC example datasets to data/input/examples/."""
    from just_biomarkers.geo_download import GEO_EXAMPLE_DATASETS, download_geo_example

    epic_datasets = {
        k: v for k, v in GEO_EXAMPLE_DATASETS.items() if v["format"] == "IlluminaEPIC"
    }

    if list_datasets:
        _print_dataset_table(epic_datasets, "IlluminaEPIC Datasets")
        return

    console.print(
        f"Downloading [bold]{len(epic_datasets)}[/bold] EPIC datasets "
        f"→ [cyan]{output_dir}[/cyan]"
    )

    results: list[tuple[str, Path | Exception]] = []
    for ds_id, info in epic_datasets.items():
        console.print(f"  [{ds_id}] {info['description']} …", end=" ")
        dest = download_geo_example(ds_id, output_dir, force=force)
        console.print(f"[green]✓[/green] {dest.name}")
        results.append((ds_id, dest))

    console.print(f"\n[bold green]Done.[/bold green] {len(results)} files in [cyan]{output_dir}[/cyan]:")
    for ds_id, dest in results:
        console.print(f"  {dest}")


@app.command("download-all-examples")
def download_all_examples(
    output_dir: Path = typer.Option(
        Path("data/input/examples"),
        "--output-dir",
        "-o",
        help="Directory to write CSV files into.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Re-download and re-parse even if output CSVs already exist.",
    ),
    list_datasets: bool = typer.Option(
        False,
        "--list",
        "-l",
        help="List all datasets that would be downloaded, then exit.",
    ),
) -> None:
    """Download every GEO example dataset (EPIC + 450k) to data/input/examples/."""
    from just_biomarkers.geo_download import GEO_EXAMPLE_DATASETS, download_geo_example

    if list_datasets:
        _print_dataset_table(GEO_EXAMPLE_DATASETS, "All GEO Example Datasets")
        return

    console.print(
        f"Downloading [bold]{len(GEO_EXAMPLE_DATASETS)}[/bold] datasets "
        f"→ [cyan]{output_dir}[/cyan]"
    )

    results: list[tuple[str, Path]] = []
    for ds_id, info in GEO_EXAMPLE_DATASETS.items():
        console.print(f"  [{ds_id}] {info['description']} …", end=" ")
        dest = download_geo_example(ds_id, output_dir, force=force)
        console.print(f"[green]✓[/green] {dest.name}")
        results.append((ds_id, dest))

    console.print(f"\n[bold green]Done.[/bold green] {len(results)} files in [cyan]{output_dir}[/cyan]:")
    for ds_id, dest in results:
        console.print(f"  {dest}")


@app.command("download-computage")
def download_computage(
    dataset_id: Optional[str] = typer.Argument(
        None,
        help=(
            "GEO dataset ID from ComputAgeBench (e.g. GSE100264). "
            "If omitted, downloads metadata only."
        ),
    ),
    output_dir: Path = typer.Option(
        Path("data/input/computage"),
        "--output-dir",
        "-o",
        help="Directory to write the output files.",
    ),
    split: str = typer.Option(
        "benchmark",
        "--split",
        "-s",
        help="ComputAgeBench split: 'benchmark' or 'train'.",
    ),
    max_samples: Optional[int] = typer.Option(
        None,
        "--max-samples",
        help="Keep only the first N sample columns.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Re-download and re-convert even if the output already exists.",
    ),
    list_datasets: bool = typer.Option(
        False,
        "--list",
        "-l",
        help="List available datasets in the split and exit.",
    ),
    meta_only: bool = typer.Option(
        False,
        "--meta",
        help="Download only the metadata TSV for the split.",
    ),
) -> None:
    """Download ComputAgeBench methylation data from Hugging Face.

    Downloads individual studies or metadata from the ComputAgeBench
    epigenetic aging clocks benchmark (computage/computage_bench on HF).
    """
    from just_biomarkers.computage_download import (
        download_computage_dataset,
        download_computage_meta,
        list_computage_datasets,
        load_computage_meta,
    )

    if list_datasets:
        console.print(f"Fetching ComputAgeBench [cyan]{split}[/cyan] dataset list ...")
        ids = list_computage_datasets(split=split)
        table = Table(title=f"ComputAgeBench {split} datasets ({len(ids)})")
        table.add_column("Dataset ID", style="cyan")
        for ds_id in ids:
            table.add_row(ds_id)
        console.print(table)
        return

    if meta_only or dataset_id is None:
        console.print(
            f"Downloading ComputAgeBench [cyan]{split}[/cyan] metadata ..."
        )
        dest = download_computage_meta(output_dir, split=split, force=force)
        console.print(f"[green]Metadata saved to:[/green] {dest}")

        meta = load_computage_meta(split=split)
        console.print(
            f"  {meta.height} samples, "
            f"{meta['DatasetID'].n_unique()} studies"
        )
        return

    console.print(
        f"Downloading ComputAgeBench [cyan]{split}/{dataset_id}[/cyan] ..."
    )
    dest = download_computage_dataset(
        dataset_id,
        output_dir,
        split=split,
        max_samples=max_samples,
        force=force,
    )
    console.print(f"[green]Saved to:[/green] {dest}")


@app.command("download-computage-all")
def download_computage_all(
    output_dir: Path = typer.Option(
        Path("data/input/computage"),
        "--output-dir",
        "-o",
        help="Directory to write the output CSV files.",
    ),
    split: str = typer.Option(
        "benchmark",
        "--split",
        "-s",
        help="ComputAgeBench split: 'benchmark' or 'train'.",
    ),
    max_samples: Optional[int] = typer.Option(
        None,
        "--max-samples",
        help="Keep only the first N sample columns per dataset.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Re-download and re-convert even if outputs already exist.",
    ),
) -> None:
    """Download ALL ComputAgeBench datasets for a split as CpGmarker CSVs."""
    from just_biomarkers.computage_download import (
        download_computage_dataset,
        download_computage_meta,
        list_computage_datasets,
    )

    console.print(
        f"Downloading all ComputAgeBench [cyan]{split}[/cyan] datasets ..."
    )

    meta_dest = download_computage_meta(output_dir, split=split, force=force)
    console.print(f"  [green]Metadata:[/green] {meta_dest}")

    ids = list_computage_datasets(split=split)
    console.print(f"  [bold]{len(ids)}[/bold] datasets to process")

    results: list[tuple[str, Path]] = []
    for ds_id in ids:
        console.print(f"  [{ds_id}] ...", end=" ")
        dest = download_computage_dataset(
            ds_id,
            output_dir,
            split=split,
            max_samples=max_samples,
            force=force,
        )
        console.print(f"[green]OK[/green] {dest.name}")
        results.append((ds_id, dest))

    console.print(
        f"\n[bold green]Done.[/bold green] "
        f"{len(results)} files in [cyan]{output_dir}[/cyan]"
    )


@app.command("compare-epic")
def compare_epic(
    dataset: Optional[str] = typer.Option(
        None,
        "--dataset",
        "-d",
        help="Single dataset ID to compare (default: all EPIC datasets).",
    ),
    clocks: Optional[str] = typer.Option(
        None,
        "--clocks",
        "-c",
        help="Comma-separated clock names (default: all supported).",
    ),
    max_samples: Optional[int] = typer.Option(
        None,
        "--max-samples",
        "-n",
        help="Limit to first N samples per dataset (faster for exploration).",
    ),
    full_match_only: bool = typer.Option(
        False,
        "--full-match-only",
        help="Only run clocks with 100%% CpG coverage on EPIC.",
    ),
) -> None:
    """Compare just-biomarkers vs biolearn on downloaded EPIC datasets.

    Requires biolearn (dev dependency). Download datasets first with:
        uv run biomarkers download-epic
    """
    from just_biomarkers.compare import (
        ALL_CLOCKS,
        EPIC_DATASETS,
        FULL_MATCH_CLOCKS,
        compare_dataset,
    )

    datasets = [dataset] if dataset else EPIC_DATASETS

    if clocks:
        clock_list = [c.strip() for c in clocks.split(",")]
    elif full_match_only:
        clock_list = FULL_MATCH_CLOCKS
    else:
        clock_list = ALL_CLOCKS

    console.print(
        f"\n[bold]just-biomarkers vs biolearn[/bold] — EPIC parity comparison\n"
        f"Datasets: {datasets}\n"
        f"Clocks:   {clock_list}\n"
    )

    for ds in datasets:
        compare_dataset(ds, clock_list, max_samples)

    console.print("[bold green]Comparison complete.[/bold green]")


def run() -> None:
    """Entry point for the CLI."""
    app()
