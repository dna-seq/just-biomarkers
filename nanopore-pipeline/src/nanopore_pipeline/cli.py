"""CLI for nanopore-pipeline: launch or run the Dagster nanopore methylation pipeline."""

import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Annotated, Optional

import typer
from dotenv import load_dotenv
from rich.console import Console

load_dotenv()

app = typer.Typer(help="Nanopore methylation data pipeline (Dagster)")
console = Console()

_DEFAULT_HOST = os.environ.get("NANOPORE_PIPELINE_HOST", "0.0.0.0")
_DEFAULT_PORT = int(os.environ.get("NANOPORE_PIPELINE_PORT", "3020"))
_GRPC_CLEANUP_TIMEOUT = 5


def _kill_stale_dagster_children() -> None:
    """SIGKILL any leftover dagster gRPC server child processes."""
    import psutil

    current = psutil.Process()
    for child in current.children(recursive=True):
        cmdline = " ".join(child.cmdline())
        if "dagster" in cmdline and "grpc" in cmdline:
            console.print(f"[dim]Cleaning up stale gRPC subprocess (PID {child.pid})...[/dim]")
            child.kill()
            child.wait(timeout=_GRPC_CLEANUP_TIMEOUT)


def _execute_job(resolved_job: "dagster.JobDefinition") -> "dagster.ExecuteInProcessResult":  # type: ignore[name-defined]
    """Run a Dagster job in-process, handling gRPC server cleanup timeouts."""
    from dagster import DagsterInstance

    with DagsterInstance.get() as instance:
        result = resolved_job.execute_in_process(instance=instance)

    _kill_stale_dagster_children()
    return result


def _find_project_root() -> Path:
    """Walk upward from cwd to find the uv workspace root (contains [tool.uv.workspace])."""
    current = Path.cwd().resolve()
    for candidate in [current, *current.parents]:
        pyproject = candidate / "pyproject.toml"
        if pyproject.exists() and "[tool.uv.workspace]" in pyproject.read_text():
            return candidate
    return current


def _setup_dagster_home() -> Path:
    """Set DAGSTER_HOME to a project-relative path and create dagster.yaml."""
    project_root = _find_project_root()
    dagster_home = project_root / "data" / "output" / "dagster-nanopore"
    dagster_home.mkdir(parents=True, exist_ok=True)
    os.environ["DAGSTER_HOME"] = str(dagster_home)

    yaml_path = dagster_home / "dagster.yaml"
    if not yaml_path.exists():
        yaml_path.write_text("telemetry:\n  enabled: false\n")
        console.print(f"[dim]Created {yaml_path}[/dim]")

    return dagster_home


def _kill_port(port: int) -> None:
    """Kill any process listening on the given TCP port."""
    result = subprocess.run(
        ["lsof", "-t", f"-iTCP:{port}"],
        capture_output=True, text=True,
    )
    pids = [int(p) for p in result.stdout.strip().splitlines() if p.strip()]
    if pids:
        console.print(f"[yellow]Port {port} in use by PIDs {pids}, terminating...[/yellow]")
        for pid in pids:
            os.kill(pid, signal.SIGTERM)
        time.sleep(1)
        result2 = subprocess.run(
            ["lsof", "-t", f"-iTCP:{port}"],
            capture_output=True, text=True,
        )
        for pid in [int(p) for p in result2.stdout.strip().splitlines() if p.strip()]:
            os.kill(pid, signal.SIGKILL)


def _cancel_orphaned_runs() -> None:
    """Cancel any runs stuck in STARTED or NOT_STARTED from a previous session."""
    from dagster import DagsterInstance, DagsterRunStatus, RunsFilter

    with DagsterInstance.get() as instance:
        stuck = instance.get_run_records(
            filters=RunsFilter(statuses=[DagsterRunStatus.STARTED, DagsterRunStatus.NOT_STARTED])
        )
        for record in stuck:
            console.print(
                f"[yellow]Cancelling orphaned run {record.dagster_run.run_id[:8]}...[/yellow]"
            )
            instance.report_run_canceled(
                record.dagster_run,
                message="Orphaned run from previous session",
            )


@app.command()
def run(
    job: Annotated[str, typer.Option(
        help="Job to run: full_pipeline, download_pipeline."
    )] = "full_pipeline",
    headless: Annotated[bool, typer.Option(
        "--headless", help="Run in-process without Dagster UI."
    )] = False,
    host: Annotated[str, typer.Option(
        help="Bind address for the Dagster webserver (UI mode only)."
    )] = _DEFAULT_HOST,
    port: Annotated[int, typer.Option(
        help="Port for the Dagster webserver (UI mode only)."
    )] = _DEFAULT_PORT,
) -> None:
    """Run the nanopore pipeline with Dagster UI by default.

    \b
    Default behavior launches Dagster UI so you can monitor execution
    live. The startup sensor submits the selected job automatically.

    Use ``--headless`` to run in-process without UI.
    """
    dagster_home = _setup_dagster_home()
    os.environ["NANOPORE_PIPELINE_STARTUP_JOB"] = job
    os.environ["NANOPORE_PIPELINE_FORCE_RUN"] = "1"

    if headless:
        _cancel_orphaned_runs()
        console.print(f"[dim]DAGSTER_HOME={dagster_home}[/dim]")
        console.print("[yellow]HEADLESS MODE — no Dagster UI.[/yellow]")

        from nanopore_pipeline.definitions import defs

        resolved_job = defs.get_job_def(job)
        console.print(f"\n[bold]Running job: {job}[/bold]")
        console.print(f"[dim]{resolved_job.description or ''}[/dim]\n")

        result = _execute_job(resolved_job)

        if result.success:
            console.print(f"\n[green bold]Job '{job}' completed successfully.[/green bold]")
        else:
            console.print(f"\n[red bold]Job '{job}' failed.[/red bold]")
            for event in result.all_events:
                if event.is_failure:
                    console.print(f"  [red]{event.message}[/red]")
            raise typer.Exit(code=1)
        return

    _kill_port(port)
    _cancel_orphaned_runs()
    console.print(f"[dim]DAGSTER_HOME={dagster_home}[/dim]")
    console.print(f"[bold green]Dagster UI:[/bold green] http://{host}:{port}")
    console.print(f"[bold]Job '{job}' will be submitted automatically on startup.[/bold]\n")

    dagster_bin = str(Path(sys.executable).parent / "dagster")
    os.execvp(dagster_bin, [
        "dagster", "dev",
        "-m", "nanopore_pipeline.definitions",
        "--host", host,
        "--port", str(port),
    ])


@app.command()
def launch(
    host: Annotated[str, typer.Option(
        help="Bind address for the Dagster webserver."
    )] = _DEFAULT_HOST,
    port: Annotated[int, typer.Option(
        help="Port for the Dagster webserver."
    )] = _DEFAULT_PORT,
) -> None:
    """Start the Dagster UI webserver without pre-selecting a job.

    \b
    Launches the Dagster dev server with full monitoring UI.
    The startup sensor will automatically submit the full_pipeline job
    if any assets are unmaterialized. You can trigger jobs manually from the UI.
    """
    dagster_home = _setup_dagster_home()
    _kill_port(port)
    _cancel_orphaned_runs()

    console.print(f"[dim]DAGSTER_HOME={dagster_home}[/dim]")
    console.print(f"[bold green]Dagster UI:[/bold green] http://{host}:{port}")
    console.print("[dim]The startup sensor will submit jobs only if assets are missing.[/dim]")
    console.print("[dim]Use the UI to trigger jobs, or run 'nanopore run' in another terminal.[/dim]\n")

    dagster_bin = str(Path(sys.executable).parent / "dagster")
    os.execvp(dagster_bin, [
        "dagster", "dev",
        "-m", "nanopore_pipeline.definitions",
        "--host", host,
        "--port", str(port),
    ])


@app.command()
def status(
    cache_dir: Annotated[Optional[str], typer.Option(
        "--cache-dir", help="Override cache directory."
    )] = None,
) -> None:
    """Show the status of downloaded and processed nanopore data."""
    import polars as pl
    from nanopore_pipeline.resources import CacheDirResource

    resource = CacheDirResource(cache_dir=cache_dir or "")
    nanopore_dir = resource.nanopore_dir()
    downloads_dir = resource.downloads_dir()
    processed_dir = resource.processed_dir()

    console.print(f"\n[bold]Nanopore pipeline status[/bold]")
    console.print(f"  Cache dir:      {nanopore_dir}")
    console.print(f"  Downloads dir:  {downloads_dir}")
    console.print(f"  Processed dir:  {processed_dir}")

    index_path = nanopore_dir / "file_index.parquet"
    if index_path.exists():
        index_df = pl.read_parquet(index_path)
        console.print(f"\n[bold]File index:[/bold] {index_df.height} files discovered")
        for row in index_df.head(20).iter_rows(named=True):
            console.print(f"  {row['name']}")
        if index_df.height > 20:
            console.print(f"  [dim]... and {index_df.height - 20} more[/dim]")
    else:
        console.print("\n[yellow]No file index found. Run 'nanopore run' first.[/yellow]")

    bed_files = list(downloads_dir.glob("*.bed.gz")) + list(downloads_dir.glob("*.bed"))
    console.print(f"\n[bold]Downloaded BED files:[/bold] {len(bed_files)}")
    for f in sorted(bed_files)[:20]:
        size_mb = f.stat().st_size / (1024 * 1024)
        console.print(f"  {f.name} ({size_mb:.1f} MB)")
    if len(bed_files) > 20:
        console.print(f"  [dim]... and {len(bed_files) - 20} more[/dim]")

    matrix_path = processed_dir / "nanopore_methylation_matrix.parquet"
    if matrix_path.exists():
        matrix_df = pl.read_parquet(matrix_path, n_rows=0)
        n_samples = len(matrix_df.columns) - 1
        schema = matrix_df.schema
        console.print(f"\n[bold]Methylation matrix:[/bold]")
        console.print(f"  Path:     {matrix_path}")
        console.print(f"  Samples:  {n_samples}")
        console.print(f"  Columns:  {list(matrix_df.columns[:10])}")
        if len(matrix_df.columns) > 10:
            console.print(f"            [dim]... and {len(matrix_df.columns) - 10} more[/dim]")

        full_df = pl.scan_parquet(matrix_path)
        row_count = full_df.select(pl.len()).collect().item()
        console.print(f"  CpG sites: {row_count:,}")
    else:
        console.print("\n[yellow]No methylation matrix found. Run the full pipeline first.[/yellow]")


@app.command()
def clean(
    dry_run: Annotated[bool, typer.Option(
        "--dry-run", help="Show what would be cleaned without doing it."
    )] = False,
) -> None:
    """Cancel stuck runs from the Dagster DB."""
    from dagster import DagsterInstance, DagsterRunStatus, RunsFilter

    _setup_dagster_home()

    with DagsterInstance.get() as instance:
        cancel_statuses = [
            DagsterRunStatus.QUEUED,
            DagsterRunStatus.NOT_STARTED,
            DagsterRunStatus.STARTED,
        ]
        cancelled = 0
        for status_val in cancel_statuses:
            records = instance.get_run_records(
                filters=RunsFilter(statuses=[status_val]),
            )
            for rec in records:
                if not dry_run:
                    instance.report_run_canceled(
                        rec.dagster_run,
                        message="Cancelled by nanopore pipeline clean command",
                    )
                cancelled += 1

        if cancelled:
            label = "Would cancel" if dry_run else "Cancelled"
            console.print(f"[yellow]{label} {cancelled} queued/in-progress runs.[/yellow]")
        else:
            console.print("[green]No stuck runs found.[/green]")

        if dry_run:
            console.print("\n[dim]Run without --dry-run to execute cleanup.[/dim]")
