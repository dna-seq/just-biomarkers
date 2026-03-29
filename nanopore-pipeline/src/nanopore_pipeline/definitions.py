"""Dagster Definitions for the nanopore methylation pipeline."""

import dagster as dg
from dagster import in_process_executor

from nanopore_pipeline.assets import (
    nanopore_bed_files,
    nanopore_methylation_matrix,
    synology_file_index,
    synology_nanopore_source,
)
from nanopore_pipeline.resources import CacheDirResource, SourceURLResource
from nanopore_pipeline.sensors import make_all_sensors
from nanopore_pipeline.utils import resource_summary_hook

download_pipeline = dg.define_asset_job(
    name="download_pipeline",
    selection=dg.AssetSelection.assets(
        "synology_file_index", "nanopore_bed_files",
    ),
    description="Download nanopore BED files from the Synology/GoFile source.",
    hooks={resource_summary_hook},
    executor_def=in_process_executor,
)

full_pipeline = dg.define_asset_job(
    name="full_pipeline",
    selection=dg.AssetSelection.assets(
        "synology_file_index",
        "nanopore_bed_files",
        "nanopore_methylation_matrix",
    ),
    description=(
        "Full nanopore pipeline: discover files, download BED data, "
        "and build a unified CpG methylation beta-value matrix."
    ),
    hooks={resource_summary_hook},
    executor_def=in_process_executor,
)

_assets = [
    synology_nanopore_source,
    synology_file_index,
    nanopore_bed_files,
    nanopore_methylation_matrix,
]
_resources = {
    "cache_dir_resource": CacheDirResource(),
    "source_url_resource": SourceURLResource(),
}
_unresolved_jobs = [
    download_pipeline,
    full_pipeline,
]


def _build_definitions() -> dg.Definitions:
    """Resolve unresolved asset jobs and build the final Definitions.

    The temporary Definitions used for resolution is a local variable so
    Dagster's module scanner only finds one Definitions object (the returned one).
    """
    tmp = dg.Definitions(assets=_assets, resources=_resources)
    asset_graph = tmp.resolve_asset_graph()
    resolved_jobs = [
        uj.resolve(asset_graph=asset_graph, resource_defs=_resources)
        for uj in _unresolved_jobs
    ]
    jobs_by_name = {j.name: j for j in resolved_jobs}

    return dg.Definitions(
        assets=_assets,
        sensors=make_all_sensors(
            full_pipeline_job=jobs_by_name["full_pipeline"],
            download_job=jobs_by_name["download_pipeline"],
        ),
        resources=_resources,
        jobs=resolved_jobs,
    )


defs = _build_definitions()


def main() -> None:
    """Launch the Dagster webserver for the nanopore pipeline."""
    import subprocess
    import sys

    subprocess.run(
        [sys.executable, "-m", "dagster", "dev", "-m", "nanopore_pipeline.definitions"],
        check=True,
    )
