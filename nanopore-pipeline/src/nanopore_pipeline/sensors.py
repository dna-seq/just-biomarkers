"""Dagster sensors for the nanopore methylation pipeline.

Provides a startup sensor that auto-submits the full pipeline when
assets are unmaterialized.
"""

import time

import dagster as dg


def make_all_sensors(
    full_pipeline_job: "dg.JobDefinition",
    download_job: "dg.JobDefinition",
) -> list[dg.SensorDefinition]:
    """Create all pipeline sensors.

    Returns a list of sensor definitions for registration in Definitions.
    """

    @dg.sensor(
        job=full_pipeline_job,
        default_status=dg.DefaultSensorStatus.RUNNING,
        minimum_interval_seconds=30,
        name="startup_sensor",
        description=(
            "Checks whether key nanopore pipeline assets are materialized. "
            "Submits the full pipeline job if any are missing."
        ),
    )
    def startup_sensor(context: dg.SensorEvaluationContext) -> dg.SensorResult | dg.SkipReason:
        import os

        force_run = os.environ.get("NANOPORE_PIPELINE_FORCE_RUN", "")
        startup_job = os.environ.get("NANOPORE_PIPELINE_STARTUP_JOB", "full_pipeline")
        target_job_name = startup_job if startup_job else "full_pipeline"

        check_keys = [
            dg.AssetKey("synology_file_index"),
            dg.AssetKey("nanopore_bed_files"),
            dg.AssetKey("nanopore_methylation_matrix"),
        ]

        missing = [
            k for k in check_keys
            if context.instance.get_latest_materialization_event(k) is None
        ]

        if not missing and not force_run:
            return dg.SkipReason("All nanopore pipeline assets already materialized.")

        active = context.instance.get_runs(
            filters=dg.RunsFilter(
                job_name=target_job_name,
                statuses=[
                    dg.DagsterRunStatus.STARTED,
                    dg.DagsterRunStatus.NOT_STARTED,
                    dg.DagsterRunStatus.QUEUED,
                ],
            )
        )
        if active:
            return dg.SkipReason(f"Already in progress (run {active[0].run_id[:8]}).")

        last_runs = context.instance.get_runs(
            filters=dg.RunsFilter(job_name=target_job_name), limit=1
        )
        if last_runs and last_runs[0].status == dg.DagsterRunStatus.FAILURE:
            run_key = f"nanopore_startup_retry_{int(time.time())}"
        elif force_run:
            run_key = f"nanopore_force_{int(time.time())}"
        else:
            run_key = "nanopore_startup"

        missing_names = [k.to_user_string() for k in missing]
        reason = f"Missing: {missing_names}" if missing else "Force run requested"
        context.log.info(f"Submitting {target_job_name}: {reason}")

        return dg.SensorResult(
            run_requests=[dg.RunRequest(run_key=run_key, job_name=target_job_name)],
        )

    return [startup_sensor]
