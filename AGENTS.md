# Python Coding Standards & Best Practices

## Project Architecture

This project is a **uv workspace** with a non-published root wrapper and three subprojects:

- **`just-biomarkers/src/just_biomarkers/`** — Core library: Polars-based methylation clock computation engine. Supports linear clocks (Horvath, Hannum, PhenoAge, DunedinPoAm38, etc.), imputation, and I/O for methylation beta-value matrices. CLI entrypoint via Typer. Published to PyPI as `just-biomarkers`.
- **`biomarkers-ui/`** — Reflex web app for interactive methylation clock computation: upload methylation CSV, select clocks, compute scores, download results. Has its own `pyproject.toml` and depends on `just_biomarkers`. Run with `uv run ui` from workspace root. Published to PyPI as `biomarkers-ui`.
- **`nanopore-pipeline/`** — Dagster pipeline for downloading and processing nanopore methylation data. Discovers files from a Synology/GoFile source, downloads BED files, and builds a unified CpG methylation beta-value matrix. Has its own `pyproject.toml` and depends on `just_biomarkers`. All pipeline commands (`run`, `launch`) default to launching the Dagster UI for monitoring. Use `--headless` on `run` for in-process execution without UI.

The workspace root (`pyproject.toml` at repo root) is a non-published wrapper named `just-biomarkers-workspace`. It depends on all three subprojects and **must re-export all CLI entry points** from subprojects so that every command is available via `uv run <name>` from the workspace root.

**CRITICAL: All subproject CLI entry points must be registered in the workspace root `pyproject.toml` `[project.scripts]`.** When adding a new CLI entry point to any subproject, always add it to the root `pyproject.toml` as well.

| Script | Entry point | Subproject |
|--------|-------------|------------|
| `biomarkers` | `just_biomarkers.cli:app` | `just-biomarkers` |
| `just-biomarkers` | `just_biomarkers.cli:app` | `just-biomarkers` |
| `ui` | `biomarkers_ui.cli:launch_ui` | `biomarkers-ui` |
| `nanopore` | `nanopore_pipeline.cli:app` | `nanopore-pipeline` |
| `nanopore-pipeline` | `nanopore_pipeline.definitions:main` | `nanopore-pipeline` |

### Key modules

| Module | Purpose |
|--------|---------|
| `just_biomarkers.scoring` | Core scoring engine: `score_clock()` and `score_clocks()`. Loads coefficients, joins with methylation matrix, weighted sum, optional transform. |
| `just_biomarkers.registry` | Clock registry: `CLOCK_DEFINITIONS` dict of `ClockDefinition` models. `list_clocks()`, `get_clock()`, `search_clocks()`. |
| `just_biomarkers.models` | Pydantic v2 models: `ClockDefinition`, `ClockResult`, `BatchClockResult`. |
| `just_biomarkers.transforms` | Named transform functions (Horvath anti-trafo, sigmoid, offsets). `TRANSFORM_REGISTRY` maps string keys to callables. |
| `just_biomarkers.imputation` | Imputation strategies: `impute_from_average()`, `impute_from_reference()`, `hybrid_impute()`. All operate on Polars DataFrames. |
| `just_biomarkers.io` | I/O helpers: `read_methylation_csv()`, `read_methylation_parquet()`, `validate_methylation_matrix()`. |
| `just_biomarkers.data_utils` | Package data loading: `get_data_dir()`, `load_coefficients()`. |
| `just_biomarkers.cli` | Typer CLI: `list-clocks`, `compute`, `download-example`, `download-epic`, `download-all-examples`, `download-computage`, `download-computage-all`, `compare-epic`. |
| `just_biomarkers.compare` | EPIC parity comparison engine: `compare_dataset()`, clock/dataset constants. Compares just-biomarkers scores vs biolearn (dev dependency) side-by-side. |
| `just_biomarkers.geo_download` | GEO downloader: `GEO_EXAMPLE_DATASETS` catalogue, `download_geo_example()`. Raw `.gz` files cached under `platformdirs.user_cache_dir("just-biomarkers")/geo_downloads/`. Pure stdlib + polars — no biolearn dependency. |
| `just_biomarkers.computage_download` | ComputAgeBench downloader: `list_computage_datasets()`, `load_computage_dataset()`, `load_computage_meta()`, `download_computage_dataset()`, `download_computage_meta()`. Uses `huggingface_hub` to download from `computage/computage_bench` HF dataset. Cached under `platformdirs.user_cache_dir("just-biomarkers")/computage_bench/`. |
| `biomarkers_ui.state` | Reflex `AppState`: upload handling, clock selection, compute dispatch, results display, GEO example loading. |
| `biomarkers_ui.biomarkers_ui` | Reflex app: upload section, clock selector, compute button, results table. |
| `biomarkers_ui.cli` | UI launcher: `launch_ui()` starts Reflex dev server. |
| `nanopore_pipeline.assets` | Dagster SDAs: `synology_file_index`, `nanopore_bed_files`, `nanopore_methylation_matrix`. SourceAsset: `synology_nanopore_source`. |
| `nanopore_pipeline.definitions` | Dagster `Definitions`: jobs (`full_pipeline`, `download_pipeline`), resources, sensors. `defs` built via `_build_definitions()`. |
| `nanopore_pipeline.cli` | Typer CLI: `run`, `launch`, `status`, `clean`. Dagster UI on port 3020 by default. |
| `nanopore_pipeline.resources` | `CacheDirResource` (cache dir via `BIOMARKERS_CACHE_DIR` / platformdirs), `SourceURLResource` (GoFile URL). |
| `nanopore_pipeline.sensors` | Startup sensor: auto-submits pipeline when assets are unmaterialized. |
| `nanopore_pipeline.runtime` | `ResourceReport` model, `resource_tracker` context manager (CPU/RAM/duration via psutil). |
| `nanopore_pipeline.utils` | `resource_summary_hook` — Dagster `@success_hook` for run-level resource aggregation. |

### Methylation data standard

Following the [biolearn methylation standard](https://github.com/bio-learn/biolearn):
- **Matrix orientation**: rows = CpG site IDs (e.g. `cg00075967`), columns = sample IDs
- **Values**: beta values in `[0, 1]` or NaN
- **Internal column name**: `CpGmarker` (first column of DataFrame)
- **Coefficient files**: CSV with `CpGmarker` (index) and `CoefficientTraining` or `Weight` column

### Running the UI

```bash
uv sync --all-packages
uv run ui
```

### Computing clocks via CLI

```bash
# List available clocks
uv run biomarkers list-clocks
uv run biomarkers list-clocks --tissue Blood

# Compute clocks from a methylation CSV
uv run biomarkers compute methylation_matrix.csv --clocks Horvathv1,Hannum,PhenoAge
uv run biomarkers compute methylation_matrix.csv -o results.csv

# Compare just-biomarkers vs biolearn on EPIC datasets (requires biolearn dev dep)
uv run biomarkers compare-epic
uv run biomarkers compare-epic --dataset GSE112618 --max-samples 10
uv run biomarkers compare-epic --clocks PhenoAge,Horvathv2,DunedinPACE
```

### Running the nanopore pipeline

```bash
uv sync --all-packages

# Launch Dagster UI (default, port 3020)
uv run nanopore run

# Launch Dagster UI only (no job pre-selected)
uv run nanopore launch

# Run headless (no UI)
uv run nanopore run --headless

# Run only the download step
uv run nanopore run --job download_pipeline

# Check status of downloads and processed data
uv run nanopore status
```

Environment variables: `NANOPORE_PIPELINE_HOST` (default `0.0.0.0`), `NANOPORE_PIPELINE_PORT` (default `3020`), `BIOMARKERS_CACHE_DIR`, `NANOPORE_SOURCE_URL` (default `http://gofile.me/76Mv5/irFIggw48`).

### Python API (for downstream reuse in just-dna-lite)

```python
from just_biomarkers import (
    score_clock, score_clocks,
    list_clocks, get_clock, search_clocks,
    read_methylation_csv,
    ClockDefinition, ClockResult, BatchClockResult,
    __version__,
)

dnam = read_methylation_csv("methylation_matrix.csv")
result = score_clocks(dnam, clock_names=["Horvathv1", "Hannum", "PhenoAge"])
for r in result.results:
    print(f"{r.sample_id}: {r.clock_name} = {r.score:.4f} ({r.match_rate:.1%})")
```

### Reference projects (READ-ONLY -- do not edit)

- **`just-prs`** (`/home/antonkulaga/sources/just-prs/`) — Workspace layout, CLI patterns, Reflex UI architecture, and PyPI publishing conventions are modeled after this project.
- **`biolearn`** (`/home/antonkulaga/sources/biolearn/`) — Clock definitions, coefficient data, and methylation standard are derived from this project. `biolearn` is used as a **dev-only parity oracle** in integration tests (not a runtime dependency).
- **`biolearn` must NEVER be added as a runtime dependency.** It is heavy (pandas, matplotlib, seaborn). Any functionality that could be borrowed from biolearn must be re-implemented using stdlib + polars + platformdirs instead.

---

## Data Directory Conventions

**Data must be strictly separated from code.** Generated data, downloaded files, uploaded files, and computation outputs must NEVER be written to the project root or source tree.

### Input data (`data/input/`)

User-provided input files (uploaded methylation CSVs) go to `data/input/`. This directory is gitignored. The Reflex UI writes uploaded files here, never inside `biomarkers-ui/` or any source directory.

### Output data (`data/output/`)

CLI commands that produce data default to writing under `data/output/`.

### Cache directory (cross-platform via `platformdirs`)

Long-lived cached data goes to the OS-appropriate user cache directory, resolved via `platformdirs.user_cache_dir("just-biomarkers")`. Override with `BIOMARKERS_CACHE_DIR` environment variable.

### Rules

- **NEVER commit large data files.** Parquet, CSV data files, and gzipped data must NEVER be added to git.
- **Library code** must use `platformdirs` or accept explicit paths. Never hardcode OS-specific cache paths.
- **UI uploaded files** must go to `data/input/`, never inside `biomarkers-ui/` or any source directory.
- **Never add data directories** to git. The `.gitignore` blocks `data/`, `output/`, and `**/uploaded_files/`.

---

## uv Project Management

- **Dependency Management**: Use `uv sync` and `uv add`. NEVER use `uv pip install`.
- **Project Configuration**: Use `pyproject.toml` as the single source of truth for dependencies and project metadata.
- **Versioning**: Do not hardcode versions in `__init__.py`; use `importlib.metadata`.

---

## Coding Standards

- **Type Hints**: Mandatory for all Python code.
- **Pathlib**: Always use `pathlib.Path` for file path operations.
- **Imports**: Always use absolute imports. No relative imports.
- **Error Handling**: Avoid nested `try-catch` blocks.
- **CLI Tools**: Use the `Typer` library for all CLI tools.
- **Data Classes**: Use `Pydantic 2` for all data validation and settings.
- **No Placeholders**: Never use temporary or custom local paths in committed code.
- **Prefer Polars over Pandas**: Use Polars for all data manipulation in the core library.

---

## Data Engineering (Polars)

- **Prefer Polars over Pandas**: Use `Polars` for data manipulation.
- **Efficiency**: Use `LazyFrame` (`scan_parquet`) and streaming (`sink_parquet`) for large datasets.
- **Memory Optimization**: Pre-filter dataframes before performing joins.

---

## Dagster 1.12.x+ API Notes & Gotchas

Many older Dagster tutorials use deprecated APIs. Keep these rules in mind for modern Dagster versions:

**Context Access:** `get_dagster_context()` does NOT exist. You must pass `context: AssetExecutionContext` explicitly to your functions.

**Metadata Logging:** `context.log.info()` does NOT accept a `metadata` keyword argument. Use `context.add_output_metadata()` separately.

**Run Logs:** `EventRecordsFilter` does NOT have a `run_ids` parameter. Instead, use `instance.all_logs(run_id, of_type=...)`.

**Asset Materializations:** Use `EventLogEntry.asset_materialization` (which returns `Optional[AssetMaterialization]`), not `DagsterEvent.asset_materialization`.

**Job Hooks:** The `hooks` parameter in `define_asset_job` must be a set, not a list (e.g., `hooks={my_hook}`).

**Asset Resolution:** Use `defs.resolve_all_asset_specs()` instead of the deprecated `defs.get_all_asset_specs()`.

**Asset Job Config:** Asset job config uses the `"ops"` key, not `"assets"`. Using `"assets"` causes a `DagsterInvalidConfigError`.

**CLI deprecation:** `dagster dev` is superseded by `dg dev`. The old command still works but emits a `SupersessionWarning`.

**Deprecation policy (CRITICAL):** treat all deprecation warnings as blockers for new changes in touched code. Investigate current upstream docs/APIs, update the implementation to non-deprecated APIs, and update `AGENTS.md` rules/examples so future changes do not reintroduce deprecated patterns.

**Definitions jobs deprecation warning fix:** do not pass unresolved asset jobs directly in `Definitions(jobs=[...])` when warning suggests resolution; resolve jobs by creating a temporary `Definitions` inside a function (local variable) to get the asset graph, then build the final `Definitions` with the resolved jobs. **CRITICAL: Dagster 1.12+ rejects multiple `Definitions` objects at module scope** — never assign a temporary `Definitions` to a module-level variable (even prefixed with `_`).

**Automation (CRITICAL — DO NOT USE `AutomationCondition`):** `AutomationCondition.on_missing()` and `.eager()` are **broken** for triggering initial materializations in Dagster 1.12. `on_missing()` on root assets silently produces 0 runs on every tick due to `InitialEvaluationCondition` canceling `SinceCondition`. `eager()` on root assets never fires (no upstream updates). `AutomationConditionSensorDefinition` starts `STOPPED` by default, and even when forced to `RUNNING`, the underlying conditions still produce 0 runs. The `dagster.yaml` `auto_materialize: enabled: true` is the legacy daemon and has no effect on the sensor system. **For startup, use a run-once bootstrap sensor** (`@dg.sensor` with `default_status=RUNNING`) that checks `instance.get_latest_materialization_event()` and submits a `RunRequest` with a `run_key` for deduplication. **For ongoing correctness, add a separate recompute sensor that triggers when upstream assets are newer than downstream outputs.**

### Resource Tracking (MANDATORY)

**Always track CPU and RAM consumption** for all compute-heavy assets using `resource_tracker` from `nanopore_pipeline.runtime`:

```python
from nanopore_pipeline.runtime import resource_tracker

@asset
def my_asset(context: AssetExecutionContext) -> Output[Path]:
    with resource_tracker("my_asset", context=context):
        # ... compute-heavy code ...
        pass
```

**Important:** Always pass `context=context` to enable Dagster UI metadata. Without it, metrics only go to the Dagster logger.
This automatically logs to Dagster UI: `duration_sec`, `cpu_percent`, `peak_memory_mb`, `memory_delta_mb`.

### Run-Level Resource Summaries (MANDATORY)

All jobs must include the `resource_summary_hook` from `nanopore_pipeline.utils` to provide aggregated resource metrics at the run level:

```python
from nanopore_pipeline.utils import resource_summary_hook

my_job = define_asset_job(
    name="my_job",
    selection=AssetSelection.assets(...),
    hooks={resource_summary_hook},  # Note: must be a set, not a list
)
```

This hook logs a summary at the end of each successful run: Total Duration, Max Peak Memory, and Top memory consumers.

### Key files for resource tracking

| File | What it does |
|------|-------------|
| `nanopore_pipeline/runtime.py` | `ResourceReport` model, `resource_tracker` context manager (uses `psutil`) |
| `nanopore_pipeline/utils.py` | `resource_summary_hook` — aggregates per-asset metrics into a run-level summary |

### Best Practices for Assets & IO

**Declarative Assets:** Prioritize Software-Defined Assets (SDA) over imperative ops. Include all assets in `Definitions(assets=[...])` for complete lineage visibility in the UI.

**Polars Integration:** Use `dagster-polars` with `PolarsParquetIOManager` for `pl.LazyFrame` assets to automatically get schema and row counts in the Dagster UI.

**Large Data / Streaming:** Use `lazy_frame.sink_parquet()` and NEVER `.collect().write_parquet()` on large data to avoid out-of-memory errors.

**Path Assets:** When returning a `Path` from an asset, add `"dagster/column_schema": polars_schema_to_table_schema(path)` to ensure schema visibility in the UI.

**Asset Checks:** Use `@asset_check` for validation and include them in your job via `AssetSelection.checks_for_assets(...)`.

### Execution & Concurrency Patterns

**Concurrency Limits:** Use `op_tags={"dagster/concurrency_key": "name"}` to limit parallel execution for resource-intensive assets.

**Timestamps:** Timestamps are on `RunRecord`, not `DagsterRun`. `run.start_time` will raise an `AttributeError`. Retrieve `instance.get_run_records()` and use `record.start_time`/`record.end_time` (Unix floats) or `record.create_timestamp` (datetime).

**Partition Keys for Runs:** `create_run_for_job` doesn't accept a direct `partition_key` parameter. Pass it via tags instead: `tags={"dagster/partition": partition_key}`.

**Dynamic Partitions Pattern:**
- Create partition def: `PARTS = DynamicPartitionsDefinition(name="files")`
- Discovery asset registers partitions: `context.instance.add_dynamic_partitions(PARTS.name, keys)`
- Partitioned assets use: `partitions_def=PARTS` and access `context.partition_key`
- Collector depends on partitioned output via `deps=[partitioned_asset]` and scans the filesystem/storage for results.

### Web UI / Asynchronous Execution Pattern

If you are running Dagster alongside a Web UI (like Reflex, FastAPI, etc.), use the Try-Daemon-With-Fallback pattern:

**Submission vs Execution:**
Attempt to submit the run to the daemon first: `instance.submit_run(run_id, workspace=None)`. If this fails (e.g., due to missing `ExternalPipelineOrigin` in web contexts), fall back to `job.execute_in_process()`.

**Rust/PyO3 Thread Safety:**
NEVER use `asyncio.to_thread()` or `asyncio.create_task()` with Dagster objects (it causes PyO3 panics: "Cannot drop pointer into Python heap without the thread being attached"). Use `loop.run_in_executor(None, sync_execution_function, ...)` for thread-safe background execution that doesn't block your UI.

**Orphaned Run Cleanup:**
If you use `execute_in_process` inside a web server process, runs will be abandoned (stuck in `STARTED` status) if the server restarts. Add startup cleanup logic targeting `DagsterRunStatus.NOT_STARTED`. Use `atexit` or signal handlers (`SIGTERM`/`SIGINT`) to mark active in-process runs as `CANCELED` on graceful server shutdown.

### Common Anti-Patterns to Avoid

- Using `dagster job execute` CLI: This is deprecated.
- Hardcoding asset names: Resolve them dynamically using defs.
- Suspended jobs holding locks: If a job crashes while querying local DBs (like DuckDB/SQLite), it can hold file locks. Handle connections properly via context managers or resources.
- Processing failures in exception handlers: Keep business logic out of exception handlers when executing runs. Catch the exception, register the failure, and cleanly proceed to your fallback mechanism.
- **Compute-heavy assets without `resource_tracker`** — if a process gets OOM-killed, there are no metrics to diagnose it. Always wrap with `resource_tracker(name, context=context)`.
- **Jobs without `resource_summary_hook`** — without it, run-level resource consumption is invisible. Always pass `hooks={resource_summary_hook}` to `define_asset_job`.
- **NEVER use `AutomationCondition` for hands-free pipeline launch.** `on_missing()` silently produces 0 runs on root assets, `eager()` never fires on root assets, and `AutomationConditionSensorDefinition` doesn't fix either issue. Use a run-once `@dg.sensor` instead.
- **NEVER use `auto_materialize: enabled: true` in `dagster.yaml`.** This is the legacy daemon and has no effect on `AutomationCondition` or sensors.
- **NEVER use `AutomationConditionSensorDefinition`** as a fix for the above — the underlying `AutomationCondition` logic is what's broken, not the sensor wrapper.
- **Default to caching; force re-materialization only via explicit `--no-cache`.** If cached data exists, use it.
- **NEVER use `os.execvp` for headless pipeline execution.** Use `job.execute_in_process()` for `--headless` runs. `os.execvp` is only for the default UI mode (used by all commands without `--headless`).

---

## Dagster Pipeline Robustness Policy (CRITICAL)

Standard bioinformatics pipeline engines (Nextflow, WDL/Cromwell, Snakemake) provide out-of-the-box guarantees: timestamp-based cache invalidation, automatic failure retry, interrupted-run resume, and completeness validation. Dagster does NOT provide these by default — sensors only check metadata events (a database flag), not actual data state. Every Dagster asset, sensor, and job in this project MUST implement the 8 guarantees below to match or exceed these standards.

### The 8 Robustness Guarantees

**1. Input-change invalidation (mtime-based).**
If an upstream input file changes (new BED file downloaded, source data refreshed), downstream cached results that depended on the old input MUST be recomputed. Detection mechanism: compare input file mtime against output file mtime. If the input is newer than the output, the cache is stale — recompute. This matches `make`/Nextflow/Snakemake timestamp-based rebuild. In `skip_existing` checks, file existence alone is NOT sufficient — always verify that the output is at least as recent as its input.

**2. Interrupted-run recovery (gap detection).**
If a process crashes, is OOM-killed, or is interrupted mid-batch, the next run MUST detect incomplete work and resume from where it left off. Detection mechanism: compare the set of expected outputs against the set of outputs that actually exist on disk. Any gap triggers recomputation of the missing items only — not a full re-run.

**3. Automatic failure retry.**
If individual items in a batch fail (e.g. a BED file download throws an exception), those failures MUST be automatically retried on subsequent runs. After N consecutive retries with the same failure set (configurable, default 3), stop retrying and log a permanent-failure warning — do not retry forever.

**4. Completeness validation before publishing.**
Before uploading results to an external destination (HuggingFace, S3, etc.), the asset MUST validate that the output has acceptable coverage. Log `coverage_ratio`, `n_ok`, `n_total`, `n_failed`, `n_missing` in output metadata. Warn (but still upload) if coverage is below a configurable threshold. NEVER silently upload a near-empty dataset.

**5. Content-aware cache validity (skip_existing).**
`skip_existing` checks MUST NOT only test for file existence. They must also verify that the cached output is not stale relative to its inputs. Minimum check: compare output file mtime against input file mtime. If input is newer, the cache is invalid and the item must be recomputed. This is the Nextflow `-resume` / Snakemake / Make default behavior.

**6. Rich batch metadata.**
Every asset that processes a batch of items MUST emit in its Dagster output metadata: `n_total` (items attempted), `n_ok` (succeeded), `n_failed` (failed), `n_cached` (reused from cache), `coverage_ratio` (`n_ok / n_total`). This makes the Dagster UI a self-documenting data catalog where you can see at a glance whether a run was complete or partial.

**7. Sensor intelligence hierarchy.**
Sensors MUST check data state in this priority order:
  1. Is a run already active for this job? Skip — never double-submit.
  2. Was this an explicit user command (force run)? Submit unconditionally.
  3. Are there failed items from the last run that should be retried? Submit targeted retry.
  4. Are there missing items (expected vs on-disk)? Submit full job with `skip_existing`.
  5. Have upstream inputs changed (fingerprint/mtime changed)? Submit full job.
  6. Are all items present and up-to-date? Skip — everything is healthy.

**8. Corrupted file detection and recovery.**
A file that exists on disk but cannot be parsed (truncated writes, OOM-killed mid-flush, invalid thrift headers, etc.) is **worse than a missing file** — existence checks pass, recomputation is skipped, but reading the file later raises an exception that crashes the entire run. Every place that reads a cached parquet MUST catch parse errors (`polars.exceptions.ComputeError` and similar) and treat the file as missing: delete the corrupt file and fall through to recompute.

### Robustness Anti-Patterns (NEVER)

- **NEVER** use Dagster materialization events as the sole indicator of completeness. Materialization means "the asset function returned" — it says nothing about how many items succeeded, failed, or were skipped.
- **NEVER** upload to an external destination without logging coverage metrics. Silent uploads of near-empty datasets are worse than no upload.
- **NEVER** skip a cached result without checking that its input is older than the output. Existence-only checks lead to stale data that is never refreshed.
- **NEVER** ignore failures from a previous run. If a quality report exists with `status == "failed"` entries, the next run must attempt to resolve them.
- **NEVER** rely on a single sensor for all orchestration concerns. Separate sensors for: startup/force-run, failure-retry, completeness-gap, and upstream-freshness. Each has different check intervals and submission logic.
- **NEVER** treat file existence as proof of a valid cache. A truncated write or OOM-killed process leaves a corrupt parquet that passes the existence check but crashes the run when read. Always validate readability (`scan_parquet` + `collect_schema()`) before skipping recomputation. Delete and recompute any file that fails to parse.

### Quick Checklist — Before Finishing Any Dagster Asset

- [ ] `skip_existing` checks compare mtimes, not just existence.
- [ ] `skip_existing` cache readers validate parquet readability (`collect_schema()`), delete corrupt files, and fall through to recompute.
- [ ] Output metadata includes `n_total`, `n_ok`, `n_failed`, `n_cached`, `coverage_ratio`.
- [ ] Upload assets validate coverage and log `n_scored`, `n_total`, `n_missing`.
- [ ] Failures are recorded in a quality parquet that sensors can read for retry.
- [ ] The asset handles partial prior results (resumes from where it left off).

---

## Dagster Asset Lineage Rules (CRITICAL)

Incomplete lineage is a recurring mistake. Every time you write a Dagster asset, apply ALL of the rules below before finishing.

### 1. SourceAssets for every external data source

Any data that Dagster **downloads but does not produce** — an FTP tarball, a remote API, a Synology/GoFile directory, a HuggingFace dataset — must be declared as a `SourceAsset`. This gives it a visible node in the lineage graph.

```python
from dagster import SourceAsset

synology_nanopore_source = SourceAsset(
    key="synology_nanopore_source",
    group_name="external",
    description="Synology/GoFile directory hosting nanopore methylation BED files.",
    metadata={"url": "http://gofile.me/76Mv5/irFIggw48"},
)
```

Register every `SourceAsset` in `Definitions(assets=[...])`. Without this, the left edge of the graph is a dangling node with no visible origin.

### 1b. NEVER put SourceAsset keys in computed asset `deps`

SourceAssets are visualization-only nodes — Dagster never materializes them. If a computed asset lists a SourceAsset in `deps`, Dagster will see the SourceAsset as permanently missing, which can block job execution and sensor-based automation.

```python
# WRONG — synology_nanopore_source is a SourceAsset → deps check will fail
@asset(deps=["synology_nanopore_source"])
def synology_file_index(...): ...

# CORRECT — no SourceAsset in deps; document the source URL in output metadata instead
@asset()
def synology_file_index(context, ...):
    ...
    context.add_output_metadata({"source_url": MetadataValue.url(url)})
```

### 2. Always declare `deps` when an asset reads filesystem side effects of another asset

If an asset scans a directory that was populated by another asset (common pattern: partitioned writer -> non-partitioned aggregator), it MUST declare a `deps` relationship. Without `deps`, Dagster draws NO edge between them even though there is a real data dependency.

Rule: if you write "scan for files produced by X" anywhere in a docstring or comment, you MUST add `deps=[AssetDep("X")]` to the decorator.

### 3. Use `AssetIn` only when data is passed through the IOManager

`AssetIn` means "load the output of asset X via the IOManager and inject it as a function argument". Use it when:
- The upstream returns a `pl.DataFrame` / `pl.LazyFrame` and you use `PolarsParquetIOManager`
- The upstream returns a picklable Python object and the default IOManager is acceptable

Do NOT use `AssetIn` with `Path` or `Output[Path]` unless you have a custom `UPathIOManager`. For `Path`-based data flow, use `deps` for lineage and reconstruct the path inside the downstream asset from a shared resource (e.g. `CacheDirResource`).

### 4. Separate assets into named groups by pipeline stage

Every `@asset` and `SourceAsset` must have a `group_name`. Use a consistent stage-based taxonomy so the Dagster UI shows a clear left-to-right graph:

| group_name | What belongs here |
|---|---|
| `external` | `SourceAsset`s for remote data (FTP, Synology, HF, S3, API) |
| `download` | Assets that fetch external data into local cache |
| `compute` | Transformation / parsing / aggregation assets |
| `upload` | Assets that push results to external destinations (HF, S3, DB) |

### 5. SourceAssets must be registered in `Definitions`

Omitting a `SourceAsset` from `Definitions` makes it invisible in the UI even if assets declare `deps` on it.

### 6. Universal Dagster Principles (Mindset & Architecture)

- **Assets over Tasks (Declarative Mindset)**: Focus on *what data should exist*, not *how to run a function*. Dependencies are expressed as data assets, not task wiring. Lineage is automatic based on these data dependencies.
- **Dynamic Partitions for Entities**: When processing many independent entities (like samples), use Dynamic Partitions for targeted materialization, backfills, and deletions. However, when each entity processes in seconds and shares expensive setup, a batch approach in a single asset with in-process iteration and error tracking is more efficient than thousands of partitions + sensor orchestration.
- **Assets vs. Jobs Separation**: Use Software-Defined Assets (`@asset`) for the declarative data graph and lineage. Use Jobs (`define_asset_job` / `ops`) strictly as operational entry points (for sensors, schedules, CLI, or UI triggers).
- **Abstracted Storage (IO Managers & Resources)**: Avoid hardcoding paths inside asset logic. Either return a value and let an IO Manager write it, or reconstruct the path from a shared resource (like `CacheDirResource`). This prevents path conflicts and separates business logic from storage concerns.
- **Rich Metadata**: Always add meaningful output metadata (`context.add_output_metadata(...)`), such as row counts, file sizes, schema details, or external URLs. This turns Dagster into an inspectable data catalog rather than just a task runner.
- **Freshness Over Presence**: "Asset exists" is not sufficient. Sensors/schedules must compare upstream vs downstream materialization recency and trigger recompute when lineage indicates stale outputs.

### Quick checklist before finishing any Dagster asset file

- [ ] Every remote data origin has a `SourceAsset` with `group_name="external"` and a `metadata={"url": ...}`.
- [ ] Every `SourceAsset` is listed in `Definitions(assets=[...])`.
- [ ] Every asset that scans a directory written by another asset has `deps=[AssetDep("that_asset")]`.
- [ ] No `AssetIn` is used with `Output[Path]` unless a `UPathIOManager` is configured.
- [ ] Every `@asset` has `group_name` set (never omit it).
- [ ] The sensor uses `run_key` to prevent duplicate submissions across ticks (fresh key only on retry after failure).
- [ ] Every asset checks on-disk cache and short-circuits if data already exists — no redundant downloads or computations.
- [ ] Every compute-heavy asset is wrapped with `resource_tracker(name, context=context)`.
- [ ] Every job has `hooks={resource_summary_hook}` for run-level resource aggregation.

---

## Dagster Pipeline Execution Model (CRITICAL)

### Core principle: never force re-materialization

Every pipeline command respects caches. If data is already on disk, it is reused. No CLI command or sensor should ever force a full re-run when cached data exists. The only triggers for re-computation are: (1) assets have never been materialized, or (2) an upstream asset has been freshly materialized (making downstream stale).

### CLI commands

All pipeline commands launch the Dagster UI by default. Headless mode requires explicit `--headless`.

| Command | What it does | When it re-materializes |
|---------|-------------|------------------------|
| `nanopore run` | Launches Dagster UI; startup sensor submits `full_pipeline` if assets are missing. | Each asset checks disk cache and short-circuits if data exists. |
| `nanopore run --headless` | Executes `full_pipeline` in-process (no UI). | Same cache-respecting behavior, but no UI monitoring. |
| `nanopore launch` | Launches Dagster UI (no specific job pre-selected). | The sensor submits a job only if key assets are unmaterialized. |

### Why `os.execvp` (not `subprocess.run` or `Popen`)

Dagster's daemon uses complex internal signal handling. When trapped inside a `subprocess.run()` or `Popen()`, SIGINT/SIGTERM do not propagate correctly and the daemon does not shut down cleanly. `os.execvp` **replaces** the current Python process with `dagster`, so Dagster becomes the primary process and owns all signal handling. Ctrl+C works correctly.

### Why startup sensor (DO NOT use `AutomationCondition`)

Because `os.execvp` replaces the current process, there is no opportunity to submit a job "after dagster starts". **Do NOT use `AutomationCondition`** (`on_missing()`, `eager()`, or `AutomationConditionSensorDefinition`) for hands-free pipeline launch — they are broken for initial materialization in Dagster 1.12.

**Reliable pattern is a run-once bootstrap sensor** that checks materialization status and submits a job only when assets are missing:

```python
import dagster as dg

@dg.sensor(
    job_name="full_pipeline",
    default_status=dg.DefaultSensorStatus.RUNNING,
    minimum_interval_seconds=60,
)
def startup_sensor(context: dg.SensorEvaluationContext) -> dg.SensorResult | dg.SkipReason:
    check_keys = [
        dg.AssetKey("synology_file_index"),
        dg.AssetKey("nanopore_bed_files"),
        dg.AssetKey("nanopore_methylation_matrix"),
    ]
    missing = [k for k in check_keys if context.instance.get_latest_materialization_event(k) is None]
    if not missing:
        return dg.SkipReason("All pipeline assets already materialized.")

    active = context.instance.get_runs(
        filters=dg.RunsFilter(job_name="full_pipeline", statuses=[
            dg.DagsterRunStatus.STARTED, dg.DagsterRunStatus.NOT_STARTED, dg.DagsterRunStatus.QUEUED,
        ])
    )
    if active:
        return dg.SkipReason(f"Already in progress (run {active[0].run_id[:8]}).")

    last_runs = context.instance.get_runs(filters=dg.RunsFilter(job_name="full_pipeline"), limit=1)
    if last_runs and last_runs[0].status == dg.DagsterRunStatus.FAILURE:
        run_key = f"nanopore_startup_retry_{int(time.time())}"
    else:
        run_key = "nanopore_startup"

    return dg.SensorResult(
        run_requests=[dg.RunRequest(run_key=run_key, job_name="full_pipeline")],
    )
```

The `run_key` prevents duplicate submissions. The sensor only generates a fresh key after a failure, allowing retries.

### Anti-patterns to NEVER use

- **Timestamp-based run keys** (e.g. `f"startup_{int(time.time())}"`) — defeats Dagster's deduplication, forces re-runs every time.
- **Implicit force-run on startup** — the sensor should only submit when assets are missing.
- **`os.execvp` for headless execution** — use `execute_in_process()` only for `--headless` runs. `os.execvp` is for the default UI mode (all commands without `--headless`).

### `dagster.yaml` template

```yaml
telemetry:
  enabled: false
```

Generate this file at `{DAGSTER_HOME}/dagster.yaml` if it does not exist. The `telemetry: enabled: false` setting prevents `RotatingFileHandler` crashes in Dagster's event log writer. Do NOT add `auto_materialize: enabled: true` — that is the legacy daemon approach and does not work with the sensor pattern above.

### Port cleanup helper

```python
def _kill_port(port: int) -> None:
    import subprocess, signal, time, os
    result = subprocess.run(["lsof", "-t", f"-iTCP:{port}"], capture_output=True, text=True)
    pids = [int(p) for p in result.stdout.strip().splitlines() if p.strip()]
    for pid in pids:
        os.kill(pid, signal.SIGTERM)
    if pids:
        time.sleep(1)
        result2 = subprocess.run(["lsof", "-t", f"-iTCP:{port}"], capture_output=True, text=True)
        for pid in [int(p) for p in result2.stdout.strip().splitlines() if p.strip()]:
            os.kill(pid, signal.SIGKILL)
```

### Orphaned run cleanup helper

```python
def _cancel_orphaned_runs() -> None:
    from dagster import DagsterInstance, DagsterRunStatus, RunsFilter
    with DagsterInstance.get() as instance:
        stuck = instance.get_run_records(
            filters=RunsFilter(statuses=[DagsterRunStatus.STARTED, DagsterRunStatus.NOT_STARTED])
        )
        for record in stuck:
            instance.report_run_canceled(record.dagster_run, message="Orphaned run from previous session")
```

### Thread-safety note for Web UI contexts

If a Web UI (Reflex, FastAPI) needs to trigger a Dagster job in the background, **never use `asyncio.to_thread()`**. Dagster's Rust/PyO3 internals panic when the GIL is released across asyncio threads. Use `loop.run_in_executor(None, sync_func)` instead.

---

## Learned Dagster Workspace Facts

- Dagster 1.12+ rejects multiple `Definitions` objects at module scope — never assign a temporary `Definitions` to a module-level variable (even prefixed with `_`); wrap resolution in a function
- Dagster's SQLite storage produces `database is locked` errors under concurrent writes (multiple runs + daemon + sensors) — these are transient and do not stop runs, but PostgreSQL (`dagster-postgres`) is the production alternative
- Dagster's `AutomationCondition.on_missing()` and `.eager()` are broken for initial materialization in Dagster 1.12 — use run-once `@dg.sensor` instead
- Dagster UI mode uses `os.execvp` because daemon signal handling breaks under `subprocess.run()`/`Popen()`
- Dagster `dagster dev` defaults to the **multiprocess executor**, which forks child processes. Arrow/numpy C extensions are not fork-safe — SIGABRT (signal 6) crashes result. All pipeline jobs must use `executor_def=in_process_executor` to avoid forking.
- Dagster asset materialization (green checkmark) is separate from asset checks — checks can fail while assets stay green
- `pyo3_runtime.PanicException` inherits from `BaseException`, not `Exception` — the standard `except Exception` in batch loops will NOT catch polars/pyo3 Rust panics. Use `except (KeyboardInterrupt, SystemExit): raise` then `except BaseException` to catch panics without swallowing interrupts.
- `polars collect_schema()` reads only parquet footer metadata — data-page corruption (truncated writes, OOM-killed flushes) is only detected when `collect()` reads row groups; always wrap the full read sequence in try/except for corrupt parquet handling

---

## Test Generation Guidelines

- **Real Data + Ground Truth**: Use actual source data and compute expected values at runtime.
- **biolearn parity**: Integration tests compare our scores against biolearn for ~10 linear clocks on the same input data with tolerance-based assertions (1e-4 to 1e-3).
- **No Mocking**: Do not mock data transformations. Run real compute pipelines.
- **Meaningful Assertions**: Prefer relationship and aggregate checks over simple existence checks.

---

## Learned User Preferences

- Always avoid adding biolearn as a runtime dependency; re-implement any needed functionality using stdlib + polars + platformdirs instead.
- When adding CLI commands for example dataset downloads, use `data/input/examples/` as the default output directory convention.
- Parity tests against biolearn must use real data and real pipelines — do not mock data transformations.

## Learned Workspace Facts

- `biolearn.quantile_normalize_using_target` has a read-only numpy array bug in newer pandas/numpy versions; parity tests for DunedinPACE (and any clock using quantile normalization) must monkeypatch this function in `conftest.py` to make the input array writable before calling biolearn as the oracle.
- `download-epic` downloads all EPIC-format GEO datasets to `data/input/examples/`; `download-all-examples` downloads all datasets (EPIC + 450k); both are subcommands of the `biomarkers` CLI.
- `uv run biomarkers compare-epic` is the CLI command for comparing just-biomarkers scores vs biolearn scores on downloaded EPIC datasets side-by-side (logic lives in `just_biomarkers.compare`).
- ComputAgeBench HF parquet files (`computage/computage_bench`, paths `data/benchmark/` and `data/train/`) are already oriented rows=CpGs, cols=samples, matching our standard; the first column is the CpG site ID and must be renamed to `CpGmarker` on load.
- The full ComputAgeBench HF snapshot is many GB; `huggingface_hub.snapshot_download` is used for caching; individual dataset parquet files are loaded on demand via `load_computage_dataset()`.
