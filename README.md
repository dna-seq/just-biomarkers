# just-biomarkers workspace

A [uv workspace](https://docs.astral.sh/uv/concepts/workspaces/) for methylation biomarker computation — from raw data (Illumina arrays, nanopore BED files) to epigenetic clock scores.

## Subprojects

| Package | Directory | Status | Description |
|---------|-----------|--------|-------------|
| **[just-biomarkers](just-biomarkers/)** | `just-biomarkers/` | **Core library — production-ready** | Polars-based methylation clock scoring engine. 44 clocks, 11 parity-tested against biolearn. CSV/TSV/Parquet I/O, GEO downloader, CLI. |
| **[biomarkers-ui](biomarkers-ui/)** | `biomarkers-ui/` | **Functional, no tests** | Reflex web app for interactive clock computation. Upload methylation data, select clocks, view results. |
| **[nanopore-pipeline](nanopore-pipeline/)** | `nanopore-pipeline/` | **Scaffolded, never run** | Dagster pipeline for downloading nanopore BED files and building a methylation matrix. Code exists but has never been executed. |

## Quick Start

```bash
# Install all workspace packages
uv sync --all-packages

# List available methylation clocks
uv run biomarkers list-clocks

# Compute clock scores from a methylation CSV
uv run biomarkers compute methylation_matrix.csv --clocks Horvathv1,Hannum,PhenoAge

# Download an example GEO dataset
uv run biomarkers download-example GSE112618

# Launch the web UI
uv run ui
```

## Python API

```python
from just_biomarkers import (
    score_clock, score_clocks,
    list_clocks, get_clock, search_clocks,
    read_methylation_csv,
)

dnam = read_methylation_csv("methylation_matrix.csv")
result = score_clocks(dnam, clock_names=["Horvathv1", "Hannum", "PhenoAge"])

for r in result.results:
    print(f"{r.sample_id}: {r.clock_name} = {r.score:.4f} ({r.match_rate:.1%})")
```

## Readiness Overview

### just-biomarkers (core library)

The scoring engine is the most mature part of the project. It supports 44 methylation clocks with a Polars-native compute pipeline (load coefficients, join, weighted sum, transform). Key readiness details:

| Component | Status | Notes |
|-----------|--------|-------|
| Scoring engine | **Production-ready** | Parity-tested against biolearn on 10 clocks across 3 real GEO datasets (6, 49, 143 samples) |
| DunedinPACE preprocessing | **Production-ready** | Quantile normalization reimplemented in pure NumPy, verified against biolearn to < 1e-3 |
| I/O (CSV/TSV/Parquet) | **Production-ready** | Unit-tested with column renaming, type casting, auto-detection |
| I/O (IDAT) | **Early stage** | Code exists using `pylluminator-modern` but no end-to-end test on real IDAT data |
| Imputation (3 strategies) | **Early stage** | Code written, zero tests — all parity tests bypass imputation |
| GEO downloader | **Functional** | 5 curated datasets; download logic works (parity tests depend on it) |
| ComputAgeBench downloader | **Early stage** | Download/load functions exist, no tests — correctness unverified |
| `compare-epic` CLI | **Functional** | Interactive comparison tool, works when biolearn dev dependency is available |

11 of 44 clocks are parity-verified against biolearn (Horvathv1, Hannum, PhenoAge, Lin, Horvathv2, DunedinPoAm38, DunedinPACE, Zhang_10, YingCausAge, YingDamAge, YingAdaptAge). The remaining 33 use the same scoring engine and are likely correct but unverified.

See [just-biomarkers/README.md](just-biomarkers/README.md) for the full clock catalogue, architecture details, and test matrix.

### biomarkers-ui (web app)

Reflex web app with Plotly charts. Upload a methylation CSV or load a GEO example dataset, select clocks from a checkbox list, compute scores, and view/download results. Depends on `just-biomarkers` for all computation.

**No tests** — this is expected for a UI layer. The app is functional but has not been extensively user-tested.

```bash
uv sync --all-packages
uv run ui
```

### nanopore-pipeline (Dagster)

Dagster pipeline with 3 assets (`synology_file_index`, `nanopore_bed_files`, `nanopore_methylation_matrix`), a startup sensor, CLI with UI/headless modes, and resource tracking. The code follows the project's Dagster conventions (in-process executor, resource tracking, rich metadata, orphaned run cleanup).

**This pipeline has never been launched or executed against real data.** It is purely scaffolded code:
- The GoFile/Synology scraper may not work (page may require authentication or have a different HTML structure)
- The BED parser assumes nanopore modkit output format but has never parsed a real file
- No tests of any kind exist

```bash
uv run nanopore run             # Dagster UI (never tested)
uv run nanopore run --headless  # In-process (never tested)
uv run nanopore status          # Check cache state
```

## CLI Entry Points

All commands are available via `uv run <name>` from the workspace root:

| Command | Entry Point | Subproject |
|---------|-------------|------------|
| `biomarkers` / `just-biomarkers` | `just_biomarkers.cli:app` | just-biomarkers |
| `ui` | `biomarkers_ui.cli:launch_ui` | biomarkers-ui |
| `nanopore` | `nanopore_pipeline.cli:app` | nanopore-pipeline |

### biomarkers CLI commands

```bash
uv run biomarkers list-clocks                    # List all 44 clocks
uv run biomarkers list-clocks --tissue Blood     # Filter by tissue
uv run biomarkers compute matrix.csv             # Score all clocks
uv run biomarkers compute matrix.csv --clocks Horvathv1,PhenoAge -o results.csv
uv run biomarkers download-example GSE112618     # Download a GEO dataset
uv run biomarkers download-example --list        # List available datasets
uv run biomarkers download-epic                  # Download all EPIC datasets
uv run biomarkers download-all-examples          # Download all datasets
uv run biomarkers download-computage GSE40279    # Download ComputAgeBench dataset
uv run biomarkers compare-epic                   # Compare vs biolearn (dev dep)
```

## Running Tests

```bash
# Unit tests (no external data needed)
uv run pytest just-biomarkers/tests/test_io.py -v

# Parity tests (require biolearn dev dep + downloaded GEO data)
uv run biomarkers download-epic
uv run pytest just-biomarkers/tests/ -v
```

## Development

```bash
# Sync workspace with dev dependencies
uv sync --all-packages --group dev

# Run all tests
uv run pytest just-biomarkers/tests/ -v

# Compare scores against biolearn interactively
uv run biomarkers compare-epic --dataset GSE112618 --clocks PhenoAge,Horvathv2
```

## Data Conventions

- **Input data**: `data/input/` (gitignored) — uploaded files, downloaded GEO examples
- **Output data**: `data/output/` (gitignored) — CLI output, Dagster home
- **Cache**: `~/.cache/just-biomarkers/` via `platformdirs` (override with `BIOMARKERS_CACHE_DIR`)
- **Methylation standard**: rows = CpG sites (`CpGmarker`), columns = samples, values = beta [0,1]

## License

MIT
