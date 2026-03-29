# nanopore-pipeline

Dagster pipeline for downloading and processing nanopore methylation data from Synology/GoFile sources.

## Pipeline Assets

| Asset | Group | Description |
|-------|-------|-------------|
| `synology_file_index` | `download` | Discover available `.bed.gz` files from the Synology GoFile directory listing |
| `nanopore_bed_files` | `download` | Download nanopore methylation BED files to local cache |
| `nanopore_methylation_matrix` | `compute` | Parse BED files and build a unified CpG methylation beta-value matrix |

## Usage

```bash
# From workspace root
uv sync --all-packages

# Launch Dagster UI (default)
uv run nanopore run

# Launch Dagster UI only (no job pre-selected)
uv run nanopore launch

# Run headless (no UI)
uv run nanopore run --headless

# Check download status
uv run nanopore status
```

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `NANOPORE_PIPELINE_HOST` | `0.0.0.0` | Dagster webserver bind address |
| `NANOPORE_PIPELINE_PORT` | `3020` | Dagster webserver port |
| `BIOMARKERS_CACHE_DIR` | OS-specific | Override cache directory |
| `NANOPORE_SOURCE_URL` | `http://gofile.me/76Mv5/irFIggw48` | Source URL for nanopore data |
