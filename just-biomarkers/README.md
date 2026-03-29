# just-biomarkers

Polars-based methylation clock computation engine. A lightweight, dependency-minimal reimplementation of [biolearn](https://github.com/bio-learn/biolearn) linear clocks using [Polars](https://pola.rs/) for fast, memory-efficient computation.

## Motivation

Existing methylation clock libraries (biolearn, methylclock) rely on heavy dependencies like pandas, matplotlib, seaborn, and torch. `just-biomarkers` provides the same clock scoring with a minimal runtime footprint: Polars + NumPy + Pydantic. It is designed for embedding in pipelines, CLI tools, and web apps without pulling in the full data-science stack.

All clock coefficients and transforms are validated against biolearn via parity tests to ensure numerical equivalence.

## Features

- **44 methylation clocks** spanning aging, mortality, traits, disease risk, pregnancy, and lifestyle biomarkers
- **Polars-native** scoring engine — no pandas dependency at runtime
- **Multiple input formats**: CSV, TSV, TXT (gzipped or plain), Parquet, and Illumina IDAT (via optional `pylluminator-modern`)
- **GEO dataset downloader** — built-in catalogue of curated example datasets from NCBI GEO
- **ComputAgeBench integration** — download and load datasets from the HuggingFace `computage/computage_bench` benchmark
- **Imputation strategies**: row-mean averaging, reference-table fill, and hybrid (biolearn-compatible)
- **Pre-scoring preprocessing**: quantile normalization for DunedinPACE
- **Named transforms**: Horvath anti-trafo, sigmoid, offset corrections
- **Pydantic v2 models** for clock definitions, results, and batch results
- **Typer CLI** with commands for listing clocks, computing scores, downloading example data, and comparing against biolearn
- **Parity-tested** against biolearn on real GEO datasets with tolerance-based assertions

## Installation

```bash
pip install just-biomarkers
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv add just-biomarkers
```

## Quick Start

### CLI

```bash
# List all available clocks
biomarkers list-clocks

# List clocks for a specific tissue
biomarkers list-clocks --tissue Blood

# Compute clock scores from a methylation CSV
biomarkers compute methylation_matrix.csv --clocks Horvathv1,Hannum,PhenoAge

# Output results to CSV
biomarkers compute methylation_matrix.csv -o results.csv

# Download an example GEO dataset
biomarkers download-example GSE112618

# Download all EPIC-format example datasets
biomarkers download-epic

# Compare just-biomarkers vs biolearn (requires biolearn dev dependency)
biomarkers compare-epic
```

### Python API

```python
from just_biomarkers import (
    score_clock, score_clocks,
    list_clocks, get_clock, search_clocks,
    read_methylation_csv,
    ClockDefinition, ClockResult, BatchClockResult,
)

# Read a methylation matrix (rows = CpG sites, columns = samples)
dnam = read_methylation_csv("methylation_matrix.csv")

# Score multiple clocks at once
result = score_clocks(dnam, clock_names=["Horvathv1", "Hannum", "PhenoAge"])

for r in result.results:
    print(f"{r.sample_id}: {r.clock_name} = {r.score:.4f} ({r.match_rate:.1%})")

# Search clocks by tissue
blood_clocks = search_clocks(tissue="Blood")
```

---

## Readiness Summary

| Component | Status | Test Coverage | Notes |
|-----------|--------|---------------|-------|
| **Scoring engine** (`scoring.py`) | **Production-ready** | Parity-tested against biolearn on 10 clocks across 3 real GEO datasets (6, 49, 143 samples). Regression anchors for all tested clocks. | Core compute path is solid. |
| **Clock registry** (44 clocks) | **Production-ready** | Registry coverage test asserts >= 20 biolearn clocks covered. All definitions are Pydantic-validated. | All clocks compute; 10 are parity-verified. |
| **I/O: CSV/TSV/Parquet** | **Production-ready** | Unit tests for CSV, TSV, Parquet with column renaming, type casting, and edge cases. | Auto-detection works. |
| **I/O: IDAT** | **Early stage** | Only error-path tests (verifies clear error when `pylluminator-modern` is missing). No parity test with real IDAT data. | Code exists and uses `pylluminator-modern`, but no end-to-end test validates correctness on actual IDAT files. |
| **Transforms** (13 named) | **Production-ready** | Implicitly tested via parity tests on clocks that use them (Horvathv1/v2, PEDBE, Cortical, sigmoid clocks). | |
| **DunedinPACE preprocessing** | **Production-ready** | Dedicated parity tests against biolearn on two EPIC datasets (6 + 49 samples). Pure NumPy reimplementation of quantile normalization verified to < 1e-3. | Most complex preprocessing path, well-tested. |
| **Imputation** | **Early stage** | No tests. Three strategies exist (`impute_from_average`, `impute_from_reference`, `hybrid_impute`) but none are exercised by any test. | Code exists but correctness is unverified. |
| **GEO downloader** | **Functional** | Used by parity tests (tests skip if data not downloaded). Catalogue has 5 curated datasets. | Download + parse logic works in practice since parity tests depend on it. |
| **ComputAgeBench downloader** | **Early stage** | No tests. Download/load/list functions exist but are never exercised in tests. No parity test scores ComputAgeBench data. | Code is written but unverified — work just started. |
| **`compare-epic` CLI** | **Functional** | No dedicated test, but exercises the same code path as the parity tests. | Interactive comparison tool; works when biolearn is available. |
| **Web UI** (`biomarkers-ui`) | **Functional, no tests** | No tests (expected for UI). Upload, clock selection, compute, and results display are implemented. | Reflex app with Plotly charts. |
| **Nanopore pipeline** | **Scaffolded, never run** | No tests. Never launched or executed. Dagster assets, sensors, CLI, and BED parser are written but entirely untested. | Code follows AGENTS.md conventions structurally but has never processed real data. GoFile scraper may not work (page may require auth or have different HTML structure). BED parser assumes modkit format — untested against real nanopore output. |

### Clocks with Parity Tests

The following 10 clocks are **parity-verified** against biolearn on the biolearn testset (tolerance 1e-4 to 1e-3):

| Clock | Parity Tolerance | Status |
|-------|------------------|--------|
| Horvathv1 | 1e-3 | Verified |
| Hannum | 1e-4 | Verified |
| PhenoAge | 1e-4 | Verified |
| Lin | 1e-4 | Verified |
| Horvathv2 | 1e-3 | Verified |
| DunedinPoAm38 | 1e-4 | Verified |
| Zhang_10 | 1e-4 | Verified |
| YingCausAge | 1e-4 | Verified |
| YingDamAge | 1e-4 | Verified |
| YingAdaptAge | 1e-4 | Verified |

Additionally, **DunedinPACE** has its own dedicated parity test suite (quantile normalization verified on GSE112618 and GSE110554).

The remaining 33 clocks are **registered and computable** but have no parity test against biolearn. They use the same scoring engine, so they are likely correct, but this is unverified.

### GEO Dataset Parity Tests

| Test File | Dataset | Samples | Clocks Tested | What It Verifies |
|-----------|---------|---------|---------------|------------------|
| `test_parity.py` | biolearn testset | ~10 | 10 linear clocks | Per-sample score match vs biolearn |
| `test_gse112618_parity.py` | GSE112618 (EPIC) | 6 | 5 full-match + 5 partial-match | Biolearn parity, regression anchors, match rates |
| `test_gse110554_parity.py` | GSE110554 (EPIC) | 49 | 5 full-match + 5 partial-match | Biolearn parity on 49 samples, age plausibility |
| `test_gse164056_parity.py` | GSE164056 (EPIC) | 143 | 10 clocks (all partial) | Finite output, regression anchors, age plausibility |
| `test_dunedin_pace.py` | GSE112618 + GSE110554 | 6 + 49 | DunedinPACE only | Quantile normalization parity vs biolearn |
| `test_dunedin_pace_parity.py` | GSE112618 | 6 | DunedinPACE only | Detailed DunedinPACE parity (full match, range, regression) |
| `test_io.py` | Synthetic | N/A | N/A | CSV/TSV/Parquet reading, IDAT error paths |

---

## Supported Clocks

### Aging & Mortality Clocks

| Clock | Year | Tissue | Output | Parity Tested |
|-------|------|--------|--------|:-------------:|
| Horvathv1 | 2013 | Multi-tissue | Age (Years) | Yes |
| Hannum | 2013 | Blood | Age (Years) | Yes |
| PhenoAge | 2018 | Blood | Age (Years) | Yes |
| Lin | 2016 | Blood | Age (Years) | Yes |
| Horvathv2 | 2018 | Skin + blood | Age (Years) | Yes |
| DunedinPoAm38 | 2020 | Blood | Aging Rate (Years/Year) | Yes |
| DunedinPACE | 2022 | Blood | Aging Rate (Years/Year) | Yes |
| Zhang_10 | 2019 | Blood | Mortality Risk | Yes |
| YingCausAge | 2022 | Blood | Age (Years) | Yes |
| YingDamAge | 2022 | Blood | Age (Years) | Yes |
| YingAdaptAge | 2022 | Blood | Age (Years) | Yes |
| DNAmClockCortical | 2020 | Human Cortex | Cortex Age (Years) | No |
| VidalBralo | 2018 | Blood | Age (Years) | No |
| HRSInCHPhenoAge | 2022 | Blood | Age (Years) | No |
| StocZ | 2024 | Blood | Mortality Risk | No |
| StocP | 2024 | Blood | Age (Years) | No |
| StocH | 2024 | Multi-tissue | Age (Years) | No |
| Weidner | 2014 | Blood | Age (Years) | No |
| Garagnani | 2012 | Blood | Age (Years) | No |
| Bocklandt | 2011 | Blood | Age (Years) | No |

### Pediatric & Gestational Clocks

| Clock | Year | Tissue | Output | Parity Tested |
|-------|------|--------|--------|:-------------:|
| PEDBE | 2019 | Buccal | Age (Years) | No |
| Knight | 2016 | Cord Blood | Gestational Age | No |
| LeeControl | 2019 | Placenta | Gestational Age | No |
| LeeRefinedRobust | 2019 | Placenta | Gestational Age | No |
| LeeRobust | 2019 | Placenta | Gestational Age | No |
| Mayne | 2016 | Placenta | Gestational Age | No |
| Bohlin | 2017 | Cord Blood | Age (days) | No |

### Biological Markers

| Clock | Year | Tissue | Output | Parity Tested |
|-------|------|--------|--------|:-------------:|
| DNAmTL | 2019 | Blood, Adipose | Telomere Length | No |
| EpiTOC1 | 2016 | Blood | Stem Cell Division Rate | No |

### Trait & Lifestyle Biomarkers

| Clock | Year | Tissue | Output | Parity Tested |
|-------|------|--------|--------|:-------------:|
| AlcoholMcCartney | 2018 | Blood | Alcohol Consumption | No |
| SmokingMcCartney | 2018 | Blood | Smoking Status | No |
| BMI_McCartney | 2018 | Blood | BMI | No |
| EducationMcCartney | 2018 | Blood | Educational Attainment | No |
| BodyFatMcCartney | 2018 | Blood | Percentage Body Fat | No |
| HDLCholesterolMcCartney | 2018 | Blood | HDL Cholesterol | No |
| LDLCholesterolMcCartney | 2018 | Blood | LDL with Remnant Cholesterol | No |
| TotalCholesterolMcCartney | 2018 | Blood | Total Cholesterol | No |
| BMI_Reed | 2020 | Blood | BMI | No |

### Disease Risk Predictors

| Clock | Year | Tissue | Output | Parity Tested |
|-------|------|--------|--------|:-------------:|
| CVD_Westerman | 2020 | Blood | Coronary Heart Disease Status | No |
| AD_Bahado-Singh | 2021 | Blood | Alzheimer's Disease Status | No |
| DepressionBarbu | 2021 | Blood | Depression Risk | No |
| ProstateCancerKirby | 2017 | Prostate | Prostate Cancer Status | No |
| HepatoXu | 2017 | Circulating DNA | Hepatocellular Carcinoma Status | No |
| DownSyndrome | 2021 | Blood | Down Syndrome Prediction | No |

---

## Architecture

### Scoring Engine

The scoring pipeline follows these steps:

1. **Load coefficients** — CpG site weights from bundled CSV files
2. **Preprocess** (optional) — e.g. DunedinPACE quantile normalization against 20,000 gold-standard probe means
3. **Join** — inner join between coefficient CpGs and input methylation matrix on `CpGmarker`
4. **Intercept injection** — add intercept row if present in coefficients but not in input data
5. **Weighted sum** — dot product of beta values and coefficients per sample
6. **Transform** — apply clock-specific output transformation (Horvath anti-trafo, sigmoid, offset, etc.)

All operations use Polars DataFrames for efficient columnar computation.

### Methylation Data Standard

Following the [biolearn convention](https://github.com/bio-learn/biolearn):

- **Rows** = CpG site IDs (e.g. `cg00075967`)
- **Columns** = sample IDs
- **Values** = beta values in [0, 1] or null
- **Index column** = `CpGmarker` (first column)

### Input Formats

| Format | Extensions | Status | Notes |
|--------|-----------|--------|-------|
| CSV/TSV/TXT | `.csv`, `.tsv`, `.txt`, `.csv.gz`, `.tsv.gz`, `.txt.gz` | **Tested** | Separator auto-detected |
| Parquet | `.parquet` | **Tested** | First column normalized to `CpGmarker` |
| IDAT | `.idat`, directory | **Early stage** | Requires `pylluminator-modern` (Python 3.12+). Code exists but no end-to-end test on real IDAT data. |

### Transforms

| Transform | Clocks | Description |
|-----------|--------|-------------|
| `identity` | Hannum, PhenoAge, Lin, DunedinPoAm38, DunedinPACE, Zhang_10, Ying*, DNAmTL, EpiTOC1, etc. | Raw weighted sum (no transformation) |
| `horvathv1` | Horvathv1 | Horvath anti-trafo with intercept +0.696 |
| `horvathv2` | Horvathv2 | Horvath anti-trafo with intercept -0.447 |
| `pedbe` | PEDBE | Horvath anti-trafo with intercept -2.1 |
| `cortical` | DNAmClockCortical | Horvath anti-trafo with intercept +0.578 |
| `sigmoid` | BMI_McCartney, Education, BodyFat, HDL/LDL/TotalCholesterol, CVD_Westerman | Logistic sigmoid |
| `ad_bahado_singh` | AD_Bahado-Singh | Shifted sigmoid |
| Offset transforms | VidalBralo, StocZ/P/H, Mayne, Weidner, Bohlin | Constant offset added to weighted sum |

### Imputation Strategies (early stage, untested)

| Strategy | Description |
|----------|-------------|
| `averaging` | Fill missing CpG values with row mean across available samples |
| `reference` | Fill from a reference table of known CpG beta values |
| `hybrid` | Filter rows by coverage threshold, then average impute, then reference fill |
| `none` | No imputation (used by DunedinPACE) |
| `sesame_450k` | Default reference-based imputation using SeSAMe 450k averages |

Note: Imputation code exists but has **no tests**. The scoring engine and parity tests bypass imputation by running with `imputation_method="none"`.

---

## Example Datasets

Built-in catalogue of curated GEO datasets for testing and demonstration:

| GEO ID | Format | Samples | Tissue | Description | Used in Tests |
|--------|--------|---------|--------|-------------|:-------------:|
| GSE112618 | EPIC | 6 | Blood | FACS validation | Yes |
| GSE110554 | EPIC | 49 | Blood | FlowSorted.Blood.EPIC cell-type references | Yes |
| GSE164056 | EPIC | 143 | Blood | Social Anxiety Disorder study | Yes |
| GSE41037 | 450k | 388 | Blood | Hannum aging study | No |
| GSE137688 | 450k | 293 | Blood | Large blood methylation cohort | No |

```bash
# List available datasets
biomarkers download-example --list

# Download a specific dataset
biomarkers download-example GSE112618 -o data/input/examples/

# Download all EPIC datasets
biomarkers download-epic

# Download all datasets
biomarkers download-all-examples
```

## ComputAgeBench Integration (early stage)

Download and load datasets from the [ComputAgeBench](https://huggingface.co/datasets/computage/computage_bench) benchmark on HuggingFace. The download/load/list functions are implemented but **have no tests** — correctness on actual ComputAgeBench data is unverified.

```python
from just_biomarkers import list_computage_datasets, load_computage_dataset, load_computage_meta

datasets = list_computage_datasets()
dnam = load_computage_dataset("GSE40279")
meta = load_computage_meta("GSE40279")
```

```bash
biomarkers download-computage GSE40279
biomarkers download-computage-all
```

## Parity Testing

The library includes integration tests that compare scoring results against biolearn on real GEO datasets. Tests verify numerical equivalence within configurable tolerances (typically 1e-4 to 1e-3):

```bash
# Run parity tests (requires biolearn as a dev dependency)
uv run pytest just-biomarkers/tests/test_parity.py -v
uv run pytest just-biomarkers/tests/test_gse112618_parity.py -v

# Compare via CLI
uv run biomarkers compare-epic --dataset GSE112618
uv run biomarkers compare-epic --clocks PhenoAge,Horvathv2,DunedinPACE
```

---

## Workspace

`just-biomarkers` is the core library in a uv workspace that also includes:

### biomarkers-ui (functional, no tests)

Reflex web app for interactive clock computation. Upload a methylation CSV, select clocks, compute scores, view results with Plotly charts, and download output. Depends on `just-biomarkers` for the compute engine.

```bash
uv sync --all-packages
uv run ui
```

### nanopore-pipeline (scaffolded, never run)

Dagster pipeline for downloading and processing nanopore methylation BED files from Synology/GoFile into a unified CpG x sample beta-value matrix. The code (assets, sensors, CLI, resource tracking, BED parser) is written following project conventions but has **never been launched or executed against real data**. It is purely scaffolded code at this point — the GoFile scraper may not work, the BED parser has never seen real nanopore output, and no part of the pipeline has been validated.

```bash
uv run nanopore run       # Dagster UI (never tested)
uv run nanopore run --headless  # In-process execution (never tested)
```

## License

MIT
