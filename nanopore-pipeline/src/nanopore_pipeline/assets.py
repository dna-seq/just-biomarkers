"""Dagster Software-Defined Assets for the nanopore methylation pipeline.

Pipeline stages:
1. synology_file_index — scrape the Synology/GoFile directory listing for .bed.gz files
2. nanopore_bed_files — download discovered BED files to local cache
3. nanopore_methylation_matrix — parse BED files into a unified CpG methylation matrix
"""

from pathlib import Path
from typing import Any

import dagster as dg
import polars as pl
from dagster import (
    AssetExecutionContext,
    MetadataValue,
    Output,
    SourceAsset,
    asset,
)

from nanopore_pipeline.resources import CacheDirResource, SourceURLResource
from nanopore_pipeline.runtime import resource_tracker

synology_nanopore_source = SourceAsset(
    key="synology_nanopore_source",
    group_name="external",
    description=(
        "Synology/GoFile directory hosting nanopore methylation BED files. "
        "URL: http://gofile.me/76Mv5/irFIggw48"
    ),
    metadata={"url": "http://gofile.me/76Mv5/irFIggw48"},
)


def _scrape_file_links(url: str) -> list[dict[str, str]]:
    """Scrape a Synology/GoFile directory listing for downloadable file links.

    Returns a list of dicts with keys: name, url, size (if available).
    """
    import requests
    from bs4 import BeautifulSoup

    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")

    files: list[dict[str, str]] = []
    for link in soup.find_all("a", href=True):
        href = link["href"]
        name = link.get_text(strip=True)
        if name and (
            name.endswith(".bed.gz")
            or name.endswith(".bed")
            or name.endswith(".tsv.gz")
            or name.endswith(".tsv")
            or name.endswith(".csv.gz")
            or name.endswith(".csv")
        ):
            full_url = href if href.startswith("http") else f"{url.rstrip('/')}/{href}"
            files.append({"name": name, "url": full_url})

    return files


@asset(
    group_name="download",
    description="Discover available nanopore methylation files from the Synology/GoFile directory listing.",
)
def synology_file_index(
    context: AssetExecutionContext,
    source_url_resource: SourceURLResource,
    cache_dir_resource: CacheDirResource,
) -> Output[Path]:
    """Scrape the directory listing and persist the file index as a parquet."""
    with resource_tracker("synology_file_index", context=context):
        url = source_url_resource.get_url()
        context.log.info(f"Scraping file index from: {url}")

        files = _scrape_file_links(url)
        context.log.info(f"Found {len(files)} downloadable files")

        if not files:
            context.log.warning(
                "No files found. The Synology/GoFile page may require authentication "
                "or the page structure may have changed. You can manually place .bed.gz "
                "files in the downloads directory and re-run."
            )

        index_df = pl.DataFrame(files) if files else pl.DataFrame(
            schema={"name": pl.Utf8, "url": pl.Utf8}
        )

        output_dir = cache_dir_resource.nanopore_dir()
        index_path = output_dir / "file_index.parquet"
        index_df.write_parquet(index_path)

        context.add_output_metadata({
            "source_url": MetadataValue.url(url),
            "n_files": MetadataValue.int(len(files)),
            "file_names": MetadataValue.text(
                "\n".join(f["name"] for f in files[:50])
            ),
            "index_path": MetadataValue.path(str(index_path)),
        })

    return Output(index_path)


def _download_file(url: str, dest: Path, context: AssetExecutionContext) -> bool:
    """Download a single file with progress logging. Returns True on success."""
    import requests

    if dest.exists() and dest.stat().st_size > 0:
        context.log.info(f"Already cached: {dest.name} ({dest.stat().st_size:,} bytes)")
        return True

    context.log.info(f"Downloading: {dest.name} from {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)

    resp = requests.get(url, stream=True, timeout=300)
    resp.raise_for_status()

    tmp = dest.with_suffix(dest.suffix + ".tmp")
    total = 0
    with open(tmp, "wb") as f:
        for chunk in resp.iter_content(chunk_size=8192):
            f.write(chunk)
            total += len(chunk)

    if total == 0:
        tmp.unlink(missing_ok=True)
        context.log.warning(f"Empty download for {dest.name}, skipping")
        return False

    tmp.rename(dest)
    context.log.info(f"Downloaded: {dest.name} ({total:,} bytes)")
    return True


@asset(
    group_name="download",
    deps=["synology_file_index"],
    description="Download nanopore methylation BED files to local cache.",
)
def nanopore_bed_files(
    context: AssetExecutionContext,
    cache_dir_resource: CacheDirResource,
) -> Output[Path]:
    """Download all BED files from the file index to the local cache directory."""
    with resource_tracker("nanopore_bed_files", context=context):
        index_path = cache_dir_resource.nanopore_dir() / "file_index.parquet"
        if not index_path.exists():
            context.log.warning("No file index found. Run synology_file_index first.")
            return Output(cache_dir_resource.downloads_dir())

        index_df = pl.read_parquet(index_path)
        downloads_dir = cache_dir_resource.downloads_dir()
        n_total = index_df.height
        n_ok = 0
        n_cached = 0
        n_failed = 0

        for row in index_df.iter_rows(named=True):
            name = row["name"]
            url = row["url"]
            dest = downloads_dir / name

            if dest.exists() and dest.stat().st_size > 0:
                n_cached += 1
                n_ok += 1
                continue

            success = _download_file(url, dest, context)
            if success:
                n_ok += 1
            else:
                n_failed += 1

        context.add_output_metadata({
            "n_total": MetadataValue.int(n_total),
            "n_ok": MetadataValue.int(n_ok),
            "n_cached": MetadataValue.int(n_cached),
            "n_failed": MetadataValue.int(n_failed),
            "downloads_dir": MetadataValue.path(str(downloads_dir)),
            "coverage_ratio": MetadataValue.float(n_ok / max(n_total, 1)),
        })

    return Output(downloads_dir)


def _parse_nanopore_bed(bed_path: Path, sample_id: str) -> pl.DataFrame:
    """Parse a single nanopore modkit BED file into a CpG-level DataFrame.

    Nanopore modkit BED files typically have columns:
    chrom, start, end, mod_code, score, strand, start2, end2, color,
    n_valid, fraction_modified, n_mod, n_canonical, n_other, n_delete, n_fail, ...

    We extract: chrom, start, strand, fraction_modified (as the beta value).
    The CpG marker ID is constructed as chrom:position.
    """
    import gzip

    open_fn = gzip.open if bed_path.suffix == ".gz" else open

    rows: list[dict[str, Any]] = []
    with open_fn(bed_path, "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue

            chrom = parts[0]
            start = int(parts[1])
            strand = parts[5] if len(parts) > 5 else "."
            n_valid = int(parts[9]) if len(parts) > 9 else 0
            fraction = float(parts[10]) / 100.0 if len(parts) > 10 else 0.0

            if n_valid < 1:
                continue

            cpg_id = f"{chrom}:{start}"
            rows.append({
                "CpGmarker": cpg_id,
                "chrom": chrom,
                "position": start,
                "strand": strand,
                "n_valid": n_valid,
                sample_id: fraction,
            })

    if not rows:
        return pl.DataFrame(
            schema={
                "CpGmarker": pl.Utf8,
                "chrom": pl.Utf8,
                "position": pl.Int64,
                "strand": pl.Utf8,
                "n_valid": pl.Int64,
                sample_id: pl.Float64,
            }
        )

    return pl.DataFrame(rows)


@asset(
    group_name="compute",
    deps=["nanopore_bed_files"],
    description="Parse downloaded BED files and build a unified CpG methylation beta-value matrix.",
)
def nanopore_methylation_matrix(
    context: AssetExecutionContext,
    cache_dir_resource: CacheDirResource,
) -> Output[Path]:
    """Build the methylation matrix from all downloaded BED files.

    The output matrix follows the biolearn methylation standard:
    rows = CpG site IDs, columns = sample IDs, values = beta values [0,1].
    """
    with resource_tracker("nanopore_methylation_matrix", context=context):
        downloads_dir = cache_dir_resource.downloads_dir()
        output_dir = cache_dir_resource.processed_dir()

        bed_files = sorted(
            p for p in downloads_dir.iterdir()
            if p.suffix in (".gz", ".bed") and "bed" in p.name
        )

        if not bed_files:
            context.log.warning(
                "No BED files found in downloads directory. "
                "Place .bed.gz files manually or run the download pipeline."
            )
            empty_path = output_dir / "nanopore_methylation_matrix.parquet"
            pl.DataFrame(schema={"CpGmarker": pl.Utf8}).write_parquet(empty_path)
            return Output(empty_path)

        context.log.info(f"Processing {len(bed_files)} BED files")

        merged: pl.DataFrame | None = None
        n_processed = 0

        for bed_path in bed_files:
            sample_id = bed_path.stem
            if sample_id.endswith(".bed"):
                sample_id = sample_id[:-4]

            context.log.info(f"Parsing: {bed_path.name} (sample: {sample_id})")

            sample_df = _parse_nanopore_bed(bed_path, sample_id)
            if sample_df.height == 0:
                context.log.warning(f"No CpG sites parsed from {bed_path.name}")
                continue

            beta_df = sample_df.select(["CpGmarker", sample_id])

            if merged is None:
                merged = beta_df
            else:
                merged = merged.join(beta_df, on="CpGmarker", how="full", coalesce=True)

            n_processed += 1
            context.log.info(
                f"  {sample_id}: {sample_df.height:,} CpG sites"
            )

        if merged is None:
            merged = pl.DataFrame(schema={"CpGmarker": pl.Utf8})

        output_path = output_dir / "nanopore_methylation_matrix.parquet"
        merged.write_parquet(output_path)

        n_cpgs = merged.height
        n_samples = len(merged.columns) - 1

        context.log.info(
            f"Methylation matrix: {n_cpgs:,} CpG sites x {n_samples} samples"
        )
        context.add_output_metadata({
            "n_cpg_sites": MetadataValue.int(n_cpgs),
            "n_samples": MetadataValue.int(n_samples),
            "n_bed_files_processed": MetadataValue.int(n_processed),
            "n_bed_files_total": MetadataValue.int(len(bed_files)),
            "output_path": MetadataValue.path(str(output_path)),
            "sample_ids": MetadataValue.text(
                "\n".join(c for c in merged.columns if c != "CpGmarker")
            ),
        })

    return Output(output_path)
