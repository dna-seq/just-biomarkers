"""Tests for methylation input readers."""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from just_biomarkers.io import read_methylation_matrix


def test_read_methylation_matrix_csv_renames_first_column(tmp_path: Path) -> None:
    csv_path = tmp_path / "matrix.csv"
    csv_path.write_text("probe,s1,s2\ncg0001,0.1,0.2\ncg0002,0.3,0.4\n", encoding="utf-8")

    df = read_methylation_matrix(csv_path)

    assert df.columns[0] == "CpGmarker"
    assert df.shape == (2, 3)
    assert df["s1"].dtype == pl.Float64


def test_read_methylation_matrix_tsv(tmp_path: Path) -> None:
    tsv_path = tmp_path / "matrix.tsv"
    tsv_path.write_text(
        "CpGmarker\ts1\ts2\ncg0001\t0.1\t0.2\ncg0002\t0.3\t0.4\n",
        encoding="utf-8",
    )

    df = read_methylation_matrix(tsv_path)

    assert df.columns == ["CpGmarker", "s1", "s2"]
    assert df["s2"].dtype == pl.Float64


def test_read_methylation_matrix_parquet(tmp_path: Path) -> None:
    parquet_path = tmp_path / "matrix.parquet"
    source = pl.DataFrame(
        {
            "marker": ["cg0001", "cg0002"],
            "sample_1": [0.25, 0.75],
        }
    )
    source.write_parquet(str(parquet_path))

    df = read_methylation_matrix(parquet_path)

    assert df.columns == ["CpGmarker", "sample_1"]
    assert df.shape == (2, 2)


def test_read_methylation_matrix_idat_file_reports_clear_error(tmp_path: Path) -> None:
    idat_path = tmp_path / "207435510088_R08C01_Grn.idat"
    idat_path.write_bytes(b"IDAT")

    with pytest.raises((ImportError, ValueError), match="pylluminator|IDAT"):
        read_methylation_matrix(idat_path)


def test_read_methylation_matrix_idat_directory_reports_clear_error(
    tmp_path: Path,
) -> None:
    (tmp_path / "207435510088_R08C01_Grn.idat").write_bytes(b"IDAT")
    (tmp_path / "207435510088_R08C01_Red.idat").write_bytes(b"IDAT")

    with pytest.raises((ImportError, ValueError), match="pylluminator|IDAT"):
        read_methylation_matrix(tmp_path)
