"""Reflex state for the biomarkers UI."""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Literal

import plotly.graph_objects as go

import reflex as rx

from just_biomarkers.registry import CLOCK_DEFINITIONS

InputMode = Literal["matrix", "idat"]
IdatPreprocess = Literal["raw", "noob"]
MatchRatePolicy = Literal["none", "hide", "fail"]


class AppState(rx.State):
    """Root application state for the methylation clock UI."""

    status_message: str = ""
    input_mode: InputMode = "matrix"
    is_downloading_example: bool = False
    available_examples: list[dict[str, str]] = []
    idat_preprocess: IdatPreprocess = "raw"
    min_match_rate_percent: float = 10.0
    match_rate_policy: MatchRatePolicy = "hide"

    uploaded_filenames: list[str] = []
    _uploaded_paths: list[str] = []
    _input_root_path: str = ""
    _idat_pairs_count: int = 0
    detected_array_type: str = "unknown"
    detected_probe_count: str = "unknown"
    compatibility_warning: str = ""

    clock_options: list[dict[str, str]] = []
    selected_clock_names: list[str] = []

    results_rows: list[dict] = []
    low_match_rows: list[dict] = []
    preprocessing_rows: list[dict[str, str]] = []
    clock_warnings: list[str] = []
    is_computing: bool = False
    active_sample_tab: str = ""
    active_sample_figure: go.Figure = go.Figure()

    def initialize(self) -> None:
        """Populate clock options and example dataset list on page load."""
        self.clock_options = [
            {
                "name": defn.name,
                "year": str(defn.year),
                "tissue": defn.tissue,
                "output": defn.output,
            }
            for defn in sorted(CLOCK_DEFINITIONS.values(), key=lambda d: d.name)
        ]
        self.selected_clock_names = [d["name"] for d in self.clock_options]

        from just_biomarkers.geo_download import GEO_EXAMPLE_DATASETS

        self.available_examples = [
            {
                "id": ds_id,
                "description": info["description"],
                "format": info["format"],
                "samples": info["samples"],
                "tissue": info["tissue"],
            }
            for ds_id, info in GEO_EXAMPLE_DATASETS.items()
        ]

    def set_input_mode(self, value: str) -> None:
        """Set active input mode."""
        if value in {"matrix", "idat"}:
            self.input_mode = value
            self.status_message = f"Input mode: {value}"

    def set_idat_preprocess(self, value: str) -> None:
        """Set IDAT preprocessing mode."""
        if value in {"raw", "noob"}:
            self.idat_preprocess = value
            self.status_message = f"IDAT preprocessing: {value}"

    def set_match_rate_policy(self, value: str) -> None:
        """Set low-match handling policy."""
        if value in {"none", "hide", "fail"}:
            self.match_rate_policy = value
            self.status_message = f"Match-rate policy: {value}"

    def set_min_match_rate_percent(self, value: str) -> None:
        """Set minimum allowed match rate (0..100)."""
        try:
            parsed = float(value)
        except ValueError:
            return
        if parsed < 0:
            parsed = 0.0
        if parsed > 100:
            parsed = 100.0
        self.min_match_rate_percent = parsed

    def toggle_clock(self, clock_name: str) -> None:
        if clock_name in self.selected_clock_names:
            self.selected_clock_names = [
                n for n in self.selected_clock_names if n != clock_name
            ]
        else:
            self.selected_clock_names = self.selected_clock_names + [clock_name]

    def select_all_clocks(self) -> None:
        self.selected_clock_names = [d["name"] for d in self.clock_options]

    def deselect_all_clocks(self) -> None:
        self.selected_clock_names = []

    async def load_geo_example(self, dataset_id: str) -> None:
        """Download a GEO example dataset from NCBI and set it as the current input."""
        self.is_downloading_example = True
        self.status_message = f"Downloading {dataset_id} from NCBI GEO…"
        self.results_rows = []
        self.low_match_rows = []
        self.clock_warnings = []
        self.preprocessing_rows = []
        yield

        from just_biomarkers.geo_download import download_geo_example

        output_dir = Path("data/input/examples")
        dest = download_geo_example(dataset_id, output_dir)

        self._input_root_path = str(dest)
        self._uploaded_paths = [str(dest)]
        self.uploaded_filenames = [dest.name]
        self.input_mode = "matrix"
        self.status_message = f"Loaded example: {dest.name}"
        self.is_downloading_example = False

    @staticmethod
    def _is_idat_file(filename: str) -> bool:
        lower = filename.lower()
        return lower.endswith(".idat") or lower.endswith(".idat.gz")

    @staticmethod
    def _idat_base_name(filename: str) -> str:
        lower = filename.lower()
        if lower.endswith(".idat.gz"):
            name = filename[:-8]
        elif lower.endswith(".idat"):
            name = filename[:-5]
        else:
            return filename
        if name.endswith("_Grn") or name.endswith("_Red"):
            return name[:-4]
        return name

    @staticmethod
    def _count_idat_pairs(filenames: list[str]) -> tuple[int, list[str]]:
        channels: dict[str, set[str]] = {}
        for filename in filenames:
            if not AppState._is_idat_file(filename):
                continue
            lower = filename.lower()
            base = AppState._idat_base_name(filename)
            channel = "grn" if "_grn.idat" in lower else "red" if "_red.idat" in lower else ""
            if not channel:
                continue
            if base not in channels:
                channels[base] = set()
            channels[base].add(channel)
        complete = sum(1 for c in channels.values() if {"grn", "red"}.issubset(c))
        incomplete = sorted(base for base, c in channels.items() if c != {"grn", "red"})
        return complete, incomplete

    @staticmethod
    def _detect_idat_platform(input_root: str) -> tuple[str, str]:
        """Detect IDAT platform using first available green IDAT file."""
        root = Path(input_root)
        green_files = sorted(root.glob("*_Grn.idat")) + sorted(root.glob("*_Grn.idat.gz"))
        if not green_files:
            return "unknown", "unknown"

        try:
            from pylluminator.annotations import detect_array  # type: ignore
            from pylluminator.read_idat import IdatDataset  # type: ignore
        except Exception:
            return "unknown", "unknown"

        try:
            idat = IdatDataset(str(green_files[0]))
            probe_count = str(int(idat.n_snps_read))
            array_type = str(detect_array(int(idat.n_snps_read)))
            return array_type, probe_count
        except Exception:
            return "unknown", "unknown"

    async def handle_upload(self, files: list[rx.UploadFile]) -> None:
        """Handle matrix or IDAT uploads."""
        if not files:
            self.status_message = "No files selected for upload."
            return
        self.results_rows = []
        self.low_match_rows = []
        self.clock_warnings = []
        self.preprocessing_rows = []
        self.detected_array_type = "unknown"
        self.detected_probe_count = "unknown"
        self.compatibility_warning = ""

        if self.input_mode == "matrix":
            upload_file = files[0]
            upload_dir = Path("data/input/matrix")
            upload_dir.mkdir(parents=True, exist_ok=True)
            dest = upload_dir / upload_file.filename
            content = await upload_file.read()
            dest.write_bytes(content)

            self.uploaded_filenames = [upload_file.filename]
            self._uploaded_paths = [str(dest)]
            self._input_root_path = str(dest)
            self._idat_pairs_count = 0
            self.status_message = f"Uploaded matrix file: {upload_file.filename}"
            return

        upload_root = Path("data/input/idat") / datetime.now().strftime("%Y%m%d_%H%M%S")
        upload_root.mkdir(parents=True, exist_ok=True)
        saved_filenames: list[str] = []
        saved_paths: list[str] = []
        for upload_file in files:
            dest = upload_root / upload_file.filename
            content = await upload_file.read()
            dest.write_bytes(content)
            saved_filenames.append(upload_file.filename)
            saved_paths.append(str(dest))

        complete_pairs, incomplete_bases = self._count_idat_pairs(saved_filenames)
        if complete_pairs == 0:
            self.status_message = (
                "No complete IDAT pairs found. Upload matching *_Grn.idat and *_Red.idat files."
            )
            self.uploaded_filenames = saved_filenames
            self._uploaded_paths = saved_paths
            self._input_root_path = str(upload_root)
            self._idat_pairs_count = 0
            return

        self.uploaded_filenames = saved_filenames
        self._uploaded_paths = saved_paths
        self._input_root_path = str(upload_root)
        self._idat_pairs_count = complete_pairs
        detected_array_type, idat_probe_count = self._detect_idat_platform(
            str(upload_root)
        )
        self.detected_array_type = detected_array_type
        self.detected_probe_count = idat_probe_count
        self.compatibility_warning = (
            "Platform compatibility warning: detected Mammal40 IDAT platform, "
            "but current clock registry is human-focused and not Mammal40-specific. "
            "Expect very low match rates and potentially unreliable scores."
            if detected_array_type == "Mammal40"
            else ""
        )
        self.status_message = (
            f"Uploaded {len(saved_filenames)} files, detected {complete_pairs} IDAT pair(s)"
            + (f", {len(incomplete_bases)} incomplete pair(s)" if incomplete_bases else "")
        )

    def compute_clocks(self) -> None:
        """Run selected clocks against uploaded matrix or IDAT input."""
        if not self._input_root_path:
            self.status_message = "Please upload input files first"
            return
        if not self.selected_clock_names:
            self.status_message = "Please select at least one clock"
            return

        self.is_computing = True
        self.status_message = "Computing..."
        self.clock_warnings = []
        self.preprocessing_rows = []
        self.low_match_rows = []
        self.compatibility_warning = ""

        from just_biomarkers.io import read_methylation_matrix, validate_methylation_matrix
        from just_biomarkers.scoring import score_clocks

        input_path = self._input_root_path
        try:
            dnam = read_methylation_matrix(
                input_path,
                idat_preprocess=self.idat_preprocess,
            )
            input_warnings = validate_methylation_matrix(dnam)
            batch = score_clocks(dnam, clock_names=self.selected_clock_names)
        except Exception as exc:
            self.status_message = f"Compute failed: {exc}"
            self.is_computing = False
            return

        threshold = self.min_match_rate_percent / 100.0
        low_match_results = [r for r in batch.results if r.match_rate < threshold]
        kept_results = [r for r in batch.results if r.match_rate >= threshold]

        if self.match_rate_policy == "none":
            selected_results = batch.results
        elif self.match_rate_policy == "hide":
            selected_results = kept_results
            if low_match_results:
                self.clock_warnings = self.clock_warnings + [
                    f"{len(low_match_results)} clock scores hidden due to match rate < {self.min_match_rate_percent:.1f}%."
                ]
        else:
            selected_results = []
            if low_match_results:
                self.clock_warnings = self.clock_warnings + [
                    f"Rejected because {len(low_match_results)} clocks are below {self.min_match_rate_percent:.1f}% match rate."
                ]

        self.low_match_rows = [
            {
                "sample_id": r.sample_id,
                "clock_name": r.clock_name,
                "output": r.output,
                "score": round(r.score, 4),
                "match_rate": f"{r.match_rate:.1%}",
            }
            for r in low_match_results
        ]

        self.results_rows = [
            {
                "sample_id": r.sample_id,
                "clock_name": r.clock_name,
                "output": r.output,
                "score": round(r.score, 4),
                "match_rate": f"{r.match_rate:.1%}",
            }
            for r in selected_results
        ]
        self.clock_warnings = batch.warnings + self.clock_warnings

        seen: set[str] = set()
        first_sample = ""
        for r in selected_results:
            if r.sample_id not in seen:
                if not first_sample:
                    first_sample = r.sample_id
                seen.add(r.sample_id)
        self.active_sample_tab = first_sample
        self.active_sample_figure = self._build_figure_for_sample(first_sample)

        n_samples = max(0, len(dnam.columns) - 1)
        detected_array_type = "n/a"
        idat_probe_count = "n/a"
        if self.input_mode == "idat":
            detected_array_type, idat_probe_count = self._detect_idat_platform(
                input_path
            )
            self.detected_array_type = detected_array_type
            self.detected_probe_count = idat_probe_count
            if detected_array_type == "Mammal40":
                self.compatibility_warning = (
                    "Platform compatibility warning: detected Mammal40 IDAT platform, "
                    "but current clock registry is human-focused and not Mammal40-specific. "
                    "Expect very low match rates and potentially unreliable scores."
                )
                self.clock_warnings = self.clock_warnings + [
                    "Detected platform Mammal40: current registry has no Mammal40-specific clocks, so low match rates are expected."
                ]

        self.preprocessing_rows = [
            {"metric": "input_mode", "value": self.input_mode},
            {"metric": "idat_preprocess", "value": self.idat_preprocess if self.input_mode == "idat" else "n/a"},
            {"metric": "detected_array_type", "value": detected_array_type},
            {"metric": "idat_probe_count", "value": idat_probe_count},
            {"metric": "min_match_rate_percent", "value": f"{self.min_match_rate_percent:.1f}"},
            {"metric": "match_rate_policy", "value": self.match_rate_policy},
            {"metric": "uploaded_files", "value": str(len(self.uploaded_filenames))},
            {"metric": "idat_pairs", "value": str(self._idat_pairs_count) if self.input_mode == "idat" else "n/a"},
            {"metric": "matrix_cpg_rows", "value": str(dnam.height)},
            {"metric": "matrix_samples", "value": str(n_samples)},
            {"metric": "total_clock_scores", "value": str(len(batch.results))},
            {"metric": "low_match_scores", "value": str(len(low_match_results))},
            {"metric": "returned_scores", "value": str(len(self.results_rows))},
            {"metric": "input_path", "value": input_path},
        ]

        warning_count = len(self.clock_warnings) + len(input_warnings)
        n_results = len(self.results_rows)
        if self.match_rate_policy == "fail" and low_match_results:
            self.status_message = (
                f"Failed: {len(low_match_results)} scores below "
                f"{self.min_match_rate_percent:.1f}% match rate"
                + (f" ({warning_count} warnings)" if warning_count else "")
            )
            self.results_rows = []
            self.is_computing = False
            return

        self.status_message = (
            f"Done: {n_results} scores computed"
            + (f" ({warning_count} warnings)" if warning_count else "")
        )
        self.is_computing = False

    def download_results_csv(self) -> rx.Component:
        """Generate a CSV download of results."""
        import io
        import polars as pl

        if not self.results_rows:
            return rx.download(data="", filename="empty.csv")
        df = pl.DataFrame(self.results_rows)
        buf = io.StringIO()
        df.write_csv(buf)
        return rx.download(data=buf.getvalue(), filename="clock_results.csv")

    def set_active_sample_tab(self, tab: str) -> None:
        """Set the currently visible sample tab."""
        self.active_sample_tab = tab
        self.active_sample_figure = self._build_figure_for_sample(tab)

    def _build_figure_for_sample(self, sample_id: str) -> go.Figure:
        from plotly.subplots import make_subplots

        rows = [r for r in self.results_rows if r["sample_id"] == sample_id]
        if not rows:
            return go.Figure()

        # Group clocks by output unit — each group gets its own subplot row
        groups: dict[str, list[dict]] = {}
        for r in rows:
            groups.setdefault(r["output"], []).append(r)

        n_groups = len(groups)
        # Height per subplot row scales with number of clocks in it
        row_heights = [max(len(g), 1) for g in groups.values()]
        total_weight = sum(row_heights)
        row_fractions = [h / total_weight for h in row_heights]

        vertical_spacing = min(0.08, 0.9 / max(n_groups - 1, 1))
        fig = make_subplots(
            rows=n_groups,
            cols=1,
            row_heights=row_fractions,
            shared_xaxes=False,
            subplot_titles=list(groups.keys()),
            vertical_spacing=vertical_spacing,
        )

        teal = "#00897b"
        for row_idx, (output_label, group_rows) in enumerate(groups.items(), start=1):
            clock_names = [r["clock_name"] for r in group_rows]
            scores = [float(r["score"]) for r in group_rows]
            match_rates = [r["match_rate"] for r in group_rows]
            fig.add_trace(
                go.Bar(
                    x=clock_names,
                    y=scores,
                    marker_color=teal,
                    customdata=list(zip(match_rates, [output_label] * len(clock_names))),
                    hovertemplate=(
                        "<b>%{x}</b><br>"
                        f"Output: {output_label}<br>"
                        "Score: %{y:.4f}<br>"
                        "Match rate: %{customdata[0]}"
                        "<extra></extra>"
                    ),
                    showlegend=False,
                ),
                row=row_idx,
                col=1,
            )
            fig.update_xaxes(tickangle=-45, automargin=True, row=row_idx, col=1)
            fig.update_yaxes(title_text=output_label, row=row_idx, col=1)

        # Dynamic total height: 220px base per group + 18px per bar
        n_bars_total = sum(len(g) for g in groups.values())
        chart_height = max(400, n_groups * 220 + n_bars_total * 18)

        fig.update_layout(
            height=chart_height,
            margin={"t": 40, "b": 60, "l": 80, "r": 20},
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    @rx.var
    def has_results(self) -> bool:
        return len(self.results_rows) > 0

    @rx.var
    def sample_names(self) -> list[str]:
        seen: list[str] = []
        seen_set: set[str] = set()
        for row in self.results_rows:
            sid = row["sample_id"]
            if sid not in seen_set:
                seen.append(sid)
                seen_set.add(sid)
        return seen

    @rx.var
    def active_sample_rows(self) -> list[dict]:
        return [r for r in self.results_rows if r["sample_id"] == self.active_sample_tab]

    @rx.var
    def has_examples(self) -> bool:
        return len(self.available_examples) > 0

    @rx.var
    def has_upload(self) -> bool:
        return len(self._uploaded_paths) > 0

    @rx.var
    def has_preprocessing(self) -> bool:
        return len(self.preprocessing_rows) > 0

    @rx.var
    def has_clock_warnings(self) -> bool:
        return len(self.clock_warnings) > 0

    @rx.var
    def has_low_match_rows(self) -> bool:
        return len(self.low_match_rows) > 0

    @rx.var
    def has_compatibility_warning(self) -> bool:
        return self.compatibility_warning != ""

    @rx.var
    def result_count(self) -> int:
        return len(self.results_rows)

    @rx.var
    def upload_count(self) -> int:
        return len(self.uploaded_filenames)
