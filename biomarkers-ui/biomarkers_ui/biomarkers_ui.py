"""Main Reflex app for methylation clock computation."""
from __future__ import annotations

import reflex as rx

from biomarkers_ui.state import AppState


def clock_selector() -> rx.Component:
    """Clock selection panel with checkboxes."""
    return rx.vstack(
        rx.hstack(
            rx.heading("Select Clocks", size="4"),
            rx.button("All", on_click=AppState.select_all_clocks, size="1", variant="outline"),
            rx.button("None", on_click=AppState.deselect_all_clocks, size="1", variant="outline"),
            rx.badge(f"{AppState.selected_clock_names.length()} selected"),
            align="center",
            spacing="2",
        ),
        rx.hstack(
            rx.text("Min match rate %", size="2"),
            rx.input(
                value=AppState.min_match_rate_percent,
                on_change=AppState.set_min_match_rate_percent,
                width="90px",
            ),
            rx.text("Policy", size="2"),
            rx.select(
                ["none", "hide", "fail"],
                value=AppState.match_rate_policy,
                on_change=AppState.set_match_rate_policy,
                size="2",
            ),
            spacing="2",
            align="center",
        ),
        rx.box(
            rx.foreach(
                AppState.clock_options,
                lambda clock: rx.hstack(
                    rx.checkbox(
                        clock["name"],
                        checked=AppState.selected_clock_names.contains(clock["name"]),
                        on_change=lambda _val: AppState.toggle_clock(clock["name"]),
                    ),
                    rx.text(clock["tissue"], color="gray", size="1"),
                    rx.text(clock["output"], color="gray", size="1"),
                    spacing="2",
                    align="center",
                ),
            ),
            max_height="400px",
            overflow_y="auto",
            border="1px solid var(--gray-6)",
            border_radius="8px",
            padding="3",
            width="100%",
        ),
        spacing="3",
        width="100%",
    )


def example_datasets_section() -> rx.Component:
    """Panel listing downloadable GEO example datasets."""
    return rx.vstack(
        rx.heading("Or Load an Example Dataset", size="4"),
        rx.text(
            "Click a dataset to download it from NCBI GEO and load it as input.",
            size="2",
            color="gray",
        ),
        rx.cond(
            AppState.has_examples,
            rx.vstack(
                rx.foreach(
                    AppState.available_examples,
                    lambda ex: rx.hstack(
                        rx.vstack(
                            rx.text(ex["id"], weight="bold", size="2"),
                            rx.text(ex["description"], size="1", color="gray"),
                            rx.hstack(
                                rx.badge(ex["format"], size="1", variant="surface"),
                                rx.badge(ex["samples"] + " samples", size="1", variant="surface"),
                                rx.badge(ex["tissue"], size="1", variant="surface"),
                                spacing="1",
                            ),
                            align_items="start",
                            spacing="1",
                        ),
                        rx.spacer(),
                        rx.button(
                            "Load",
                            on_click=AppState.load_geo_example(ex["id"]),
                            size="1",
                            variant="soft",
                            loading=AppState.is_downloading_example,
                            disabled=AppState.is_downloading_example,
                        ),
                        width="100%",
                        align="center",
                        padding="2",
                        border="1px solid var(--gray-5)",
                        border_radius="6px",
                    ),
                ),
                spacing="2",
                width="100%",
            ),
        ),
        spacing="3",
        width="100%",
    )


def upload_section() -> rx.Component:
    """File upload section."""
    return rx.vstack(
        rx.heading("Upload Methylation Matrix", size="4"),
        rx.hstack(
            rx.text("Input mode:", size="2"),
            rx.select(
                ["matrix", "idat"],
                value=AppState.input_mode,
                on_change=AppState.set_input_mode,
                size="2",
            ),
            rx.cond(
                AppState.input_mode == "idat",
                rx.hstack(
                    rx.text("IDAT preprocess:", size="2"),
                    rx.select(
                        ["raw", "noob"],
                        value=AppState.idat_preprocess,
                        on_change=AppState.set_idat_preprocess,
                        size="2",
                    ),
                    spacing="2",
                    align="center",
                ),
            ),
            spacing="3",
            align="center",
        ),
        rx.text(
            rx.cond(
                AppState.input_mode == "matrix",
                "Matrix mode: upload CSV/TSV/TXT/Parquet (rows=CpG, cols=samples).",
                "IDAT mode: upload matched *_Grn.idat + *_Red.idat files for each sample.",
            ),
            size="2",
            color="gray",
        ),
        rx.cond(
            AppState.input_mode == "matrix",
            rx.vstack(
                rx.upload(
                    rx.vstack(
                        rx.button("Select Matrix File", variant="outline"),
                        rx.text("or drag and drop here"),
                        align="center",
                        spacing="1",
                    ),
                    id="matrix_upload",
                    accept={
                        ".csv": ["text/csv"],
                        ".tsv": ["text/tab-separated-values"],
                        ".txt": ["text/plain"],
                        ".parquet": ["application/octet-stream"],
                    },
                    max_files=1,
                    border="2px dashed var(--gray-6)",
                    border_radius="8px",
                    padding="4",
                    width="100%",
                ),
                rx.button(
                    "Upload",
                    on_click=AppState.handle_upload(
                        rx.upload_files(upload_id="matrix_upload")
                    ),
                    size="2",
                ),
                spacing="2",
                width="100%",
            ),
            rx.vstack(
                rx.upload(
                    rx.vstack(
                        rx.button("Select IDAT Files", variant="outline"),
                        rx.text("select matching *_Grn.idat and *_Red.idat files"),
                        align="center",
                        spacing="1",
                    ),
                    id="idat_upload",
                    max_files=500,
                    border="2px dashed var(--gray-6)",
                    border_radius="8px",
                    padding="4",
                    width="100%",
                ),
                rx.button(
                    "Upload",
                    on_click=AppState.handle_upload(
                        rx.upload_files(upload_id="idat_upload")
                    ),
                    size="2",
                ),
                spacing="2",
                width="100%",
            ),
        ),
        rx.cond(
            AppState.has_upload,
            rx.callout(
                f"Loaded files: {AppState.upload_count}",
                icon="check",
                color_scheme="green",
                size="1",
            ),
        ),
        spacing="3",
        width="100%",
    )


def _sample_chart() -> rx.Component:
    """Bar chart for a single sample's clock scores using rx.plotly."""
    return rx.plotly(
        data=AppState.active_sample_figure,
        use_resize_handler=True,
        width="100%",
    )


def _sample_tab_content() -> rx.Component:
    """Content panel for the active sample tab: chart + scrollable table."""
    return rx.vstack(
        _sample_chart(),
        rx.box(
            rx.data_table(
                data=AppState.active_sample_rows,
                columns=["clock_name", "output", "score", "match_rate"],
                pagination=False,
                search=True,
                sort=True,
            ),
            width="100%",
            overflow_y="auto",
            max_height="500px",
            border="1px solid var(--gray-5)",
            border_radius="8px",
        ),
        spacing="3",
        width="100%",
    )


def results_section() -> rx.Component:
    """Results display section."""
    return rx.vstack(
        rx.cond(
            AppState.has_preprocessing,
            rx.vstack(
                rx.heading("Preprocessing Report", size="4"),
                rx.box(
                    rx.data_table(
                        data=AppState.preprocessing_rows,
                        columns=["metric", "value"],
                        pagination=False,
                        search=False,
                        sort=False,
                    ),
                    width="100%",
                    overflow_y="auto",
                    max_height="400px",
                    border="1px solid var(--gray-5)",
                    border_radius="8px",
                ),
                spacing="2",
                width="100%",
            ),
        ),
        rx.cond(
            AppState.has_low_match_rows,
            rx.vstack(
                rx.heading("Low Match-Rate Scores", size="4"),
                rx.text(
                    "These scores are below threshold and were hidden/failed by policy.",
                    size="2",
                    color="gray",
                ),
                rx.box(
                    rx.data_table(
                        data=AppState.low_match_rows,
                        columns=["sample_id", "clock_name", "output", "score", "match_rate"],
                        pagination=False,
                        search=True,
                        sort=True,
                    ),
                    width="100%",
                    overflow_y="auto",
                    max_height="400px",
                    border="1px solid var(--gray-5)",
                    border_radius="8px",
                ),
                spacing="2",
                width="100%",
            ),
        ),
        rx.cond(
            AppState.has_results,
            rx.vstack(
                rx.hstack(
                    rx.heading("Results", size="4"),
                    rx.badge(f"{AppState.result_count} scores"),
                    rx.button(
                        "Download CSV",
                        on_click=AppState.download_results_csv,
                        size="1",
                        variant="outline",
                    ),
                    align="center",
                    spacing="2",
                ),
                rx.tabs.root(
                    rx.tabs.list(
                        rx.foreach(
                            AppState.sample_names,
                            lambda name: rx.tabs.trigger(
                                name,
                                value=name,
                            ),
                        ),
                        overflow_x="auto",
                    ),
                    rx.tabs.content(
                        _sample_tab_content(),
                        value=AppState.active_sample_tab,
                    ),
                    value=AppState.active_sample_tab,
                    on_change=AppState.set_active_sample_tab,
                    width="100%",
                ),
                spacing="3",
                width="100%",
            ),
        ),
        rx.cond(
            AppState.has_clock_warnings,
            rx.vstack(
                rx.heading("Clock Warnings", size="4"),
                rx.box(
                    rx.foreach(
                        AppState.clock_warnings,
                        lambda w: rx.text(f"- {w}", size="2"),
                    ),
                    max_height="220px",
                    overflow_y="auto",
                    border="1px solid var(--gray-6)",
                    border_radius="8px",
                    padding="3",
                    width="100%",
                ),
                spacing="2",
                width="100%",
            ),
        ),
        spacing="4",
        width="100%",
    )


def index() -> rx.Component:
    """Main page layout."""
    return rx.box(
        rx.vstack(
            rx.hstack(
                rx.heading("Methylation Clock Calculator", size="6"),
                rx.cond(
                    AppState.status_message != "",
                    rx.badge(AppState.status_message, variant="surface"),
                ),
                align="center",
                spacing="3",
                width="100%",
            ),
            rx.cond(
                AppState.has_compatibility_warning,
                rx.callout(
                    AppState.compatibility_warning,
                    icon="triangle_alert",
                    color_scheme="red",
                    size="2",
                    width="100%",
                ),
            ),
            rx.separator(),
            rx.grid(
                rx.box(upload_section(), padding="3"),
                rx.box(clock_selector(), padding="3"),
                columns="2",
                spacing="4",
                width="100%",
            ),
            rx.box(example_datasets_section(), padding="3", width="100%"),
            rx.separator(),
            rx.hstack(
                rx.button(
                    "Compute Selected Clocks",
                    on_click=AppState.compute_clocks,
                    size="3",
                    loading=AppState.is_computing,
                    disabled=~AppState.has_upload,
                ),
                spacing="3",
            ),
            results_section(),
            spacing="4",
            padding="6",
            max_width="1200px",
            margin="0 auto",
        ),
    )


app = rx.App(
    theme=rx.theme(appearance="light", accent_color="teal"),
)
app.add_page(index, title="Methylation Clock Calculator", on_load=AppState.initialize)


def main() -> None:
    """Alternate entry point."""
    app.run()
