"""just-biomarkers: Polars-based methylation clock computation engine."""
from importlib.metadata import metadata

from just_biomarkers.computage_download import (
    download_computage_dataset,
    download_computage_meta,
    list_computage_datasets,
    load_computage_dataset,
    load_computage_meta,
)
from just_biomarkers.geo_download import (
    GEO_EXAMPLE_DATASETS,
    download_geo_example,
)
from just_biomarkers.imputation import (
    hybrid_impute,
    impute_from_average,
    impute_from_reference,
)
from just_biomarkers.io import (
    read_methylation_csv,
    read_methylation_idat,
    read_methylation_matrix,
    read_methylation_parquet,
    validate_methylation_matrix,
)
from just_biomarkers.models import BatchClockResult, ClockDefinition, ClockResult
from just_biomarkers.registry import (
    CLOCK_DEFINITIONS,
    get_clock,
    list_clocks,
    search_clocks,
)
from just_biomarkers.scoring import score_clock, score_clocks

_meta = metadata("just-biomarkers")
__version__: str = _meta["Version"]
__package_name__: str = _meta["Name"]
