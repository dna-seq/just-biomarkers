"""Pydantic v2 models for methylation clock computation."""
from __future__ import annotations

from typing import Any, Callable, Optional

from pydantic import BaseModel, Field


class ClockDefinition(BaseModel):
    """Metadata and scoring configuration for a single methylation clock."""

    name: str
    year: int
    species: str = "Human"
    tissue: str = "Multi-tissue"
    source: str = ""
    output: str = "Age (Years)"
    coefficient_file: str
    intercept_name: str = "intercept"
    coefficient_column: str = "CoefficientTraining"
    transform_name: Optional[str] = None
    transform_offset: float = 0.0
    preprocess_name: Optional[str] = None
    default_imputation: str = "sesame_450k"

    model_config = {"frozen": True}


class ClockResult(BaseModel):
    """Result of a single clock scored against a methylation matrix."""

    clock_name: str
    sample_id: str
    score: float
    output: str = "Age (Years)"
    cpgs_matched: int
    cpgs_required: int
    match_rate: float = Field(ge=0.0, le=1.0)


class BatchClockResult(BaseModel):
    """Result of scoring multiple clocks and/or samples."""

    results: list[ClockResult] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    @property
    def scores_by_sample(self) -> dict[str, dict[str, float]]:
        out: dict[str, dict[str, float]] = {}
        for r in self.results:
            out.setdefault(r.sample_id, {})[r.clock_name] = r.score
        return out
