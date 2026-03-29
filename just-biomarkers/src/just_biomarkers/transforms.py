"""Transform functions used by methylation clocks.

Each transform takes a raw weighted-sum scalar/array and produces the final
prediction.  Named transforms are registered in TRANSFORM_REGISTRY so that
clock definitions can reference them by string key.
"""
from __future__ import annotations

from typing import Callable

import numpy as np


def identity(x: float | np.ndarray) -> float | np.ndarray:
    return x


def anti_trafo(x: float | np.ndarray, adult_age: float = 20.0) -> float | np.ndarray:
    """Horvath anti-transformation (inverse of the age transformation)."""
    return np.where(
        x < 0,
        (1 + adult_age) * np.exp(x) - 1,
        (1 + adult_age) * x + adult_age,
    )


def horvathv1_transform(x: float | np.ndarray) -> float | np.ndarray:
    return anti_trafo(x + 0.696)


def horvathv2_transform(x: float | np.ndarray) -> float | np.ndarray:
    return anti_trafo(x - 0.447119319)


def pedbe_transform(x: float | np.ndarray) -> float | np.ndarray:
    return anti_trafo(x - 2.1)


def cortical_transform(x: float | np.ndarray) -> float | np.ndarray:
    return anti_trafo(x + 0.577682570446177)


def sigmoid(x: float | np.ndarray) -> float | np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def ad_bahado_singh_transform(x: float | np.ndarray) -> float | np.ndarray:
    return 1.0 / (1.0 + np.exp(-x - 0.072))


def _offset_transform(offset: float) -> Callable:
    def _t(x: float | np.ndarray) -> float | np.ndarray:
        return x + offset
    return _t


TRANSFORM_REGISTRY: dict[str, Callable] = {
    "identity": identity,
    "horvathv1": horvathv1_transform,
    "horvathv2": horvathv2_transform,
    "pedbe": pedbe_transform,
    "cortical": cortical_transform,
    "sigmoid": sigmoid,
    "ad_bahado_singh": ad_bahado_singh_transform,
    "vidal_bralo": _offset_transform(84.7),
    "stocz": _offset_transform(64.8077188694894),
    "stocp": _offset_transform(92.8310813279039),
    "stoch": _offset_transform(59.8015666314217),
    "mayne": _offset_transform(24.99026),
    "weidner": _offset_transform(38.0),
    "bohlin": _offset_transform(277.2421),
}
