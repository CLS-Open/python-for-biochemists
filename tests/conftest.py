"""Shared fixtures and data factories.

Every dataset used in the tests is generated here from a seeded ``Generator`` and a
stated error model. No test depends on a file on disk, so the suite cannot pass
because a stale CSV happened to match — and no test depends on a hidden global RNG,
so a failure is reproducible from its seed alone.
"""

from __future__ import annotations

import numpy as np
import pytest

# Truth values used throughout the suite. Chosen to be non-round-ish but easy to
# recognise in a failure message.
VMAX_TRUE = 10.0  # µmol min^-1 (mg protein)^-1
KM_TRUE = 25.0  # µM
KI_TRUE = 8.0  # µM


@pytest.fixture
def mm_data():
    """Factory: Michaelis-Menten data with constant-variance Gaussian noise."""

    def make(
        n_points: int = 7,
        noise: float = 0.20,
        seed: int = 7,
        vmax: float = VMAX_TRUE,
        km: float = KM_TRUE,
    ):
        S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])[:n_points]
        rng = np.random.default_rng(seed)
        v = vmax * S / (km + S) + noise * rng.standard_normal(S.size)
        return S, v

    return make


@pytest.fixture
def inhibition_rows():
    """Factory: a full (S, I, v) series from one of the three mechanisms."""

    def make(
        mechanism: str = "competitive", noise: float = 0.15, seed: int = 11, ki: float = KI_TRUE
    ):
        from kinetics.models import INHIBITION_MODELS

        func = INHIBITION_MODELS[mechanism]["func"]
        rng = np.random.default_rng(seed)
        rows = []
        for S in (2.0, 5.0, 10.0, 25.0, 100.0):
            for I in (0.0, 4.0, 16.0):
                clean = func(S, VMAX_TRUE, KM_TRUE, I, ki)
                rows.append({"S": S, "I": I, "v": float(clean + noise * rng.standard_normal())})
        return rows

    return make


@pytest.fixture
def calibration_data():
    """Factory: standards with proportional (5 %) multiplicative error."""

    def make(
        slope: float = 0.5,
        intercept: float = 0.01,
        noise: float = 0.05,
        seed: int = 3,
        concentrations=(0.5, 1.0, 2.0, 5.0, 10.0, 20.0),
    ):
        x = np.asarray(concentrations, dtype=float)
        rng = np.random.default_rng(seed)
        y = slope * x + intercept
        y = y * (1.0 + noise * rng.standard_normal(x.size))
        return x, y

    return make


@pytest.fixture
def write_csv(tmp_path):
    """Write a CSV from a header and rows, returning its path as a string."""

    def make(name: str, header: list[str], rows: list[list[float]]) -> str:
        path = tmp_path / name
        lines = [",".join(header)] + [",".join(str(v) for v in row) for row in rows]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return str(path)

    return make
