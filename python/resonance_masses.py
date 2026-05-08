#!/usr/bin/env python3
"""Compute the resonance mass m_a = omega_pl at the B-field maximum of each solar layer.

This script uses the SolarAxionFlux Python bindings when available. It searches for the
radius where the magnetic field is maximal in each layer, evaluates the plasma frequency
there, and reports the corresponding resonance mass.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = REPO_ROOT / "data/solar_models/SolarModel_B16-AGSS09.dat"

# These are the same layer scales used in include/solaxflux/constants.hpp.
RADIUS_CZ = 0.712
SIZE_TACH = 0.02
RADIUS_OUTER = 0.96
SIZE_OUTER = 0.035

# Default magnetic-field strengths used by SolarModel.
B_RAD_T = 3.0e3
B_TACH_T = 50.0
B_OUTER_T = 4.0

# Inverse of the factor used by SolarModel::bfield to convert Tesla -> keV^2.
TESLA_TO_KEV2 = 1.0 / (1.0e6 * math.sqrt(4.0 * math.pi) * 1.4440271e-3)

# Same profile-shape constants as in SolarModel::bfield.
LAMBDA = 10.0 * RADIUS_CZ + 1.0
LAMBDA_FACTOR = (1.0 + LAMBDA) * (1.0 + 1.0 / LAMBDA) ** LAMBDA


def load_pyaxionflux():
    """Load the compiled Python wrapper from a few common locations."""
    try:
        import pyaxionflux  # type: ignore

        return pyaxionflux
    except Exception:
        pass

    try:
        from lib import pyaxionflux  # type: ignore

        return pyaxionflux
    except Exception:
        pass

    search_roots = [REPO_ROOT / "build", REPO_ROOT]
    for root in search_roots:
        for so_path in root.rglob("pyaxionflux*.so"):
            spec = importlib.util.spec_from_file_location("pyaxionflux", so_path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module

    raise ImportError(
        "Could not import pyaxionflux. Build the project first, or add the compiled "
        "extension to PYTHONPATH."
    )


def bfield_keV2(r: float) -> float:
    """Magnetic field profile in the same units as SolarModel::bfield."""
    result = 0.0
    if r < RADIUS_CZ + SIZE_TACH:
        x = (r / RADIUS_CZ) ** 2
        if x < 1.0:
            result += B_RAD_T * LAMBDA_FACTOR * x * (1.0 - x) ** LAMBDA
        y = ((r - RADIUS_CZ) / SIZE_TACH) ** 2
        if y < 1.0:
            result += B_TACH_T * (1.0 - y)
    else:
        z = ((r - RADIUS_OUTER) / SIZE_OUTER) ** 2
        if z < 1.0:
            result += B_OUTER_T * (1.0 - z)

    return result * TESLA_TO_KEV2


def bfield_tesla(r: float) -> float:
    return bfield_keV2(r) / TESLA_TO_KEV2


def find_layer_maximum(func, r_lo: float, r_hi: float, points: int) -> tuple[float, float]:
    grid = np.linspace(r_lo, r_hi, points)
    values = np.array([func(float(r)) for r in grid], dtype=float)
    index = int(np.argmax(values))
    return float(grid[index]), float(values[index])


def format_row(layer: str, r_max: float, b_max_t: float, wpl_keV: float) -> str:
    return (
        f"{layer:<26} "
        f"r_Bmax = {r_max:8.5f} R_sun   "
        f"B_max = {b_max_t:10.4f} T   "
        f"omega_pl = {wpl_keV:9.6f} keV   "
        f"m_res = {wpl_keV * 1.0e3:9.3f} eV"
    )
