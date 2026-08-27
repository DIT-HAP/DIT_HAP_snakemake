#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for per-point 2D kernel density estimation.

Asserts the grid-interpolated fast path agrees with exact KDE by rank, stays
deterministic, and degrades to None rather than raising on degenerate clouds.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-27
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import time

import numpy as np
import pytest
from scipy.stats import spearmanr

from figure_render._point_density import DENSITY_EXACT_MAX_N, point_density


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def two_clusters() -> tuple[np.ndarray, np.ndarray]:
    """A dense cluster at the origin plus a sparse halo, so density has a known ordering."""
    rng = np.random.default_rng(42)
    core = rng.normal(0, 0.3, (2, 2000))
    halo = rng.normal(0, 4.0, (2, 200))

    return np.concatenate([core[0], halo[0]]), np.concatenate([core[1], halo[1]])


# =============================================================================
# TESTS
# =============================================================================
def test_dense_region_scores_higher_than_sparse(two_clusters: tuple[np.ndarray, np.ndarray]) -> None:
    """Assert points near the cluster centre get higher density than distant ones."""
    x, y = two_clusters
    density = point_density(x, y)

    assert density is not None
    radius = np.hypot(x, y)
    assert density[radius < 0.5].mean() > density[radius > 5].mean()


def test_output_length_matches_input_with_nan_at_nonfinite() -> None:
    """Assert non-finite rows come back as NaN so the result indexes against the source frame."""
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1, 100)
    y = rng.normal(0, 1, 100)
    x[3] = np.nan
    y[7] = np.inf

    density = point_density(x, y)

    assert density is not None
    assert len(density) == 100
    assert np.isnan(density[[3, 7]]).all()
    assert np.isfinite(np.delete(density, [3, 7])).all()


def test_grid_path_agrees_with_exact_by_rank() -> None:
    """Assert the interpolated fast path preserves the density ordering the exact path gives.

    Absolute values differ slightly (different bandwidth from the subsample),
    so rank correlation is the meaningful comparison: colour mapping only needs
    the ordering to survive.
    """
    rng = np.random.default_rng(7)
    n = DENSITY_EXACT_MAX_N * 2
    x = rng.normal(0, 1, n)
    y = x * 0.6 + rng.normal(0, 0.8, n)

    grid = point_density(x, y)
    exact = point_density(x, y, exact_max_n=n)

    assert grid is not None and exact is not None
    correlation = spearmanr(grid, exact).statistic
    assert correlation > 0.9, f"Fast path diverged from exact KDE: rho={correlation:.3f}"


def test_grid_path_is_deterministic() -> None:
    """Assert the fixed seed makes repeat calls identical, so PNG baselines stay byte-exact."""
    rng = np.random.default_rng(11)
    n = DENSITY_EXACT_MAX_N * 2
    x = rng.normal(0, 1, n)
    y = rng.normal(0, 1, n)

    first = point_density(x, y)
    second = point_density(x, y)

    assert first is not None and second is not None
    np.testing.assert_array_equal(first, second)


@pytest.mark.parametrize(
    ("x", "y", "reason"),
    [
        (np.array([1.0, 2.0]), np.array([1.0, 2.0]), "fewer than 3 finite pairs"),
        (np.full(50, 3.0), np.arange(50.0), "zero range in x"),
        (np.arange(50.0), np.full(50, 3.0), "zero range in y"),
        (np.full(50, np.nan), np.arange(50.0), "all-NaN column"),
    ],
)
def test_degenerate_input_returns_none(x: np.ndarray, y: np.ndarray, reason: str) -> None:
    """Assert degenerate clouds return None instead of raising, so callers can fall back."""
    assert point_density(x, y) is None, f"Expected None for {reason}"


def test_mismatched_shapes_raise() -> None:
    """Assert misaligned inputs fail loudly rather than silently truncating."""
    with pytest.raises(ValueError, match="same shape"):
        point_density(np.arange(10.0), np.arange(5.0))


def test_large_input_completes_quickly() -> None:
    """Assert 50k points stay fast, so a regression to the O(n^2) path fails here."""
    rng = np.random.default_rng(3)
    x = rng.normal(0, 1, 50_000)
    y = rng.normal(0, 1, 50_000)

    start = time.perf_counter()
    density = point_density(x, y)
    elapsed = time.perf_counter() - start

    assert density is not None
    # Exact KDE on 50k points takes minutes; the grid path is well under a
    # second. A generous ceiling keeps this from flaking on a loaded machine.
    assert elapsed < 10, f"Density took {elapsed:.1f}s; the fast path may have been bypassed"
