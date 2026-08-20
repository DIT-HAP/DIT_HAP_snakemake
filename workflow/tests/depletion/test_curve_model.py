"""Tests for the shared Gompertz curve model."""

import numpy as np
import pytest

from depletion.curve_model import sigmoid_function


def test_zero_amplitude_returns_zeros():
    """A == 0 must short-circuit to zeros rather than divide by zero."""
    x = np.array([0.0, 1.0, 2.0])
    assert np.array_equal(sigmoid_function(x, 0.0, 1.0, 1.0), np.zeros_like(x))


def test_monotonic_decreasing_in_x():
    """With positive amplitude the Gompertz curve decreases as x grows."""
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y = sigmoid_function(x, 2.0, 0.5, 2.0)
    assert np.all(np.diff(y) >= -1e-12)


def test_extreme_input_does_not_overflow():
    """Exponent clipping must keep a huge negative x finite."""
    y = sigmoid_function(np.array([-1e6]), 2.0, 0.5, 2.0)
    assert np.all(np.isfinite(y))


def test_matches_amplitude_at_large_x():
    """As x grows past DL the curve approaches the amplitude A."""
    y = sigmoid_function(np.array([1e3]), 2.0, 0.5, 2.0)
    assert y[0] == pytest.approx(2.0, abs=1e-9)
