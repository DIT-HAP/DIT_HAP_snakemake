"""Gompertz curve model shared by curve fitting and its figures.

Input
-----
- None (pure numerical model).

Output
------
- None (returns arrays).

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np

# =============================================================================
# CONSTANTS
# =============================================================================
# np.exp overflows past ~709; clip well inside that bound.
EXPONENT_CLIP = 700

# =============================================================================
# CORE LOGIC
# =============================================================================
def sigmoid_function(x: np.ndarray, A: float, DR: float, DL: float) -> np.ndarray:
    """Calculate sigmoid function values using gompertz function."""
    if A == 0:
        return np.zeros_like(x)
    alpha = (DR * np.e) / A
    u = alpha * (DL - x) + 1
    exponent = np.clip(u, -EXPONENT_CLIP, EXPONENT_CLIP)
    return A * np.exp(-np.exp(exponent))
