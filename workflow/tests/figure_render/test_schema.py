#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for shared figure input schema validation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
import pytest

from figure_render._schema import require_columns


# =============================================================================
# TESTS
# =============================================================================
def test_passes_when_all_columns_present() -> None:
    """Assert a frame holding every required column validates silently."""
    df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})

    require_columns(df, ["a", "b"])


def test_passes_on_empty_frame_with_correct_columns() -> None:
    """Assert an empty frame with the right schema is accepted.

    Empty-input handling is exercised by every figure's empty-data test, which
    writes a header-only TSV; those must not fail validation.
    """
    df = pd.DataFrame(columns=["a", "b"])

    require_columns(df, ["a", "b"])


def test_raises_listing_every_missing_column() -> None:
    """Assert the error names all missing columns, not just the first."""
    df = pd.DataFrame({"a": [1]})

    with pytest.raises(ValueError) as excinfo:
        require_columns(df, ["a", "b", "c"])

    message = str(excinfo.value)
    assert "b" in message
    assert "c" in message


def test_error_message_includes_context() -> None:
    """Assert the caller-supplied context appears, so the failing file is identifiable."""
    df = pd.DataFrame({"a": [1]})

    with pytest.raises(ValueError, match="strand pairs"):
        require_columns(df, ["plus"], context="strand pairs TSV")


def test_preserves_required_order_in_message() -> None:
    """Assert missing columns are reported in the order the caller declared them."""
    df = pd.DataFrame({"z": [1]})

    with pytest.raises(ValueError) as excinfo:
        require_columns(df, ["a", "b"])

    message = str(excinfo.value)
    assert message.index("'a'") < message.index("'b'")
