#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the project-wide colour contract.

A reader assumes consistent colour means the same thing, so every colour in
every figure must resolve from the one house palette. These tests are the
regression guard: they fail when a renderer names a colour of its own again,
which is how the palette silently fragmented before.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import re
from pathlib import Path

import cnsplots as cns
import matplotlib.colors as mcolors
import pytest

import figure_render
from figures import (
    DENSITY_CMAP,
    FURNITURE_COLOR,
    HOUSE_PALETTE,
    apply_house_style,
    observed_fitted_colors,
    series_colors,
)

# =============================================================================
# CONSTANTS
# =============================================================================
# Off-palette literals the renderers used to carry. Matched as whole words so a
# substring inside an identifier (e.g. "red" in "require_columns") is ignored.
FORBIDDEN_COLOR_LITERALS = (
    "red",
    "blue",
    "green",
    "black",
    "white",
    "gray",
    "grey",
    "darkred",
    "firebrick",
    "lightgray",
    "k",
    "b",
    "r",
    "g",
)

# Only these may name a colour: the registry is the single source, and viridis
# is the documented sequential override for density magnitude.
ALLOWED_COLOR_NAMES = frozenset({DENSITY_CMAP})

# Rec. 601 luma gap below which two colours collapse into each other in print.
MIN_PRINT_LUMA_GAP = 0.10

_QUOTED_LITERAL = re.compile(r"""(?:color|edgecolor|facecolor|c)\s*=\s*['"]([^'"]+)['"]""")
_HEX_LITERAL = re.compile(r"""['"](#[0-9A-Fa-f]{3,8})['"]""")


# =============================================================================
# HELPERS
# =============================================================================
def _renderer_sources() -> list[Path]:
    """Return every module file in the figure_render package."""
    package_dir = Path(figure_render.__file__).parent
    return sorted(package_dir.glob("*.py"))


def _luma(color: str) -> float:
    """Return the Rec. 601 luma of a colour, which is what greyscale printing keeps."""
    r, g, b = mcolors.to_rgb(color)
    return 0.299 * r + 0.587 * g + 0.114 * b


# =============================================================================
# TESTS
# =============================================================================
@pytest.mark.parametrize("source", _renderer_sources(), ids=lambda p: p.name)
def test_no_named_colour_literals_in_renderers(source: Path) -> None:
    """Assert no renderer assigns a colour by name instead of via the registry.

    A named literal bypasses HOUSE_PALETTE, so re-pointing the palette leaves
    that one mark behind at its old colour.
    """
    found = {
        literal
        for literal in _QUOTED_LITERAL.findall(source.read_text())
        if literal.lower() in FORBIDDEN_COLOR_LITERALS and literal not in ALLOWED_COLOR_NAMES
    }

    assert not found, f"{source.name} names colours directly: {sorted(found)}"


@pytest.mark.parametrize("source", _renderer_sources(), ids=lambda p: p.name)
def test_no_hex_colour_literals_in_renderers(source: Path) -> None:
    """Assert no renderer hardcodes a hex colour, even one copied from the palette.

    curves.py previously held '#c84c3a' and '#2f7e8f' — the palette's own first
    two entries — which drift out of sync the moment HOUSE_PALETTE changes.
    """
    found = set(_HEX_LITERAL.findall(source.read_text()))

    assert not found, f"{source.name} hardcodes hex colours: {sorted(found)}"


def test_house_palette_is_installed_globally() -> None:
    """Assert apply_house_style() sets the palette both cns.* plots and multipanel read."""
    apply_house_style()

    assert cns.settings.palette_qual == HOUSE_PALETTE


def test_registry_colours_come_from_the_house_palette() -> None:
    """Assert every data colour the registry hands out is a house-palette entry."""
    palette = {color.lower() for color in cns.get_hexcolors_from_apalette(range(10), HOUSE_PALETTE)}

    for color in (*observed_fitted_colors(), *series_colors(3)):
        assert color.lower() in palette, f"{color} is not a {HOUSE_PALETTE} entry"


def test_furniture_colour_is_outside_the_data_palette() -> None:
    """Assert reference lines and annotations cannot be mistaken for a data series."""
    data_colors = {
        color.lower()
        for color in (*observed_fitted_colors(), *series_colors(3))
    }

    assert FURNITURE_COLOR.lower() not in data_colors


def test_series_colours_survive_greyscale_printing() -> None:
    """Assert overlaid series stay distinguishable once colour is removed.

    Figures get printed and photocopied. The default palette order fails this:
    entries 0 and 1 differ by 0.026 luma, against 0.145 for the curated trio.
    """
    lumas = [_luma(color) for color in series_colors(3)]

    gaps = [
        abs(first - second)
        for index, first in enumerate(lumas)
        for second in lumas[index + 1 :]
    ]

    assert min(gaps) >= MIN_PRINT_LUMA_GAP, (
        f"Series colours collapse in greyscale: luma gaps {[round(g, 3) for g in gaps]}"
    )
