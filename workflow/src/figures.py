"""
Shared Rendering Library
========================

House style configuration and dual-artifact saving for cnsplots figures.
This is a library module (no main(), no CLI) per python-script-conventions.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-13
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import cnsplots as cns


# =============================================================================
# CONSTANTS
# =============================================================================
# Journal two-column width: 180 mm ≈ 510 px at 72 DPI base
JOURNAL_WIDTH_PX = 510
JOURNAL_HEIGHT_PX = 425  # ~150 mm, 4:3 aspect


# =============================================================================
# CORE LOGIC
# =============================================================================
def apply_house_style() -> None:
    """Configure cnsplots with project house style: Ecotyper1 palette, Arial font, PDF fonttype 42."""
    # Remove Helvetica to prevent findfont warnings; Arial is installed
    cns.settings.font_sans_serif = ("Arial", "DejaVu Sans")
    cns.settings.panel_label_fontname = "Arial"

    # Fonttype 42 embeds TrueType, editable in Illustrator
    cns.settings.pdf_fonttype = 42

    # Ecotyper1: colour-blind-safe, print-validated
    cns.setup_matplotlib(color_cycle="Ecotyper1")


def save_dual(stem: Path | str) -> None:
    """Save the current figure as a journal PDF plus a review PNG, given a stem with no extension."""
    # Append rather than replace suffixes: Path.with_suffix()/.stem would truncate
    # at the first dot, so a stem like "HD1328-4.YES0_corr" would silently collapse
    # to "HD1328-4" and two samples would overwrite each other's figures.
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    cns.savefig(stem.parent / f"{stem.name}.pdf")
    cns.savefig(stem.parent / f"{stem.name}.review.png")
