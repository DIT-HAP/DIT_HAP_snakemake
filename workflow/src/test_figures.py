"""Tests for workflow/src/figures.py rendering utilities."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
import sys
from pathlib import Path

# Add src to path so we can import figures
sys.path.insert(0, str(Path(__file__).parent))


def test_apply_house_style_sets_arial_font():
    """apply_house_style() must set Arial to avoid Helvetica findfont warnings."""
    import cnsplots as cns
    from figures import apply_house_style

    apply_house_style()

    # Helvetica should be removed, Arial should be first
    assert "Helvetica" not in cns.settings.font_sans_serif
    assert cns.settings.font_sans_serif[0] == "Arial"
    assert cns.settings.panel_label_fontname == "Arial"


def test_apply_house_style_sets_pdf_fonttype_42():
    """PDF fonttype 42 ensures Illustrator can edit text."""
    import cnsplots as cns
    from figures import apply_house_style

    apply_house_style()
    assert cns.settings.pdf_fonttype == 42


def test_save_dual_creates_both_artifacts(tmp_path):
    """save_dual() must write both PDF and PNG."""
    import cnsplots as cns
    from figures import apply_house_style, save_dual

    apply_house_style()
    cns.figure(width=90, height=90)
    plt.plot([1, 2], [1, 2])

    stem = tmp_path / "test_fig"
    save_dual(stem)

    assert (tmp_path / "test_fig.pdf").exists()
    assert (tmp_path / "test_fig.review.png").exists()
    assert (tmp_path / "test_fig.pdf").stat().st_size > 1000
    assert (tmp_path / "test_fig.review.png").stat().st_size > 1000
