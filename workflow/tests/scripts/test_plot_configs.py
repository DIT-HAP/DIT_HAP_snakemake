"""Validation tests for the figure scripts' PlotConfig dataclasses."""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pytest

# =============================================================================
# TESTS
# =============================================================================
@pytest.mark.parametrize(
    ("stem", "kwargs"),
    [
        ("plot_ma_plot", {"basemean_path": "nope_baseMean.tsv", "lfc_path": "nope_LFC.tsv"}),
        ("plot_ma_plot_replicates", {"ma_values_path": "nope_ma_values.tsv"}),
        ("plot_dispersions", {"dispersion_data_path": "nope_dispersion.tsv"}),
        ("plot_distribution_of_curve_fitting", {"fitting_stats_path": "nope_stats.tsv"}),
        ("plot_curve_fitting", {"fitting_stats_path": "nope_stats.tsv", "lfc_path": "nope_LFC.tsv"}),
    ],
)
def test_config_rejects_missing_input(load_script, tmp_path: Path, stem: str, kwargs: dict) -> None:
    """Assert each script's PlotConfig rejects a non-existent input path."""
    module = load_script(stem)
    paths = {key: tmp_path / value for key, value in kwargs.items()}

    with pytest.raises(ValueError, match="does not exist"):
        module.PlotConfig(output_stem=tmp_path / "out", **paths)
