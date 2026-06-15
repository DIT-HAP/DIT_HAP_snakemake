# ================================ Imports =================================
import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib import pyplot as plt

# ================================ Constants =================================
SCRIPT_DIR = Path(__file__).parent.resolve()
TARGET_path = str((SCRIPT_DIR / "../../config/DIT_HAP.mplstyle").resolve())
plt.style.use(TARGET_path)
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# ================================ Functions =================================
def create_scatter_correlation_plot(
    x: pd.Series | np.ndarray | list,
    y: pd.Series | np.ndarray | list,
    ax: Axes,
    xscale: None | str = None,
    yscale: None | str = None,
    show_diagonal: bool = True,
    **kwargs
) -> Axes:
    """Create correlation plot for a single file with statistics."""    
    # Plot data points

    x, y = np.array(x), np.array(y)
    x, y = x[~np.isnan(x) & ~np.isnan(y)], y[~np.isnan(x) & ~np.isnan(y)]

    mask = np.isfinite(x) & np.isfinite(y)
    if xscale == 'log':
        mask &= (x > 0)
    if yscale == 'log':
        mask &= (y > 0)

    x = x[mask]
    y = y[mask]
    
    ax.scatter(
        x, y,
        alpha=0.5,
        s=10,
        facecolor="none",
        edgecolor="gray",
        rasterized=True,
        **kwargs
    )

    # Get axis limits for diagonal line
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Plot diagonal reference line (y=x)
    if show_diagonal:
        min_val = min(min(xlim), min(ylim))
        max_val = max(max(xlim), max(ylim))
        ax.plot([min_val, max_val], [min_val, max_val], 
                'k--', alpha=0.8, linewidth=2)
    
    # Set log scale for both axes
    if xscale == 'log':
        ax.set_xscale('log')
        x_for_fitting = np.log10(x)
    else:
        x_for_fitting = x

    if yscale == 'log':
        ax.set_yscale('log')
        y_for_fitting = np.log10(y)
    else:
        y_for_fitting = y

    # Calculate correlation statistics
    # Pearson correlation coefficient
    pcc = np.corrcoef(x_for_fitting, y_for_fitting)[0, 1]
    # R-squared
    r_squared = pcc**2
    try:
        # Linear regression
        slope, intercept = np.polyfit(x_for_fitting, y_for_fitting, 1)
        # RMSE
        y_pred = intercept + slope * x_for_fitting
        rmse = np.sqrt(np.mean((y_for_fitting - y_pred)**2))
    except Exception as e:
        print("Error in linear regression:", e)
        slope, intercept, rmse = np.nan, np.nan, np.nan
    
    # Add statistics text box
    stats_text = []
    stats_text.append(f"Data points: {len(x):,}")
    stats_text.append(f"PCC: {pcc:.4f}")
    stats_text.append(f"R²: {r_squared:.4f}")
    stats_text.append(f"Slope: {slope:.4f}")
    stats_text.append(f"Intercept: {intercept:.4f}")
    stats_text.append(f"RMSE: {rmse:.4f}")

    # Add text box with statistics
    textstr = '\n'.join(stats_text)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, verticalalignment='top')
    
    return ax

def donut_chart(
    values: list[int],
    labels: list[str],
    colors: list[str],
    center_text: str = "",
    ax: Axes | None = None,
) -> Axes | Figure:
    """Create a donut chart with given values, labels, and colors."""
    # Create axis if not provided
    return_ax = True
    if ax is None:
        # Create the donut chart
        fig, ax = plt.subplots()
        return_ax = False
    
    # Create the donut chart
    ax.pie(
        values, 
        # labels=labels,
        colors=colors,
        autopct=lambda pct: f'{pct:.1f}%\n({int(round(pct/100*sum(values))):,})',
        startangle=90,
        pctdistance=0.75,
        wedgeprops=dict(width=0.5, edgecolor='white'),
        textprops={'fontsize': 22, 'weight': 'bold'},
    )

    # Add total count in center
    ax.text(
        0, 0, 
        center_text, 
        ha='center', 
        va='center', 
        fontsize=26, 
        fontweight='bold'
    )

    # Equal aspect ratio ensures that pie is drawn as a circle
    ax.axis('equal')

    if return_ax:
        return ax
    else:
        return fig
