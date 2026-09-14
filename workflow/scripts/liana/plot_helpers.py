"""
Dark-themed spatial scatter helper for the LIANA pipeline 03_spatial outputs.

Conventions:
  - Black figure + axis background
  - Zero/non-finite values plotted first in dim gray (background layer)
  - Nonzero values plotted on top, sorted low->high (brightest most visible)
  - vmax = quantile of nonzero values (default 99th) — avoids one outlier dominating
"""

import numpy as np
import matplotlib.pyplot as plt


def _to_dense_1d(x):
    """Coerce sparse/dense/scalar matrix slices to a 1d numpy array."""
    if hasattr(x, 'toarray'):
        x = x.toarray()
    return np.asarray(x).ravel()


def plot_spatial_dark(
    ax,
    x, y, values,
    *,
    title='',
    point_size=8,
    cmap='plasma',
    zero_color='#404040',
    zero_alpha=0.5,
    vmax_quantile=0.99,
    show_colorbar=True,
    title_fontsize=11,
):
    """
    Render a single spatial panel with dark theme.

    Parameters
    ----------
    ax : matplotlib Axes
    x, y : array-like (n,) cell coordinates
    values : array-like (n,) per-cell expression / inflow values
    title : str
    """
    x = np.asarray(x); y = np.asarray(y)
    v = _to_dense_1d(values)

    ax.set_facecolor('black')
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    finite = np.isfinite(v)
    nonzero = finite & (v > 0)

    # Background: zero / non-finite cells in dim gray
    bg_mask = ~nonzero
    if bg_mask.any():
        ax.scatter(x[bg_mask], y[bg_mask],
                   c=zero_color, s=point_size,
                   edgecolors='none', alpha=zero_alpha,
                   rasterized=True)

    # Foreground: nonzero cells, sorted low->high so high values draw on top
    if nonzero.any():
        order = np.argsort(v[nonzero])
        x_nz = x[nonzero][order]
        y_nz = y[nonzero][order]
        v_nz = v[nonzero][order]
        vmax = float(np.quantile(v_nz, vmax_quantile))
        if vmax <= 0:
            vmax = float(v_nz.max())  # fall back if quantile is 0

        sc = ax.scatter(x_nz, y_nz, c=v_nz, s=point_size,
                        cmap=cmap, vmin=0, vmax=vmax,
                        edgecolors='none', rasterized=True)
        if show_colorbar:
            cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(colors='white')
            cbar.outline.set_edgecolor('white')
            for label in cbar.ax.get_yticklabels():
                label.set_color('white')

    if title:
        ax.set_title(title, color='white', fontsize=title_fontsize)


def save_dark_figure(fig, path, dpi=150):
    """Save with black face color preserved (otherwise the bg goes white in the PDF)."""
    fig.savefig(path, dpi=dpi, bbox_inches='tight',
                facecolor='black', edgecolor='none', pad_inches=0.1)
    plt.close(fig)
