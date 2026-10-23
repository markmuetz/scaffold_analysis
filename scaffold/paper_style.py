"""Shared matplotlib style for the rce_org1 paper figures.

Import and call `use_paper_style()` before plotting. Figures are drawn at final print size
(single column ~8 cm, full width ~17 cm) and saved as vector PDF with `save_paper_fig`.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

SINGLE_COL_CM = 8
FULL_WIDTH_CM = 17

# Lifetime-distribution colours (colour-blind safe).
LIFETIME_COLOURS = {'all': '#4d4d4d', 'simple': '#2f7fc1', 'complex': '#d1603d'}

PAPER_RC = {
    'font.family': 'sans-serif',
    'font.size': 7.5,
    'axes.labelsize': 7.5,
    'axes.titlesize': 7.5,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.linewidth': 0.6,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'lines.linewidth': 1.2,
    'lines.markersize': 4,
    'legend.frameon': False,
    'savefig.bbox': 'tight',
    'pdf.fonttype': 42,
    'figure.autolayout': False,
}


def use_paper_style():
    mpl.rcParams.update(PAPER_RC)


def cm_to_inch(*vals):
    return [v / 2.54 for v in vals]


def panel_label(ax, label, **kwargs):
    """Panel label, e.g. "(a)", at the top-left of an axes."""
    kwargs.setdefault('fontweight', 'bold')
    ax.text(0.01, 0.99, label, transform=ax.transAxes, ha='left', va='top', **kwargs)


def panel_labels(axes, **kwargs):
    for i, ax in enumerate(axes):
        panel_label(ax, '({})'.format('abcdefghijklmnopqrstuvwxyz'[i]), **kwargs)


def save_paper_fig(filename, fig=None):
    """Save as vector PDF and as a PNG preview next to it."""
    fig = fig or plt.gcf()
    stem = filename[:-4] if filename.endswith(('.pdf', '.png')) else filename
    fig.savefig(stem + '.pdf', bbox_inches='tight')
    fig.savefig(stem + '.png', bbox_inches='tight', dpi=200)
