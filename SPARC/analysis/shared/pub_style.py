"""
Publication-quality dark-theme styling constants for empirical audit figures.

These constants are used identically across all empirical audit scripts
(DESI, SPARC, Pantheon, Planck, SLACS, Spin-Parity, Disulfide, Lattice, Geometry).
Centralised here to eliminate duplication and ensure visual consistency.

The values match the PUB_* constants that have been the standard across
the established audit figure style used in this workspace.
"""

from __future__ import annotations


# ── Publication rendering ───────────────────────────────────────────
PUB_DPI = 300

# ── Background / foreground ─────────────────────────────────────────
PUB_BG         = "#0a0a0f"
PUB_FG         = "#e0e0e0"
PUB_GRID       = "#1a1a2e"
PUB_SPINE      = "#333344"
PUB_FACE       = "#0e0e16"

# ── Font sizes ──────────────────────────────────────────────────────
PUB_TITLE_SIZE  = 13
PUB_LABEL_SIZE  = 11
PUB_TICK_SIZE   = 10
PUB_LEGEND_SIZE = 10
PUB_ANNOT_SIZE  = 9

# ── Accent colours ──────────────────────────────────────────────────
C_CYAN  = "#4cc9f0"
C_PINK  = "#f72585"
C_GREEN = "#7ef08a"
C_GOLD  = "#f6c453"
C_CORAL = "#ff6b6b"
C_LAV   = "#b57edc"


def style_ax(ax):
    """Apply the PUB_* dark theme to a single matplotlib Axes object.

    Sets the face colour, tick colours/sizes, spine colours, and grid.
    This is the standard styling applied to every audit panel.
    """
    ax.set_facecolor(PUB_FACE)
    ax.tick_params(colors=PUB_FG, labelsize=PUB_TICK_SIZE)
    for spine in ax.spines.values():
        spine.set_color(PUB_SPINE)
    ax.grid(True, color=PUB_GRID, alpha=0.5, which="both")
