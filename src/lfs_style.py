"""
lfs_style — stile matplotlib riutilizzabile di Lorenzo.

Basato sulle Fig. 2, 3 e 7 del paper sul protein self-entanglement
(J. Chem. Phys. 155, 115101, 2021):
  - font: DejaVu Serif (Fig. 2/3)
  - griglia grigia tratteggiata, sempre attiva
  - coppia dicotomica before/after: blu puro / rosso puro (Fig. 2)
  - mappa di colore per superfici di energia libera: jet_r (Fig. 7)

Uso:
    import lfs_style
    lfs_style.apply()                 # una volta, all'inizio dello script
    ax.bar(x - w/2, y0, w, color=lfs_style.BEFORE, label="Before optimization")
    ax.bar(x + w/2, y1, w, color=lfs_style.AFTER,  label="After optimization")

Le impostazioni stanno in lfs.mplstyle (stessa cartella): modifica li'.
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

STYLE_FILE = Path(__file__).with_name("lfs.mplstyle")

# Palette (modificabile)
BEFORE = "#0000FF"   # blu Fig. 2
AFTER = "#FF0000"    # rosso Fig. 2
PAIR = (BEFORE, AFTER)
FES_CMAP = "jet_r"   # Fig. 7


def apply(**overrides):
    """Applica lo stile globalmente. Gli override sono rcParams, es. apply(**{'font.size': 14})."""
    plt.style.use(STYLE_FILE)
    if overrides:
        mpl.rcParams.update(overrides)


def context(**overrides):
    """Stile solo dentro un blocco:  with lfs_style.context(): ..."""
    return plt.style.context([STYLE_FILE, overrides] if overrides else STYLE_FILE)


def grid(ax, axis="both"):
    """Riapplica la griglia tratteggiata (utile dopo imshow/contourf, che la nascondono)."""
    ax.grid(True, axis=axis, linestyle="--", color="#b0b0b0", linewidth=0.6)
    ax.set_axisbelow(True)


if __name__ == "__main__":
    # Demo: ricrea Fig. 2 e una FES tipo Fig. 7
    import numpy as np
    apply()
    names = ["1APS", "1FKB", "1URN", "1TEN", "1BNZ_A", "1PSF",
             "2VIK", "1UBQ", "1HRC", "1BZP", "2CI2", "2ABD"]
    rng = np.random.default_rng(0)
    b = rng.uniform(0, .6, 12); a = np.clip(b + rng.uniform(.1, .5, 12), 0, .95)
    x = np.arange(12); w = .35
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.bar(x - w/2, b, w, color=BEFORE, label="Before optimization")
    ax1.bar(x + w/2, a, w, color=AFTER, label="After optimization")
    ax1.set_xticks(x, names, rotation=45); ax1.set_ylabel("Estimated folding frequency")
    ax1.legend(loc="upper left")
    X, Y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, 5, 200))
    G = 9 * (1 - np.exp(-((X - .35)**2 / .05 + (Y - 2)**2 / 2)))
    cf = ax2.contourf(X, Y, G, levels=50, cmap=FES_CMAP)
    fig.colorbar(cf, ax=ax2, label=r"Free Energy [$\varepsilon$]")
    ax2.set_xlabel("NC"); ax2.set_ylabel("RMSD")
    fig.tight_layout(); fig.savefig(Path(__file__).with_name("demo.png"))
    print("salvato demo.png")
