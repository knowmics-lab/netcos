"""
Figure style of a case study, from the FIGURES block of its case_conf.

FIGURES['style'] = 'lfs' (default) applies lfs_style (lfs.mplstyle);
'matplotlib' keeps the matplotlib defaults. FIGURES['rc'] (rcParams) is
applied on top in both cases. FIGURES['colours'] = (control colour, case
colour), default (lfs_style.BEFORE, lfs_style.AFTER).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import lfs_style


def figures(cc):
    """The FIGURES block of the study (required: it holds the figure text)."""
    if not hasattr(cc, 'FIGURES'):
        raise ValueError(f'case_conf of {cc.P.slug}: FIGURES is required (figure text of every step)')
    return cc.FIGURES


def apply_style(cc):
    f = figures(cc)
    if f.get('style', 'lfs') == 'lfs':
        lfs_style.apply(**f.get('rc', {}))
    else:
        plt.rcParams.update(f.get('rc', {}))


def colours(cc):
    """(BLUE, RED): control / case colour of the study figures."""
    return figures(cc).get('colours', (lfs_style.BEFORE, lfs_style.AFTER))
