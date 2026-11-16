"""
Figure style and figure text of a case study.

Every figure uses the house style lfs_style (lfs.mplstyle: DejaVu Serif,
grey dashed grid, pure blue / pure red pair). The FIGURES block of the
case_conf holds only figure text (dataset, disease, axis and panel labels).

CONFIG_COLOURS gives one colour per CS configuration, the same in every
figure: the DEG family in blues, the MITHrIL node family in reds, MITHrIL
pathways in black.
"""

import matplotlib
matplotlib.use('Agg')

import lfs_style

BLUE, RED = lfs_style.BEFORE, lfs_style.AFTER

CONFIG_COLOURS = {
    'DEG': BLUE,
    'DEG_LM': '#00A0C8',
    'mith': RED,
    'mith_LM': '#FF8C00',
    'mith_pw': '#8B00C8',
}


def figures(cc):
    """The FIGURES block of the study (required: it holds the figure text)."""
    if not hasattr(cc, 'FIGURES'):
        raise ValueError(f'case_conf of {cc.P.slug}: FIGURES is required (figure text of every step)')
    return cc.FIGURES


def apply_style():
    lfs_style.apply()


def colours():
    """(BLUE, RED): control / case colour of the study figures."""
    return BLUE, RED


def config_colour(config):
    return CONFIG_COLOURS.get(config, '0.5')
