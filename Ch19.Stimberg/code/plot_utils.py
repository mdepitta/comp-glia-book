"""
plot_utils.py

Module containing several methods related to plotting, formatting axes, etc.

Methods:
- formataxis(ax,**kwargs)
  General routine for basic axis formatting for publication-quality figures.
- rgbtuple(rgb)
  Return a tuple with normalized rgb values.
 
Maurizio De Pitta', The University of Chicago, Feb 26th, 2015.
"""

import matplotlib.pyplot as plt
import numpy as np
import os, sys
# sys.path.append(os.path.join(os.path.expanduser('~'),'Dropbox/Ongoing.Projects/pycustommodules'))
# import general_utils as gu

def hexc(rgb):
    """
    Simple function to format RGB-to-HEX colors appropriately.

    Input arguments:
    - rgb  : RGB tuple

    Return:
    - color string in HEX format
    """

    return '#%02x%02x%02x' % rgb

def adjust_spines(ax, spines, position=5, smart_bounds=False):
    """
    Set custom visibility and position of axes

    ax       : Axes
     Axes handle
    spines   : List
     String list of 'left', 'bottom', 'right', 'top' spines to show
    position : Integer
     Number of points for position of axis
    """
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', position))  # outward by 10 points
            spine.set_smart_bounds(smart_bounds)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in spines:
        ax.yaxis.set_ticks_position('right')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])
        ax.tick_params(axis='y',which='both',left='off',right='off')

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    elif 'top' in spines:
        ax.xaxis.set_ticks_position('top')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
        ax.tick_params(axis='x',which='both',bottom='off',top='off')

def adjust_ylabels(ax,x_offset=0):
    '''
    Scan all ax list and identify the outmost y-axis position. Setting all the labels to that position + x_offset.
    '''

    xc = 0.0
    for a in ax:
        xc = min(xc, (a.yaxis.get_label()).get_position()[0])

    for a in ax:
        a.yaxis.set_label_coords(xc + x_offset, (a.yaxis.get_label()).get_position()[1])


def make_patch_spines_invisible(ax):
    '''
    Used in the creation of a twin axis. Activate the twin frame but make the patch and spines invisible

    ax : twin axis handle
    '''
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)

def set_axlw(ax, lw=1.0):
    '''
    Adjust axis line width
    '''
    for axis in ax.spines.keys():
        ax.spines[axis].set_linewidth(lw)

def set_axfs(ax, fs=14):
    '''
    Adjust axis label font size
    '''
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)