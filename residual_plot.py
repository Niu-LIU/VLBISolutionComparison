#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: residual_plot.py
"""
Created on Fri Jan 12 10:01:25 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import matplotlib.pyplot as plt

__all__ = {"plot_res", "errorbarplot_res", "errorbarplot_nut"}


# -----------------------------  FUNCTIONS -----------------------------
def plot_res(x, y, lab, unit, plot_dir):
    '''
    '''

    plt.figure(figsize=(10, 4))
    # plt.plot(x, y, '.', markersize=0.1)
    # plt.xlabel('MJD')
    # plt.title("%s(%s)" % (lab, unit))
    # plt.xlim([1979.0, 2018.0])

    plt.plot(x, y, '.', markersize=0.1)
    plt.xlabel('MJD (year)',
               fontsize=15)
    plt.title("%s(%s)" % (lab, unit),
              fontsize=18)
    plt.xlim([1979.5, 2018.5])

    # plt.ylim([-30, 30])
    # plt.savefig("%s/plots/%s_residual30.eps" % (plot_dir, lab))
    # plt.ylim([-10, 10])
    # plt.savefig("%s/plots/%s_residual10.eps" % (plot_dir, lab))
    plt.ylim([-5, 5])
    plt.savefig("%s/plots/%s_dif05.eps" % (plot_dir, lab))
    # plt.ylim([-1, 1])
    # plt.savefig("%s/plots/%s_residual01.eps" % (plot_dir, lab))
    plt.close()


def errorbarplot_res(x, y, err, comp, unit, plot_dir, label='',
                     slope=0, intercept=0):
    '''
    '''

    plt.figure(figsize=(10, 4))
    # plt.plot(x, y, '.', markersize=0.1)
    # plt.xlabel('MJD')
    # plt.title("%s(%s)" % (comp, unit))
    # plt.xlim([1979.0, 2018.0])

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.errorbar(x, y, yerr=err,
                fmt='b.', ms=0.5,
                ecolor='grey',
                elinewidth=0.1,
                clip_on=True)
    if slope:
        ax.plot(x, slope * (x - 2005.) + intercept,
                'r-', ms=1, clip_on=False)

    ax.set_xlabel('MJD (year)',
                  fontsize=15)
    ax.set_title("%s(%s)" % (comp, unit),
                 fontsize=18)
    ax.set_xlim([1979.5, 2018.5])

    ax.set_ylim([-50, 50])
    plt.savefig("%s/plots/%s_dif50.eps" % (plot_dir, label))
    ax.set_ylim([-30, 30])
    plt.savefig("%s/plots/%s_dif30.eps" % (plot_dir, label))
    ax.set_ylim([-10, 10])
    plt.savefig("%s/plots/%s_dif10.eps" % (plot_dir, label))
    plt.ylim([-5, 5])
    plt.savefig("%s/plots/%s_dif05.eps" % (plot_dir, label))
    # plt.ylim([-1, 1])
    # plt.savefig("%s/plots/%s_dif01.eps" % (plot_dir, label))
    plt.close()


def errorbarplot_nut(epo, dX, dXerr, dY, dYerr, plot_dir):
    '''
    '''

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, figsize=(10, 4))

    ax0.errorbar(epo, dX, yerr=dXerr,
                 fmt='.', ms=1,
                 ecolor='grey',
                 elinewidth=0.1)

    ax1.errorbar(epo, dY, yerr=dYerr,
                 fmt='.', ms=1,
                 ecolor='grey',
                 elinewidth=0.1)

    ax0.set_ylabel("dX (mas)",
                   fontsize=15)
    ax0.set_ylim([-1, 1])

    ax1.set_ylabel("dY (mas)",
                   fontsize=15)
    ax1.set_ylim([-1, 1])

    ax1.set_xlabel('MJD (year)',
                   fontsize=15)
    ax1.set_xlim([1979.5, 2018.5])

    plt.suptitle("Nutation offset w.r.t IAU 2006/2000A model",

                 fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1, top=0.92)

    plt.savefig("%s/plots/nutation_offset.eps" % plot_dir)
    plt.close()
# --------------------------------- END --------------------------------
