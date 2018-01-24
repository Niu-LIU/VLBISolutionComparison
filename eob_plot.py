#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: eob_plot.py
"""
Created on Wed Jan 10 16:32:09 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import matplotlib.pyplot as plt
from read_eob import read_eob


# -----------------------------  FUNCTIONS -----------------------------
def errorbarplot_res(x, y, err, lab, unit):
    '''
    '''

    plt.figure(figsize=(10, 4))
    # plt.plot(x, y, '.', markersize=0.1)
    # plt.xlabel('MJD')
    # plt.title("%s(%s)" % (lab, unit))
    # plt.xlim([1979.0, 2018.0])

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.errorbar(x, y, yerr=err,
                fmt='.', ms=1,
                ecolor='grey',
                elinewidth=0.1)

    ax.set_xlabel('MJD')
    ax.set_title("%s(%s)" % (lab, unit))
    ax.set_xlim([1979.5, 2018.5])

    # plt.ylim([-30, 30])
    # plt.savefig("%s_residual30.eps" % lab)
    # ax.set_ylim([-10, 10])
    # plt.savefig("%s_residual10.eps" % lab)
    plt.ylim([-5, 5])
    plt.savefig("%s_residual05.eps" % lab)
    plt.ylim([-1, 1])
    plt.savefig("%s_residual01.eps" % lab)
    plt.close()


EOB_file = "/home/nliu/solutions/opa2018a/opa2018a.eob"

[dbname, obsnum, tag_eop, Xp, Xp_err, Yp, Yp_err, U, U_err,
 XR, XR_err, YR, YR_err, UR, UR_err,
 corXY, corXU, corYU, corXUR, corYUR, corUUR,
 tag_nut, dX, dX_err, dY, dY_err, cordXdY] = read_eob(EOB_file)

epo = (tag_nut - 51544.5) / 365.25 + 2000.0

# Plot
errorbarplot_res(epo, dX, dX_err, "dX", "mas")
errorbarplot_res(epo, dY, dY_err, "dY", "mas")

# --------------------------------- END --------------------------------
