#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: check_nutation_offset.py
"""
Created on Thu Jan 11 12:37:40 2018

@author: Neo(liuniu@smail.nju.edu.cn)

History
N.Liu, 31JAN18 : change the unit of Nutation offset in function stats_calc,
                 from mas to uas
"""

import numpy as np
from os import path
from scipy import stats
import matplotlib.pyplot as plt
from linear_regression import linear_regression
from read_nut import read_nut
from read_eob import read_eob
from residual_statistics import stats_calc
from residual_plot import errorbarplot_nut

__all__ = {"nutation_offset_stat"}


# -----------------------------  FUNCTIONS -----------------------------
# def nutation_offset_stat(NUT_file):
def nutation_offset_stat(EOB_file):
    '''
    '''

    # Load data
    # dbname, epo, obsnum, dX, dXerr, dY, dYerr = read_nut(NUT_file)
    [dbname, obsnum, tag_eop, Xp, Xp_err, Yp, Yp_err, U, U_err,
     XR, XR_err, YR, YR_err, UR, UR_err,
     corXY, corXU, corYU, corXUR, corYUR, corUUR,
     tag_nut, dX, dX_err, dY, dY_err, cordXdY] = read_eob(EOB_file)
    epo = (tag_nut - 51544.5) / 365.25 + 2000.0

    data_dir, data_fil = path.split(EOB_file)

    # Plot
    # plot_res(epo, dX, "dX", "mas", data_dir)
    # plot_res(epo, dY, "dY", "mas", data_dir)

    errorbarplot_nut(epo, dX, dX_err, dY, dY_err, data_dir)

    # Statistical calculation
    # LOG_file = data_fil.replace(".nut", "_nut.log")
    LOG_file = data_fil.replace(".eob", "_nut.log")
    flog = open("%s/logs/%s" % (data_dir, LOG_file), "w")

    # units are changed on 31 JAN 2018.
    dX, dX_err = dX * 1.e3, dX_err * 1.e3
    dY, dY_err = dY * 1.e3, dY_err * 1.e3

    # Statistics for all-time series
    print("For all-time:", file=flog)
    # units are changed on 31 JAN 2018.
    # print("dX (mas):", file=flog)
    print("dX (uas):", file=flog)
    stats_calc(epo, dX, dX_err, flog)

    # units are changed on 31 JAN 2018.
    # print("dY: (mas)", file=flog)
    print("dY: (uas)", file=flog)
    stats_calc(epo, dY, dY_err, flog)

    # Statistics for post-1995.0 series
    splitpoint = 1993.
    print("For post-%.1f:" % splitpoint, file=flog)

    con = (epo >= splitpoint)
    epo = epo[con]
    dX = dX[con]
    dX_err = dX_err[con]
    dY = dY[con]
    dY_err = dY_err[con]

# units are changed on 31 JAN 2018.
    # print("dX (mas):", file=flog)
    print("dX (uas):", file=flog)
    stats_calc(epo, dX, dX_err, flog)

# units are changed on 31 JAN 2018.
    # print("dY: (mas)", file=flog)
    print("dY: (uas)", file=flog)
    stats_calc(epo, dY, dY_err, flog)

    flog.close()


# NUT_file = "/home/nliu/solutions/GalacticAberration/opa2018a_ga/opa2018a_ga.eob"
# nutation_offset_stat(NUT_file)
# --------------------------------- END --------------------------------
