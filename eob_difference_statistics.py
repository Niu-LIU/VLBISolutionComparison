#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: eob_diff_stat.py
"""
Created on Thu Dec 14 10:18:49 2017

@author: Neo(liuniu@smail.nju.edu.cn)

History:
N. Liu, 30 Jan 2018 : change the unit of Xp and Yp in function
                      "eop_diff_stats", from mas to uas.

"""

import numpy as np
from os import path
# from scipy import stats
from residual_plot import errorbarplot_res
from residual_statistics import stats_calc

__all__ = {"get_res_eop", "get_res_nut",
           "eop_diff_stats", "nut_diff_stats"}


# -----------------------------  FUNCTIONS -----------------------------
def get_res_eop(dif_file):
    '''Get the EOP residuals.
    '''

    # dif_file = "result/test.eob_diff"
    ep_mjd, dX, dXerr, dY, dYerr, dU, dUerr = np.genfromtxt(
        dif_file, unpack=True)

    return ep_mjd, dX, dXerr, dY, dYerr, dU, dUerr


def get_res_nut(dif_file):
    '''Get the Nutation offset residuals.
    '''

    ep_mjd, dX, dXerr, dY, dYerr = np.genfromtxt(
        dif_file, unpack=True)

    return ep_mjd, dX, dXerr, dY, dYerr


def get_res_eop(dif_file):
    '''Get the EOP residuals.
    '''

    # dif_file = "result/test.eob_diff"
    ep_mjd, dX, dXerr, dY, dYerr, dU, dUerr = np.genfromtxt(
        dif_file, usecols=np.arange(7), unpack=True)

    return ep_mjd, dX, dXerr, dY, dYerr, dU, dUerr


def get_res_nut(dif_file):
    '''Get the Nutation residuals.
    '''

    ep_mjd, dX, dXerr, dY, dYerr = np.genfromtxt(
        dif_file, usecols=np.arange(5), unpack=True)

    return ep_mjd, dX, dXerr, dY, dYerr


def eop_diff_stats(DIF_file):
    '''
    '''

    data_dir, data_fil = path.split(DIF_file)

    # Load data
    ep_mjd, dX, dXerr, dY, dYerr, dU, dUerr = get_res_eop(DIF_file)
    epo = (ep_mjd - 51544.5) / 365.25 + 2000.0

    # Plot
    # plot_res(epo, dX, "Xp", "mas", data_dir)
    # plot_res(epo, dY, "Yp", "mas", data_dir)
    # plot_res(epo, dU, "UT1-UTC", "ms", data_dir)

    errorbarplot_res(epo, dX, dXerr, "Xp", "mas", data_dir, label='Xp')
    errorbarplot_res(epo, dY, dYerr, "Yp", "mas", data_dir, label='Yp')
    errorbarplot_res(epo, dU, dUerr, "UT1-UTC", "ms", data_dir,
                     label='UT1')

    # Statistical calculation
    LOG_file = data_fil.replace(".eop_c04_dif", "_eop_c04_dif.log")
    flog = open("%s/logs/%s" % (data_dir, LOG_file), "w")

    # print("X pole (mas):", file=flog)
    # stats_calc(epo, dX, dXerr, flog)

    # print("Y pole: (mas)", file=flog)
    # stats_calc(epo, dY, dYerr, flog)

# Modified in 30 Jan 2018
    print("X pole (uas):", file=flog)
    stats_calc(epo, dX * 1000, dXerr * 1000, flog)

    print("Y pole (uas):", file=flog)
    stats_calc(epo, dY * 1000, dYerr * 1000, flog)

    print("UT1-UTC (us):", file=flog)
    stats_calc(epo, dU * 1000, dUerr * 1000, flog)

    flog.close()


def nut_diff_stats(DIF_file):
    '''
    '''

    # Load data
    ep_mjd, dX, dXerr, dY, dYerr = get_res_nut(DIF_file)
    epo = (ep_mjd - 51544.5) / 365.25 + 2000.0

    data_dir, data_fil = path.split(DIF_file)

    # Plot
    # plot_res(epo, dX, "dX", "mas", data_dir)
    # plot_res(epo, dY, "dY", "mas", data_dir)

    errorbarplot_res(epo, dX, dXerr, "dX", "mas", data_dir, label="dX")
    errorbarplot_res(epo, dY, dYerr, "dX", "mas", data_dir, label="dY")

    # Statistical calculation
    LOG_file = data_fil.replace(".nut_c04_dif", "_nut_c04_dif.log")
    flog = open("%s/logs/%s" % (data_dir, LOG_file), "w")

    # print("dX (mas):", file=flog)
    # stats_calc(epo, dX, dXerr, flog)

    # print("dY (mas):", file=flog)
    # stats_calc(epo, dY, dYerr, flog)

# Modified in 30 Jan 2018
    print("dX (uas):", file=flog)
    stats_calc(epo, dX * 1000, dXerr * 1000, flog)

    print("dY (uas):", file=flog)
    stats_calc(epo, dY * 1000, dYerr * 1000, flog)

    flog.close()


# --------------------------------- END --------------------------------
