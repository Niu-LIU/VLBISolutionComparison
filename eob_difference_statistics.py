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
import sys
# from scipy import stats
from residual_plot import errorbarplot_res
from residual_statistics import stats_calc

__all__ = {"get_res_eop", "get_res_nut",
           "eop_diff_stats", "nut_diff_stats"}


# -----------------------------  FUNCTIONS -----------------------------
def eop_cut(epo, Xp, Xperr, Yp, Yperr, UT, UTerr,
            tb=None, te=None):
    '''Cut the EOP series among [tb, te]

    Paramters
    ---------
    epo : array of float
        epoch in Julian year
    Xp/Yp : array of float
        X/Y component of polar motion in mas
    UT : array of float
        UT1-TAI in ms
    Xperr/Yperr : array of float
        formal error of X/Y component of polar motion in mas
    UT : array of float
        formal error of UT1-TAI in ms
    tb/te : float
        beginning/end of the time span in Julian year.

    Returns
    ---------
    epo_ts : array of float
        epoch in Julian year
    Xp_ts/Yp_ts : array of float
        X/Y component of polar motion in mas
    UT_ts : array of float
        UT1-TAI in ms
    Xperr_ts/Yperr_ts : array of float
        formal error of X/Y component of polar motion in mas
    UTerr_ts : array of float
        formal error of UT1-TAI in ms
    '''

    if tb is None:
        if te is None:
            print("Please give at least one parameters!")
            sys.exit()
        else:
            con = (epo <= te)
    elif te is None:
        con = (epo >= tb)
    else:
        con = (epo <= te) & (epo >= tb)

    # Extract data
    epo_ts = epo[con]
    Xp_ts = Xp[con]
    Xperr_ts = Xperr[con]
    Yp_ts = Yp[con]
    Yperr_ts = Yperr[con]
    UT_ts = UT[con]
    UTerr_ts = UTerr[con]

    return [epo_ts, Xp_ts, Xperr_ts, Yp_ts, Yperr_ts, UT_ts, UTerr_ts]


def nut_cut(epo, dX, dXerr, dY, dYerr, tb=None, te=None):
    '''Cut the EOP series among [tb, te]

    Paramters
    ---------
    epo : array of float
        epoch in Julian year
    dX/dY : array of float
        X/Y component of polar motion in mas
    dXerr/dYerr : array of float
        formal error of X/Y component of polar motion in mas
    tb/te : float
        beginning/end of the time span in Julian year.

    Returns
    ---------
    epo_ts : array of float
        epoch in Julian year
    dX_ts/dY_ts : array of float
        X/Y component of Nutation offset wrt IAU model in mas
    dXerr_ts/dYerr_ts : array of float
        formal error of X/Y component of polar motion in mas
    '''

    if tb is None:
        if te is None:
            print("Please give at least one parameters!")
            sys.exit()
        else:
            con = (epo <= te)
    elif te is None:
        con = (epo >= tb)
    else:
        con = (epo <= te) & (epo >= tb)

    # Extract data
    epo_ts = epo[con]
    dX_ts = dX[con]
    dXerr_ts = dXerr[con]
    dY_ts = dY[con]
    dYerr_ts = dYerr[con]

    return [epo_ts, dX_ts, dXerr_ts, dY_ts, dYerr_ts]


# def get_res_eop(dif_file):
#     '''Get the EOP residuals.
#     '''

#     # dif_file = "result/test.eob_diff"
#     ep_mjd, Xp, Xperr, Yp, Yperr, UT, UTerr = np.genfromtxt(
#         dif_file, unpack=True)

#     return ep_mjd, Xp, Xperr, Yp, Yperr, UT, UTerr


# def get_res_nut(dif_file):
#     '''Get the Nutation offset residuals.
#     '''

#     ep_mjd, Xp, Xperr, Yp, Yperr = np.genfromtxt(
#         dif_file, unpack=True)

#     return ep_mjd, Xp, Xperr, Yp, Yperr
#
#
def get_res_eop(dif_file):
    '''Get the EOP residuals.
    '''

    # dif_file = "result/test.eob_diff"
    ep_mjd, Xp, Xperr, Yp, Yperr, UT, UTerr = np.genfromtxt(
        dif_file, usecols=np.arange(7), unpack=True)

    return ep_mjd, Xp, Xperr, Yp, Yperr, UT, UTerr


def get_res_nut(dif_file):
    '''Get the Nutation residuals.
    '''

    ep_mjd, Xp, Xperr, Yp, Yperr = np.genfromtxt(
        dif_file, usecols=np.arange(5), unpack=True)

    return ep_mjd, Xp, Xperr, Yp, Yperr


def eop_diff_stats(DIF_file):
    '''
    '''

    data_dir, data_fil = path.split(DIF_file)

    # Load data
    ep_mjd, Xp, Xperr, Yp, Yperr, UT, UTerr = get_res_eop(DIF_file)
    epo = (ep_mjd - 51544.5) / 365.25 + 2000.0

    # Plot
    # plot_res(epo, Xp, "Xp", "mas", data_dir)
    # plot_res(epo, Yp, "Yp", "mas", data_dir)
    # plot_res(epo, UT, "UT1-UTC", "ms", data_dir)

    errorbarplot_res(epo, Xp, Xperr, "Xp", "mas", data_dir, label='Xp')
    errorbarplot_res(epo, Yp, Yperr, "Yp", "mas", data_dir, label='Yp')
    errorbarplot_res(epo, UT, UTerr, "UT1-UTC", "ms", data_dir,
                     label='UT1')

    # Statistical calculation
    LOG_file = data_fil.replace(".eop_c04_dif", "_eop_c04_dif.log")
    flog = open("%s/logs/%s" % (data_dir, LOG_file), "w")

    # print("X pole (mas):", file=flog)
    # stats_calc(epo, Xp, Xperr, flog)

    # print("Y pole: (mas)", file=flog)
    # stats_calc(epo, Yp, Yperr, flog)

# Modified in 30 Jan 2018
    print("X pole (uas):", file=flog)
    stats_calc(epo, Xp * 1000, Xperr * 1000, flog)

    print("Y pole (uas):", file=flog)
    stats_calc(epo, Yp * 1000, Yperr * 1000, flog)

    print("UT1-UTC (us):", file=flog)
    stats_calc(epo, UT * 1000, UTerr * 1000, flog)

    flog.close()

# Post-1993.0 data
    tb = 1993.0
    [epo_ts,
     Xp_ts, Xperr_ts, Yp_ts, Yperr_ts, UT_ts, UTerr_ts] = eop_cut(
        epo, Xp, Xperr, Yp, Yperr, UT, UTerr, tb)
    # Statistical calculation
    LOG_file = data_fil.replace(".eop_c04_dif",
                                "_eop_c04_post%.0f_dif.log" % tb)
    flog = open("%s/logs/%s" % (data_dir, LOG_file), "w")
    print("Use post %.1f data" % tb, file=flog)
    print("X pole (uas):", file=flog)
    stats_calc(epo_ts, Xp_ts * 1000, Xperr_ts * 1000, flog)

    print("Y pole (uas):", file=flog)
    stats_calc(epo_ts, Yp_ts * 1000, Yperr_ts * 1000, flog)

    print("UT1-UTC (us):", file=flog)
    stats_calc(epo_ts, UT_ts * 1000, UTerr_ts * 1000, flog)

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

    errorbarplot_res(epo, dX, dXerr, "X_p", "mas", data_dir, label="dX")
    errorbarplot_res(epo, dY, dYerr, "Y_p", "mas", data_dir, label="dY")

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

# Post-1993.0 data
    tb = 1993.0
    [epo_ts, dX_ts, dXerr_ts, dY_ts, dYerr_ts] = nut_cut(
        epo, dX, dXerr, dY, dYerr, tb)
    # Statistical calculation
    LOG_file = data_fil.replace(".nut_c04_dif",
                                "_nut_c04_post%.0f_dif.log" % tb)
    flog = open("%s/logs/%s" % (data_dir, LOG_file), "w")
    print("Use post %.1f data" % tb, file=flog)
    print("dX (uas):", file=flog)
    stats_calc(epo_ts, dX_ts * 1000, dXerr_ts * 1000, flog)

    print("dY (uas):", file=flog)
    stats_calc(epo_ts, dY_ts * 1000, dYerr_ts * 1000, flog)

    flog.close()


# --------------------------------- END --------------------------------
