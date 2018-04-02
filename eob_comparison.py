#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: EOB_diff.py
"""
Created on Tue Jan  2 15:09:14 2018

@author: Neo(liuniu@smail.nju.edu.cn)


Calculate the CPO position (dX, dY) and EOP (xp, yp, UT1) differences between two VLBI solutions.

"""

import numpy as np
import matplotlib.pyplot as plt
from os import path
from read_eob import read_eob
from linear_regression import linear_regression
from eob_difference_statistics import eop_diff_stats, nut_diff_stats
from residual_statistics import stats_calc
from residual_plot import errorbarplot_res


# -----------------------------  FUNCTIONS -----------------------------
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


def eob_diff_calc(dat1, dat2):
    '''Calculate the position difference.

    Parameters
    ----------
    dat1, dat2 : list, containing
            X/Y/U : Xp/Yp/UT1, mas/mas/msec
            X_err/Y_err/U_err : formal uncertainty of Xp/Yp/UT1, mas/mas/msec
            P/E : dX/dY, mas
            P_err/E_err : formal uncertainty of P/E, mas
            cor : correlation coefficient

        The formal uncertainty of dx = x1 -x2 is calculated as
           dx_err = sqrt(x_err1^2 + x_err2^2)

        The covriance between dx and dy is calculated as
            covdxdy = corxy1 * x_err1 * y_err1 + corxy2 * x_err2 * y_err2
                see Appendix B in Mignard et al 2016

    Returns
    ----------
    dif : list, containing:
            dX/dY/dU : difference of Xp/Yp/UT1, uas/uas/us
            dX_err/dY_err/dU_err : formal uncertainty of Xp/Yp/UT1, uas
            covdXdY, covdXdU, covdYdU : covriance between dX/dY, dX/dU, dY/dU,
                                        uas^2
            dP/dE : difference of dX/dY, uas
            dP_err/dE_err : formal uncertainty of dP/dE, uas
            covdXdY, covdXdU, covdYdU : covriance between dP/dE, uas^2
    '''

    [X1, X_err1, Y1, Y_err1, U1, U_err1,
     corXY1, corXU1, corYU1,
     P1, P_err1, E1, E_err1, corPE1] = dat1
    [X2, X_err2, Y2, Y_err2, U2, U_err2,
     corXY2, corXU2, corYU2,
     P2, P_err2, E2, E_err2, corPE2] = dat2

    dX = (X1 - X2) * 1.e3
    dX_err = np.sqrt(X_err1**2 + X_err2**2) * 1.e3

    dY = (Y1 - Y2) * 1.e3
    dY_err = np.sqrt(Y_err1**2 + Y_err2**2) * 1.e3

    dU = (U1 - U2) * 1.e3
    dU_err = np.sqrt(U_err1**2 + U_err2**2) * 1.e3

    covdXdY = (corXY1 * X_err1 * Y_err1 + corXY2 * X_err2 * Y_err2) * 1.e6
    covdXdU = (corXU1 * X_err1 * U_err1 + corXU2 * X_err2 * U_err2) * 1.e6
    covdYdU = (corYU1 * U_err1 * Y_err1 + corYU2 * U_err2 * Y_err2) * 1.e6

    dP = (P1 - P2) * 1.e3
    dP_err = np.sqrt(P_err1**2 + P_err2**2) * 1.e3

    dE = (E1 - E2) * 1.e3
    dE_err = np.sqrt(E_err1**2 + E_err2**2) * 1.e3

    covdPdE = (corPE1 * P_err1 * E_err1 + corPE2 * P_err2 * E_err2) * 1.e6

    return [dX, dX_err, dY, dY_err, dU, dU_err,
            covdXdY, covdXdU, covdYdU,
            dP, dP_err, dE, dE_err, covdPdE]


def EOB_diff(eob1, eob2, eopdif, nutdif):
    '''Calculate the EOP and nutation offsets differences.

    Here I use the convention P/E stand for the dX/dY, because I use the option 'APR_NUT' for getpar
    to fetch the result.
    '''

    [dbname1, obsnum1, tag_eop1,
     X1, X_err1, Y1, Y_err1, U1, U_err1,
     XR1, XR_err1, YR1, YR_err1, UR1, UR_err1,
     corXY1, corXU1, corYU1, corXUR1, corYUR1, corUUR1,
     tag_nut1, P1, P_err1, E1, E_err1, corPE1] = read_eob(eob1)

    [dbname2, obsnum2, tag_eop2,
     X2, X_err2, Y2, Y_err2, U2, U_err2,
     XR2, XR_err2, YR2, YR_err2, UR2, UR_err2,
     corXY2, corXU2, corYU2, corXUR2, corYUR2, corUUR2,
     tag_nut2, P2, P_err2, E2, E_err2, corPE2] = read_eob(eob2)

    [dX, dX_err, dY, dY_err, dU, dU_err,
     covdXdY, covdXdU, covdYdU,
     dP, dP_err, dE, dE_err, covdPdE] = eob_diff_calc(
        [X1, X_err1, Y1, Y_err1, U1, U_err1,
         corXY1, corXU1, corYU1,
         P1, P_err1, E1, E_err1, corPE1],
        [X2, X_err2, Y2, Y_err2, U2, U_err2,
         corXY2, corXU2, corYU2,
         P2, P_err2, E2, E_err2, corPE2])

    feopdif = open(eopdif, "w")
    fnutdif = open(nutdif, "w")

    print("# EOP differences between two catalogs:\n"
          "# %s\n# %s\n"
          "# Epoch  dXp  dXp_err  dYp  dYp_err  dU  dU_err     COV\n"
          "# Year   uas  uas      uas  uas      ms  ms    XpYp XpU YpU\n"
          % (eob1, eob2), file=feopdif)

    print("# Nutation differences between two catalogs:\n"
          "# %s\n# %s"
          "# Epoch  dX  dX_err  dY  dY_err  COV\n"
          "# Year   uas uas     uas uas     dXdY\n"
          % (eob1, eob2), file=fnutdif)

    for i in range(dX.size):
        print("%10.2f  %+8.1f  %8.1f  %+8.1f  %8.1f  %+8.1f  %8.1f  "
              "%14.1f  %14.1f  %14.1f  " %
              (tag_eop1[i],
               dX[i], dX_err[i],
               dY[i], dY_err[i],
               dU[i], dU_err[i],
               covdXdY[i], covdXdU[i], covdYdU[i]),
              file=feopdif)

        print("%10.2f  %+8.1f  %8.1f  %+8.1f  %8.1f  %14.1f" %
              (tag_nut1[i], dP[i], dP_err[i], dE[i], dE_err[i], covdPdE[i]),
              file=fnutdif)

    feopdif.close()
    fnutdif.close()


def eop_diff_stats(DIF_file, lab):
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

    # errorbarplot_res(epo, dX, dXerr, "$\Delta$Xp", "uas",
    #                  path.dirname(data_dir), "%s_dXp" % lab)
    # errorbarplot_res(epo, dY, dYerr, "$\Delta$Yp", "uas",
    #                  path.dirname(data_dir), "%s_dYp" % lab)
    # errorbarplot_res(epo, dU, dUerr, "$\Delta$UT1-UTC", "ms",
    #                  path.dirname(data_dir), "%s_dU" % lab)

    # errorbarplot_res(epo, dX, dXerr, "$\Delta$Xp", "uas",
    #                  path.dirname(data_dir), "%s_dXp" % lab,
    #                  -0.241, 1.796)
    # errorbarplot_res(epo, dY, dYerr, "$\Delta$Yp", "uas",
    #                  path.dirname(data_dir), "%s_dYp" % lab,
    #                  0.097, 2.172)
    # errorbarplot_res(epo, dU, dUerr, "$\Delta$UT1-UTC", "ms",
    #                  path.dirname(data_dir), "%s_dU" % lab,
    #                  9.492e-3, 116.766e-3)

    # Statistical calculation
    flog = open("%s/%s_eop.log" % (data_dir, lab), "w")

    # Statistics for all-time series
    print("For all-time:", file=flog)
    print("X pole (uas):", file=flog)

    print("Y pole: (uas)", file=flog)
    stats_calc(epo, dY, dYerr, flog)

    print("UT1-UTC (us):", file=flog)
    stats_calc(epo, dU * 1000, dUerr * 1000, flog)

    # Statistics for post-1995.0 series
    splitpoint = 1995.
    print("For post-%.1f:" % splitpoint, file=flog)

    con = (epo >= splitpoint)
    epon = epo[con]
    dXn = dX[con]
    dXerrn = dXerr[con]
    dYn = dY[con]
    dYerrn = dYerr[con]
    dUn = dU[con]
    dUerrn = dUerr[con]

    print("X pole (uas):", file=flog)
    stats_calc(epon, dXn, dXerrn, flog)
    slope, intercept = stats_calc(epo, dX, dXerr, flog)

    errorbarplot_res(epo, dX, dXerr, "$\Delta$Xp", "$\mu$as",
                     path.dirname(data_dir), "%s_dXp" % lab,
                     slope, intercept)

    print("Y pole: (uas)", file=flog)
    slope, intercept = stats_calc(epon, dYn, dYerrn, flog)

    errorbarplot_res(epo, dY, dYerr, "$\Delta$Yp", "$\mu$as",
                     path.dirname(data_dir), "%s_dYp" % lab,
                     slope, intercept)

    print("UT1-UTC (us):", file=flog)
    # slope, intercept = stats_calc(epon, dUn * 1000, dUerrn * 1000,
    #                               flog)
    # It seems no need to multiply a 1000 because the unit for UT1 is
    # already uas.
    slope, intercept = stats_calc(epon, dUn, dUerrn,
                                  flog)

    errorbarplot_res(epo, dU, dUerr, "$\Delta$UT1-UTC", "ms",
                     path.dirname(data_dir), "%s_dU" % lab,
                     slope/1.e3, intercept/1.e3)
    flog.close()


def nut_diff_stats(DIF_file, lab):
    '''
    '''

    # Load data
    ep_mjd, dX, dXerr, dY, dYerr = get_res_nut(DIF_file)
    epo = (ep_mjd - 51544.5) / 365.25 + 2000.0

    data_dir, data_fil = path.split(DIF_file)

    # Statistical calculation
    flog = open("%s/%s_nut.log" % (data_dir, lab), "w")

    # Statistics for all-time series
    print("For all-time:", file=flog)
    print("dX (uas):", file=flog)
    stats_calc(epo, dX, dXerr, flog)

    print("dY: (uas)", file=flog)
    stats_calc(epo, dY, dYerr, flog)

    # Statistics for post-1995.0 series
    splitpoint = 1995.
    print("For post-%.1f:" % splitpoint, file=flog)

    con = (epo >= splitpoint)
    epon = epo[con]
    dXn = dX[con]
    dXerrn = dXerr[con]
    dYn = dY[con]
    dYerrn = dYerr[con]

    print("dX (uas):", file=flog)
    slope, intercept = stats_calc(epon, dXn, dXerrn, flog)
    errorbarplot_res(epo, dX, dXerr, "$\Delta$dX", "$\mu$as",
                     path.dirname(data_dir), "%s_ddX" % lab,
                     slope, intercept)

    print("dY: (uas)", file=flog)
    slope, intercept = stats_calc(epon, dYn, dYerrn, flog)
    errorbarplot_res(epo, dY, dYerr, "$\Delta$dY", "$\mu$as",
                     path.dirname(data_dir), "%s_ddY" % lab,
                     slope, intercept)

    # Plot
    # plot_res(epo, dX, "dX", "mas", data_dir)
    # plot_res(epo, dY, "dY", "mas", data_dir)

    # errorbarplot_res(epo, dX, dXerr, "$\Delta$dX", "$\mu$as",
    #                  path.dirname(data_dir), "%s_ddX" % lab)
    # errorbarplot_res(epo, dY, dYerr, "$\Delta$dY", "$\mu$as",
    #                  path.dirname(data_dir), "%s_ddY" % lab)

    flog.close()


def eob_comparison(eob1, eob2, lab):
    '''
    '''

    # eob1 = "/home/nliu/solutions/test/GA/opa2017a_aprx.eob"
    # eob2 = "/home/nliu/solutions/test/GA/opa2017a_ga.eob"

    # log_dir = "/home/nliu/solutions/GalacticAberration/logs" # vlbi2
    log_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/logs"  # My MacOS

    eopdif = "%s/%s.eop" % (log_dir, lab)
    nutdif = "%s/%s.nut" % (log_dir, lab)

    EOB_diff(eob1, eob2, eopdif, nutdif)

    eop_diff_stats(eopdif, lab)

    nut_diff_stats(nutdif, lab)

# -------------------- MAIN ----------------------------------

# # Test
# eob1 = "/home/nliu/solutions/test/a1/result_a1.eob"
# eob2 = "/home/nliu/solutions/test/a2/result_a2.eob"
# eopdif = "/home/nliu/solutions/test/a1_a2_dif.eop"
# nutdif = "/home/nliu/solutions/test/a1_a2_dif.nut"

# eob1 = "/home/nliu/solutions/test/GA/opa2017a_aprx.eob"
# eob2 = "/home/nliu/solutions/test/GA/opa2017a_ga.eob"
# eopdif = "/home/nliu/solutions/test/GA/opa2017a_ga_dif.eop"
# nutdif = "/home/nliu/solutions/test/GA/opa2017a_ga_dif.nut"


# EOB_diff(eob1, eob2, eopdif, nutdif)

# eop_diff_stats(eopdif)

# nut_diff_stats(nutdif)
# --------------------------------- END --------------------------------
