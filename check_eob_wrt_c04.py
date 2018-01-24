#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: check_eob_wrt_c04.py
"""
Created on Tue Dec 12 16:19:31 2017

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import sys
from read_eob import read_eob
from dat import Dat
from c04_apr_calc import calc_eob
from eob_difference_statistics import eop_diff_stats, nut_diff_stats


# -----------------------------  FUNCTIONS -----------------------------
def calc_eop_diff(EOB_file):
    '''
    '''

    APR_file = EOB_file.replace("eob", "eop_c04_apr")
    DIF_file = EOB_file.replace("eob", "eop_c04_dif")

    [dbname, obsnum, tag_eop, Xp, Xp_err, Yp, Yp_err, U, U_err,
     XR, XR_err, YR, YR_err, UR, UR_err,
     corXY, corXU, corYU, corXUR, corYUR, corUUR,
     tag_nut, dX, dX_err, dY, dY_err, cordXdY] = read_eob(EOB_file)

    mjdapr, Xpapr, Ypapr, Uapr = np.genfromtxt(APR_file, unpack=True)

    if all(tag_eop != mjdapr):
        print("These two files do not match!")
        exit()

    # IAI - UTC
    dat = [Dat(epoch) for epoch in tag_eop]
    dat = np.array(dat) * 1000  # unit: ms
    U = dat + U

    # compare
    fout = open(DIF_file, "w")
    print("# EOP difference wrt C04:\n"
          "# Epoch  dXp  errXp  dYp  errYp  dU  errU\n"
          "# mjd    mas  mas    mas  mas    ms  ms",
          file=fout)

    # difference
    for (i, mjd) in enumerate(mjdapr):
        print("%13.6f" % mjd, "  %+10.3f  %7.3f" * 3 %
              (Xp[i]-Xpapr[i], Xp_err[i],
               Yp[i]-Ypapr[i], Yp_err[i],
               U[i]-Uapr[i], U_err[i]), file=fout)

    fout.close()

    return DIF_file


def calc_nut_diff(EOB_file):
    '''
    '''

    APR_file = EOB_file.replace("eob", "nut_c04_apr")
    DIF_file = EOB_file.replace("eob", "nut_c04_dif")

    [dbname, obsnum, tag_eop, Xp, Xp_err, Yp, Yp_err, U, U_err,
     XR, XR_err, YR, YR_err, UR, UR_err,
     corXY, corXU, corYU, corXUR, corYUR, corUUR,
     tag_nut, dX, dX_err, dY, dY_err, cordXdY] = read_eob(EOB_file)

    mjdapr, dXapr, dYapr = np.genfromtxt(APR_file, unpack=True)

    if all(tag_nut != mjdapr):
        print("These two files do not match!")
        exit()

    # compare
    fout = open(DIF_file, "w")
    print("# EOP difference wrt C04:\n"
          # "# %s\n# %s" % (cat1, cat2), file=fdif)
          # print(
          "# Epoch  ddX  errdX  ddY  errdY\n"
          "# mjd    mas  mas    mas  mas",
          file=fout)

    # difference
    for (i, mjd) in enumerate(mjdapr):
        print("%13.6f" % mjd, "  %+10.3f  %7.3f" * 2 %
              (dX[i] - dXapr[i], dX_err[i],
               dY[i] - dYapr[i], dY_err[i]), file=fout)

    fout.close()

    return DIF_file


def check_eob_wrt_c04(EOB_file):
    '''
    '''

    # Calculate the C04 apriori EOP and Nutation
    calc_eob(EOB_file)

    # EOP difference
    EOPDIF_file = calc_eop_diff(EOB_file)

    # # Statistical calculation
    eop_diff_stats(EOPDIF_file)

    # nutation difference
    NUTDIF_file = calc_nut_diff(EOB_file)

    # Statistical calculation
    nut_diff_stats(NUTDIF_file)

# -------------------- MAIN ----------------------------------
# if len(sys.argv) == 1:
#     EOB_file = 'result/test.eob'
# else:
#     EOB_file = sys.argv[1]
# dbname, tag, obsnum, P, P_err, E, E_err = read_nut(datafile)
# EOB_file = "/home/nliu/solutions/test/GA/opa2017a_aprx.eob"
# EOB_file = "/home/nliu/solutions/test/GA/opa2017a_ga.eob"
# EOB_file = "/home/nliu/solutions/GalacticAberration/opa2018a_ga/opa2018a_ga.eob"

# check_eob_wrt_c04(EOB_file)
# --------------------------------- END --------------------------------
