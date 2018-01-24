#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: c04_apr_calc.py
"""
Created on Sun Dec 10 18:37:34 2017

@author: Neo(liuniu@smail.nju.edu.cn)

Calculate the apriori EOP based C04 series.

"""

import numpy as np
from scipy.interpolate import CubicSpline
import sys
from read_c04 import read_c04


__all__ = {"cubspl", "itpl_eop", "itpl_nut", "calc_eob"}

# -----------------------------  FUNCTIONS -----------------------------


def cubspl(x, xs, ys):
    '''For series (xi,yi), get y(x) using cubic spline.


    Parameters
    ----------
    xs : array, float
        time series of X component
    ys : array, float
        time series of Y component
    x : float
        X position need to be interpolated

    Returns
    ---------
    Interpolated Y value.
    '''

    cs = CubicSpline(xs, ys)
    return cs(x)


def itpl_eop(epo, mjd20, Xp20, Yp20, U20):
    '''Get interpolated EOP at a certain epoch.


    Parameters
    ----------
    epo : float
        epoch to be interpolated
    mjd20 : array, float
        20-point epoch series, Julian day
    Xp20 : array, float
        20-point Xp series, mas
    Yp20 : array, float
        20-point Yp series, mas
    U20 : array, float
        20-point UT1-UTC series, ms

    Returns
    ----------
    Xpitpl : float
        interpolated Xp value, mas
    Ypitpl : float
        interpolated Yp value, mas
    Uitpl : float
        interpolated UT1-UTC value, ms
    '''

    # X pole
    Xpitpl = cubspl(epo, mjd20, Xp20)
    # Y pole
    Ypitpl = cubspl(epo, mjd20, Yp20)
    # U pole
    Uitpl = cubspl(epo, mjd20, U20)

    return Xpitpl, Ypitpl, Uitpl


def itpl_nut(epo, mjd20, dX20, dY20):
    '''Get interpolated EOP at a certain epoch.


    Parameters
    ----------
    epo : float
        epoch to be interpolated
    mjd20 : array, float
        20-point epoch series, Julian day
    dX20 : array, float
        20-point dX series, mas
    dY20 : array, float
        20-point dY series, mas

    Returns
    ----------
    dXitpl : float
        interpolated dX value, mas
    dYitpl : float
        interpolated dY value, mas
    '''

    # dX
    dXitpl = cubspl(epo, mjd20, dX20)
    # dY
    dYitpl = cubspl(epo, mjd20, dY20)

    return dXitpl, dYitpl


# def calc_eop(epoitpl, mjd, Xp, Yp, U, XpErr, YpErr, UErr):
def calc_eob(EOB_file):
    '''Use Cubic Spline Interpolation to calculate the EOP.


    For the target eopch, We use the 10 successive points in front of
    and behind it in the time series (total 20 points) to do the Cubic
    Spline Interpolation.

    Parameters
    ----------
#     epoitpl : float, array
#         itpget epoch
#     mjd : array, float
#         epoch in modified Julian date
#     Xp : array, float
#         Xp position series of CIP in ITRS, mas
#     Yp : array, float
#         Yp position series of CIP in ITRS, mas
#     U : array, float
#         UT1 - UTC series, ms
#     XErr : array, float
#         formal uncertainty of Xp, mas
#     YErr : array, float
#         formal uncertainty of Yp, mas
#     UErr : array, float
#         formal uncertainty of U, ms

    Returns
    ----------
    itp_X : array, float
        X position of CIP in ITRS at itpget epoch, mas
    itp_Y : array, float
        Y position series of CIP in ITRS, mas
    itp_U : array, float
        UT1 - UTC series, ms
#     itp_XErr : array, float
#         formal uncertainty of itp_X, mas
#     itp_YErr : array, float
#         formal uncertainty of itp_Y, mas
#     itp_UErr : array, float
#         formal uncertainty of itp_U, ms
    '''

    print("# ---------- BEGIN ----------")
    # EOB_file = "/home/nliu/solutions/test/x1/result_x1.eob"
    # EOB_file = "test.eob"
    teopitpl, tnutitpl = np.genfromtxt(EOB_file, usecols=(0, 28), unpack=True)

    # C04_file = "/Users/Neo/Astronomy/Data/EOP_IERS/eopc04_IAU2000.62-now.txt"
    # C04_file = "/home/nliu/Data/eop/eopc04_IAU2000.62-now"

    # vlbi2
    # C04_file = "/home/nliu/Data/eop/eopc04_IAU2000"

    # My MacOS
    C04_file = "/Users/Neo/Astronomy/Data/EOP_IERS/eopc04_IAU2000"

    print("C04 file: %s" % C04_file)
    mjd, Xp, Yp, U, dX, dY, XpErr, YpErr, UErr, dXErr, dYErr = read_c04(
        C04_file)
    print("Read C04 file: done!")

    OPEOP_file = EOB_file.replace("eob", "eop_c04_apr")
    print("# Output file: %s" % OPEOP_file)
    fopeop = open(OPEOP_file, "w")

    OPNUT_file = EOB_file.replace("eob", "nut_c04_apr")
    print("# Output file: %s" % OPNUT_file)
    fopnut = open(OPNUT_file, "w")

    print("\n# Begin to interpolate!")

    # For EOP
    epo_ind = np.searchsorted(mjd, teopitpl)

    for i, ind in enumerate(epo_ind):
        if ind == 0 or ind >= mjd.size:
            # normally it won't happen!!
            print("The epoch %f was too early or too late"
                  " for the C04 series." % teopitpl[i])
            sys.exit()
        elif ind < 9 or mjd.size - ind < 10:
            # In this case we will use less datapoints.
            pnum = np.min(ind, mjd.size - ind)
        else:
            pnum = 10

        mjd20 = mjd[ind - pnum + 1: ind + pnum + 1]
        Xp20 = Xp[ind - pnum + 1: ind + pnum + 1]
        Yp20 = Yp[ind - pnum + 1: ind + pnum + 1]
        U20 = U[ind - pnum + 1: ind + pnum + 1]

        Xpitpl, Ypitpl, Uitpl = itpl_eop(
            teopitpl[i], mjd20, Xp20, Yp20, U20)

        print("%12.6f  %+8.3f  %+8.3f  %+8.3f" %
              (teopitpl[i], Xpitpl, Ypitpl, Uitpl),
              file=fopeop)

    # For NUT
    epo_ind = np.searchsorted(mjd, tnutitpl)

    for i, ind in enumerate(epo_ind):
        if ind == 0 or ind == mjd.size:
            # normally it won't happen!!
            print("The epoch %f was too early or too late"
                  " for the C04 series." % tnutitpl[i])
            sys.exit()
        elif ind < 9 or mjd.size - ind < 10:
            # In this case we will use less datapoints.
            pnum = np.min(ind, mjd.size - ind)
        else:
            pnum = 10

        indmin, indmax = ind - pnum, ind + pnum

        mjd20 = mjd[indmin: indmax]
        dX20 = dX[indmin: indmax]
        dY20 = dY[indmin: indmax]

        dXitpl, dYitpl = itpl_nut(
            tnutitpl[i], mjd20, dX20, dY20)

        print("%12.6f  %+8.3f  %+8.3f" %
              (tnutitpl[i], dXitpl, dYitpl),
              file=fopnut)
        # print("%12.6f:   %+8.3f mas  %+8.3f mas  %+8.3f ms" %
        #       (epoitpl[i], Xpitpl, Ypitpl, Uitpl))

    fopeop.close()
    fopnut.close()
    print("# ----------  END  ----------")

# # read_c04()
# if len(sys.argv) == 2:
#     calc_eop(sys.argv[1])
# else:
#     calc_eop("test.eob")
# --------------------------------- END --------------------------------
