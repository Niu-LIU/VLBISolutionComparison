#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:52:21 2017

@author: Neo


Dec 8, N. Liu : modify the bugs in Sect. 3.10 of the help document of getpar_02
                25 ----> cor(U, UR)
                26 ----> cor(X, UR)
                27 ----> cor(Y, UR)

11 Jan 2018, N. Liu : set the zero-formal-uncertainty to 0.999999


Retrieve the estimates of X pole coordinate, Y pole coordinate, UT1-TAI
angle, UT1 rate, daily offsets of nutation angles as well as their formal
uncertainties and correlations from .eob file which is generated by
the program getpar.

.eob file contains series of the estimates of X pole coordinate,
Y pole coordinate, UT1-TAI angle, UT1 rate, daily offsets of nutation angles
as well as their formal uncertainties and correlations. Time tag and database
name is attached to each line. .EOB format is an extension of the IERS EOP
format.

   File contains lines of three types:
1) Comment. The first character is #. Header comments contain some information
   about solution.

2) Header. The first two symbols are blank. Header lines contain titles of the
   columns

3) Estimates.


  1    1-1    A1     ---     Usage flag
  2    3-14   F12.6  days    Modified Julian date of the TDT time tag for
                             pole coordinates and UT1
  3   16-25   A10    ---     Database name
  4   27-32   A6     ---     IVS session code (if available)
  5   34-41   F8.6   arcsec  The estimate of X pole coordinate
  6   43-50   F8.6   arcsec  The estimate of Y pole coordinate
  7   52-62   F11.7  sec     The UT1-TAI function
  8   64-71   F8.3   mas     Adjustment of the nutation in longitude angle with
                                        respect to IAU 1980 nutation expansion
  9   73-80   F8.3   mas     Adjustment of the nutation in obliquity angle with
                                        respect to IAU 1980 theory
 10   82-90   F9.6   asc/day The estimate of X pole rate
 11   92-100  F9.6   asc/day The estimate of Y pole rate
 12  102-108  F7.4   ms/day  The estimate of UT1 rate
 13  110-117  F8.6   arcsec  Formal uncertainty of X pole coordinate
 14  119-126  F8.6   arcsec  Formal uncertainty of Y pole coordinate
 15  128-136  F9.7   sec     Formal uncertainty of UT1-UTC function
 16  138-144  F7.3   mas     Formal uncertainty of nutation in longitude angle
 17  146-152  F7.3   mas     Formal uncertainty of nutation in obliquity angle
 18  154-162  F9.6   asc/day Formal uncertainty of X pole rate
 19  164-172  F9.6   asc/day Formal uncertainty of Y pole rate
 20  174-180  F7.4   asc/day Formal uncertainty of UT1 rate
 21  182-187  F6.4   --      Correlation between the estimates of X-pole
                                          positions and Y-pole position
 22  189-194  F6.4   --      Correlation between the estimates of X-pole
                                         positions and UT1-TAI angle
 23  196-201  F6.4   --      Correlation between the estimates of Y-pole
                                         positions and UT1-TAI angle
 24  203-208  F6.4   --      Correlation between the estimates of nutation in
                                         longitude and nutation in obliquity
 25  210-215  F6.4   --      Correlation between the estimates of X-pole
                                          positions and UT1 rate
 26  217-222  F6.4   --      Correlation between the estimates of Y-pole
                                         positions and UT1-TAI date
 27  224-229  F6.4   --      Correlation between the estimates of
                                         UT1-TAI angle UT1 rate
 28  231-235  F5.2   hours   Session duration
 29  237-243  F7.2   psec    Weighted root mean square of postfit residuals
 30  245-250  I6     --      Number of used observations in the session
 31  252-263  F12.6  days    Modified Julian date for nutation at TDT time
                             scale
 32  265-328  A64    --      The network configuration line. Consists of
                             two characters IVS station codes listed
                             in alphabetic order for stations that participated
                             in the experiment and supplied the data that have
                             been used in processing this experiment.

If the specific parameter was not estimated in this experiment, the field
for its value and formal uncertainty is replaced by filler: $$$$$$. The filler
takes entire field.

"""

import numpy as np
import matplotlib.pyplot as plt


__all__ = {"zero_err", "read_eob"}


# ------------------------------  FUNCTIONS  ---------------------------
def zero_err(err):
    '''set zero-error to 0.999999
    '''

    return np.where(err == 0, 0.999999, err)


def read_eob(datafile):
    '''Retrieve the result from .eob file.

    Parameters
    ----------
    datafile : string
        name of data file

    Returns
    ----------
    dbname : array, string
       database name with leading dollar sign
    tag_eop : array, float
        EOP time flag, modified Julian data
    tag_nut : array, float
        Nutation time flag, modified Julian data
    obsnum : array, int
        number of observations used
    X : array, float
        X-pole coordinate, mas
    X_err : array, float
        formal uncertainty of X, mas
    Y : array, float
        Y-pole coordinate, mas
    Y_err : array, float
        formal uncertainty of Y, mas
    U : array, float
        UT1 - TAI, msec
    U_err : array, float
        formal uncertainty of U, msec
    XR : array, float
        X-pole coordinate rate, mas/day
    XR_err : array, float
        formal uncertainty of XR, mas/day
    YR : array, float, mas/day
        Y-pole coordinate rate, mas/day
    YR_err : array, float
        formal uncertainty of YR, mas/day
    UR : array, float
        UT1 - TAI rate, msec/day
    UR_err : array, float
        formal uncertainty of UR, msec/day
    P : array, float
        nutation in longitude, mas
    P_err : array, float
        formal uncertainty of P, mas
    E : array, float
        nutation on obliquity, mas
    E_err : array, float
        formal uncertainty of E, mas
    corXY : array, float
        correlation between X and Y
    corXU : array, float
        correlation between X and U
    corYU : array, float
        correlation between Y and U
    corPE : array, float
        correlation between X and Y
    corXUR : array, float
        correlation between X and UR
    corYUR : array, float
        correlation between Y and UR
    corUUR : array, float
        correlation between U and UR
    '''

    dbname = np.genfromtxt(datafile, dtype=str, usecols=(1,))

    tag_eop, tag_nut = np.genfromtxt(datafile, usecols=(0, 28), unpack=True)

    obsnum = np.genfromtxt(datafile, dtype=int, usecols=(27,))

    X, Y, U, P, E, XR, YR, UR = np.genfromtxt(datafile,
                                              usecols=np.arange(2, 10),
                                              missing_values='*'*7,
                                              filling_values=0.0,
                                              unpack=True)

    [X_err, Y_err, U_err,
     P_err, E_err,
     XR_err, YR_err, UR_err] = np.genfromtxt(datafile,
                                             usecols=np.arange(10, 18),
                                             missing_values='*'*7,
                                             filling_values=1.0e6,
                                             unpack=True)

    # Avoid zero-error, 0.0 -> 0.999999
    X_err = zero_err(X_err)
    Y_err = zero_err(Y_err)
    U_err = zero_err(U_err)
    P_err = zero_err(P_err)
    E_err = zero_err(E_err)
    XR_err = zero_err(XR_err)
    YR_err = zero_err(XR_err)
    UR_err = zero_err(UR_err)

    [corXY, corXU, corYU, corPE, corUUR, corXUR, corYUR
     ] = np.genfromtxt(datafile, usecols=np.arange(18, 25), unpack=True)

    # as -> mas or sec -> msec
    X, Y, U = X * 1000.0, Y * 1000.0, U * 1000.0
    X_err, Y_err, U_err = X_err * 1000.0, Y_err * 1000.0, U_err * 1000.0

    # as/day -> mas/day or sec/day -> msec/day
    XR, YR = XR * 1000, YR * 1000
    XR_err, YR_err, UR_err = XR_err * 1000, YR_err * 1000, UR_err * 1000

    return [dbname, obsnum, tag_eop, X, X_err, Y, Y_err, U, U_err,
            XR, XR_err, YR, YR_err, UR, UR_err,
            corXY, corXU, corYU, corXUR, corYUR, corUUR,
            tag_nut, P, P_err, E, E_err, corPE]


def read_eops(eops_file):
    """Read data from .eops file.

    Parameters
    ----------
    datafile : string
        name of data file

    Returns
    ----------
    db_name : array, string
       database name with leading dollar sign
    mjd : array, float
        reference epoch, modified Julian data
    obs_num : array, int
        number of observations used
    sess_len : array, float
        session duration
    rms : array, float
        postfit rms delay in ps
    xp : array, float
        X-pole coordinate, mas
    xp_err : array, float
        formal uncertainty of xp, mas
    yp : array, float
        Y-pole coordinate, mas
    yp_err : array, float
        formal uncertainty of yp, mas
    ut : array, float
        UT1 - TAI, mas
    ut_err : array, float
        formal uncertainty of ut, mas
    xpr : array, float
        X-pole coordinate rate, mas/day
    xpr_err : array, float
        formal uncertainty of xpr, mas/day
    ypr : array, float, mas/day
        Y-pole coordinate rate, mas/day
    ypr_err : array, float
        formal uncertainty of ypr, mas/day
    utr : array, float
        UT1 - TAI rate, mas/day
    utr_err : array, float
        formal uncertainty of utr, mas/day
    dx : array, float
        Celestial pole offset dX wrt IAU 2006, mas
    dx_err : array, float
        formal uncertainty of dX, mas
    dy : array, float
        Celestial pole offset dY wrt IAU 2006, mas
    dy_err : array, float
        formal uncertainty of dY, mas
    xp_yp_corr : array, float
        correlation between xp and yp
    xp_ut_corr : array, float
        correlation between xp and ut
    yp_ut_corr : array, float
        correlation between yp and ut
    dx_dy_corr : array, float
        correlation between dx and dy
    """

    db_name = np.genfromtxt(eops_file, usecols=(11,), dtype=str)

    mjd, sess_len, rms = np.genfromtxt(
        eops_file, usecols=(0, 18, 29), unpack=True)

    obs_num = np.genfromtxt(
        eops_file, usecols=(16,), dtype=int, unpack=True)

    # Earth orientation parameters
    xp, yp, ut, dx, dy = np.genfromtxt(
        eops_file, usecols=range(1, 6), unpack=True)
    xp_err, yp_err, ut_err, dx_err, dy_err = np.genfromtxt(
        eops_file, usecols=range(6, 11), unpack=True)

    # Correlation
    xp_yp_corr, xp_ut_corr, yp_ut_corr, dx_dy_corr = np.genfromtxt(
        eops_file, usecols=range(12, 16), unpack=True)

    # EOP rate
    xpr, ypr, utr = np.genfromtxt(
        eops_file, usecols=range(19, 22), unpack=True)
    xpr_err, ypr_err, utr_err = np.genfromtxt(
        eops_file, usecols=range(24, 27), unpack=True)

    return [db_name, mjd, obs_num, sess_len, rms,
            xp, yp, ut, dx, dy,
            xp_err, yp_err, ut_err, dx_err, dy_err,
            xp_yp_corr, xp_ut_corr, yp_ut_corr, dx_dy_corr,
            xpr, ypr, utr, xpr_err, ypr_err, utr_err]


# ------------------------------  MAIN BODY  ---------------------------
# # Retrieve estimates.
# dat1 = read_eops("/Users/Neo/Astronomy/Data/VLBISolutions/"
#                  "vlbi2_server/GA-eop/opa2018r.eops-SL")
# dat2 = read_eops("/Users/Neo/Astronomy/Data/VLBISolutions/"
#                  "vlbi2_server/GA-eop/opa2018r.eops")
# for (dat1i, dat2i) in zip(dat1, dat2):
#     print(dat1i[0], dat2i[0])
# ------------------------------ END -----------------------------------
