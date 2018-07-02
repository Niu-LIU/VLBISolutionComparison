#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: residual_statistics.py
"""
Created on Fri Jan 12 09:56:18 2018

@author: Neo(liuniu@smail.nju.edu.cn)


N. Liu, 11/04/2018: add the statistics of residuals after removing the
                    linear trend;
"""

import numpy as np
from linear_regression import linear_regression

__all__ = {"calc_wrms", "elim_wrms", "stats_calc"}


# -----------------------------  FUNCTIONS -----------------------------
def calc_wrms(x, err=None):
    '''Calculate the (weighted) wrms and std of x series.

    Standard deviation
    std = sqrt(sum( (xi-mean)^2/erri^2 ) / sum( 1.0/erri^2 ))
         if weighted,
         = sqrt(sum( (xi-mean)^2/erri^2 ) / (N-1))
         otherwise.

    Weighted root mean square
    wrms = sqrt(sum( xi^2/erri^2 ) / sum( 1.0/erri^2 ))
         if weighted,
         = sqrt(sum( xi^2/erri^2 ) / (N-1))
         otherwise.

    Parameters
    ----------
    x : array, float
        Series
    err : array, float, default is None.

    Returns
    ----------
    mean : float
    wrms : float
        weighted rms
    std : float
        standard deviation
    '''

    if err is None:
        mean = np.mean(x)
        xn = x - mean
        wrms = np.sqrt(np.dot(xn, xn) / (xn.size - 1))
        std = np.sqrt(np.dot(x, x) / (x.size - 1))
    else:
        p = 1. / err
        mean = np.dot(x, p**2) / np.dot(p, p)
        xn = (x - mean) * p
        wrms = np.sqrt(np.dot(xn, xn) / np.dot(p, p))
        std = np.sqrt(np.dot(x*p, x*p) / np.dot(p, p))

    return mean, wrms, std


def elim_wrms(res, err=None, n=3.0):
    '''Using 3.0-sigma to eliminate the outlier and then calculate the wrms.
    '''

    mean, wrms, std = calc_wrms(res, err)

    # print(mean, wrms)

    cond = (np.fabs((res - mean)/err) <= n)
    resn = res[cond]

    # print(resn)

    if err is None:
        mean, wrms, std = calc_wrms(resn)
        return resn, mean, wrms, std, cond
    else:
        errn = err[cond]
        mean, wrms, std = calc_wrms(resn, errn)
        return resn, errn, mean, wrms, std, cond


def stats_calc(t, res, err, flog):
    '''Statistical result of residual.
    '''
    # beginning epoch for calculate the trend
    t0 = 2000.0

    resn, errn, mean, wrms, std, cond = elim_wrms(res, err)

    # slope, intercept, r_value, p_value, std_err = stats.linregress(
    #     t[cond], resn)

    tn = t[cond]

    par, parerr, outlier, cor = linear_regression(
        tn - t0, resn, errn)
    slope, intercept = par
    slperr, itperr = parerr

    # print("# weighted\n",
    #       "# Mean      : %.3f\n" % mean,
    #       "# Std       : %.3f\n" % std,
    #       "# WRMS      : %.3f\n" % wrms,
    #       "# Slope     : %.3f +/- %.3f\n" % (slope, slperr),
    #       "# Intercept : %.3f  " % intercept,
    #       file=flog)

    print("# weighted\n",
          "# Mean      : %.2f\n" % mean,
          "# Std       : %.2f\n" % std,
          "# WRMS      : %.2f\n" % wrms,
          "# Slope     : %.2f +/- %.2f\n" % (slope, slperr),
          "# Intercept : %.2f  " % intercept,
          file=flog)

    print("STAS_ALL ", mean, std, wrms, slope, slperr,
          file=flog)

    # Add the statistics after removing the linear trend;
    res1 = res - slope * (t - t0)
    # resn1, errn1, mean1, wrms1, std1, cond1 = elim_wrms(res1, err)
    # tn1 = t[cond1]
    # par1, parerr1, outlier1, cor1 = linear_regression(
    #     tn1 - t0, resn1, errn1)
    # slope1, intercept1 = par1
    # slperr1, itperr1 = parerr1

    resn1 = resn - slope * (tn - t0)
    errn1 = errn
    _, _, mean1, wrms1, std1, cond1 = elim_wrms(res1, err)
    # tn1 = t[cond1]
    par1, parerr1, outlier1, cor1 = linear_regression(
        tn - t0, resn1, errn1)
    slope1, intercept1 = par1
    slperr1, itperr1 = parerr1

    print("# After removing linear trend:\n",
          "# Mean      : %.2f\n" % mean1,
          "# Std       : %.2f\n" % std1,
          "# WRMS      : %.2f\n" % wrms1,
          "# Slope     : %.2f +/- %.2f\n" % (slope1, slperr1),
          "# Intercept : %.2f\n" % intercept1,
          file=flog)

    print("STAS_AFTER ", mean1, std1, wrms1, slope1, slperr1,
          file=flog)

    return slope, intercept

    # # No weights
    # resn, wrms, cond = elim_wrms(res)
    # mean = np.average(resn)

    # slope, intercept, r_value, p_value, std_err = stats.linregress(
    #     t[cond], resn)

    # print("# None-wgt\n",
    #       "# Mean : %f\n" % mean,
    #       "#  RMS : %f\n" % wrms,
    #       "# Slope : %f\n" % slope,
    #       "# Intercept : %f" % intercept,
    #       file=flog)

# --------------------------------- END --------------------------------
