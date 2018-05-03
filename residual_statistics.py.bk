#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: residual_statistics.py
"""
Created on Fri Jan 12 09:56:18 2018

@author: Neo(liuniu@smail.nju.edu.cn)
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

    resn, errn, mean, wrms, std, cond = elim_wrms(res, err)

    # slope, intercept, r_value, p_value, std_err = stats.linregress(
    #     t[cond], resn)

    tn = t[cond]
    t0 = 2007.0

    par, err, outlier, cor = linear_regression(
        tn - t0, resn, errn)
    slope, intercept = par

    print("# weighted\n",
          "# Mean      : %.3f\n" % mean,
          "# Std       : %f  \n" % std,
          "# WRMS      : %.3f\n" % wrms,
          "# Slope     : %.3f\n" % slope,
          "# Intercept : %.3f  " % intercept,
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
