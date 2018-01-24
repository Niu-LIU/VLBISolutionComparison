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
def calc_wrms(res, err=None):
    '''Calculate the (weighted) wrms of residuals.
    '''

    if err is None:
        mean = np.mean(res)
        resn = res - mean
        wrms = np.sqrt(np.dot(resn, resn) / (resn.size - 1))
    else:
        p = 1. / err
        mean = np.dot(res, p**2) / np.dot(p, p)
        resn = (res - mean) * p
        wrms = np.sqrt(np.dot(resn, resn) / np.dot(p, p))

    return mean, wrms


def elim_wrms(res, err=None, n=3.0):
    '''Using 3.0-sigma to eliminate the outlier and then calculate the wrms.
    '''

    mean, wrms = calc_wrms(res, err)

    # print(mean, wrms)

    cond = (np.fabs((res - mean)/err) <= n)
    resn = res[cond]

    # print(resn)

    if err is None:
        mean, wrms = calc_wrms(resn)
        return resn, mean, wrms, cond
    else:
        errn = err[cond]
        mean, wrms = calc_wrms(resn, errn)
        return resn, errn, mean, wrms, cond


def stats_calc(t, res, err, flog):
    '''Statistical result of residual.
    '''

    resn, errn, mean, wrms, cond = elim_wrms(res, err)

    # slope, intercept, r_value, p_value, std_err = stats.linregress(
    #     t[cond], resn)

    tn = t[cond]
    t0 = 2005.0

    par, err, outlier, cor = linear_regression(
        tn - t0, resn, errn)
    slope, intercept = par

    print("# weighted\n",
          "# Mean : %.3f\n" % mean,
          # "# Mean1 : %f\n" % mean1,
          "# WRMS : %.3f\n" % wrms,
          "# Slope : %.3f\n" % slope,
          "# Intercept : %.3f" % intercept,
          file=flog)

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
