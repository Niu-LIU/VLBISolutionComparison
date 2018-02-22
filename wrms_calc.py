#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: wrms_calc.py
"""
Created on Wed Feb 14 11:02:00 2018

@author: Neo(liuniu@smail.nju.edu.cn)

This script is used for calculating the pre-fit wrms,
post-fit wrms, reduced-chi square, and standard deviation.

"""

import numpy as np
from functools import reduce


# -----------------------------  FUNCTIONS -----------------------------
def calc_wrms(x, err=None):
    '''Calculate the (weighted) wrms of x series.

    wrms = sqrt(sum( (xi-mean)^2/erri^2 ) / sum( 1.0/erri^2 ))
         if weighted,
         = sqrt(sum( (xi-mean)^2/erri^2 ) / (N-1))
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
    '''

    if err is None:
        mean = np.mean(x)
        xn = x - mean
        wrms = np.sqrt(np.dot(xn, xn) / (xn.size - 1))
    else:
        p = 1. / err
        mean = np.dot(x, p**2) / np.dot(p, p)
        xn = (x - mean) * p
        wrms = np.sqrt(np.dot(xn, xn) / np.dot(p, p))

    return mean, wrms


def calc_chi2(x, err, reduced=False, deg=0):
    '''Calculate the (reduced) Chi-square.


    Parameters
    ----------
    x : array, float
        residuals
    err : array, float
        formal errors of residuals
    reduced : boolean
        True for calculating the reduced chi-square
    deg : int
        degree of freedom

    Returns
    ----------
    (reduced) chi-square
    '''

    wx = x / err
    chi2 = np.dot(wx, wx)

    if reduced:
        if deg:
            return chi2 / (x.size - deg)
        else:
            print("# ERROR: the degree of freedom cannot be 0!")
    else:
        return chi2


def calc_2Dchi2(x, errx, y, erry, covxy, reduced=False):
    '''Calculate the 2-Dimension (reduced) Chi-square.


    Parameters
    ----------
    x : array, float
        residuals of x
    errx : array, float
        formal errors of x
    x : array, float
        residuals of x
    errx : array, float
        formal errors of x
    covxy : array, float
        summation of covariance between x and y
    reduced : boolean
        True for calculating the reduced chi-square

    Returns
    ----------
    (reduced) chi-square
    '''

    Qxy = np.zeros_like(x)

    for i, (xi, errxi, yi, erryi, covxyi
            ) in enumerate(zip(x, errx, y, erry, covxy)):

        wgt = np.linalg.inv(np.mat([[errxi**2, covxyi],
                                    [covxyi, erryi**2]]))

        Xmat = np.array([xi, yi])

        Qxy[i] = np.sqrt(reduce(np.dot, (Xmat, wgt, Xmat)))

    if reduced:
        return np.sum(Qxy) / (x.size - 2)
    else:
        return np.sum(Qxy)

# --------------------------------- END --------------------------------
