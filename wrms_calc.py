#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: wrms_calc.py
"""
Created on Wed Feb 14 11:02:00 2018

@author: Neo(liuniu@smail.nju.edu.cn)

This script is used for calculating the pre-fit wrms,
post-fit wrms, reduced-chi square, and standard deviation.

3 Mar 2018, Niu : now function 'calc_wrms' also computes the standard
                  deviation

"""

import numpy as np
from functools import reduce

__all__ = ["calc_wrms", "calc_chi2", "calc_2Dchi2"]


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
        std = np.sqrt(np.dot(xn, xn) / (xn.size - 1))
        wrms = np.sqrt(np.dot(x, x) / (x.size - 1))
    else:
        p = 1. / err
        mean = np.dot(x, p**2) / np.dot(p, p)
        xn = (x - mean) * p
        std = np.sqrt(np.dot(xn, xn) / np.dot(p, p))
        wrms = np.sqrt(np.dot(x*p, x*p) / np.dot(p, p))

    return mean, wrms, std


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

        wgt = np.linalg.inv(np.array([[errxi**2, covxyi],
                                      [covxyi, erryi**2]]))

        Xmat = np.array([xi, yi])

        Qxy[i] = reduce(np.dot, (Xmat, wgt, Xmat))

        if Qxy[i] < 0:
            print(Qxy[i])

    if reduced:
        return np.sum(Qxy) / (x.size - 2)
    else:
        return np.sum(Qxy)

# --------------------------------- END --------------------------------
