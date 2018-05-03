#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: nor_sep.py
"""
Created on Mon Apr  2 16:52:06 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from functools import reduce
sqrt = np.sqrt

__all__ = ['nor_sep_calc', 'vlbi_gaia_sep']


# -----------------------------  FUNCTIONS -----------------------------
def nor_sep_calc(dRA, dRA_err, dDC, dDC_err, C):
    '''Calculate the normalized seperation.

    Parameters
    ----------
    dRA/dDC : Right Ascension / Declination differences in micro-as
    dRA_err/dDC_err : formal uncertainty of dRA/dDC in micro-as
    C : correlation coeffient between dRA and dDC.

    Returns
    ----------
    ang_sep : angular seperation, in micro-as
    X_a / X_d : normalized coordinate differences in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    # Angular seperations
    ang_sep = sqrt(dRA**2 + dDC**2)

    # Normalised coordinate difference
    X_a = dRA / dRA_err
    X_d = dDC / dDC_err

    # Normalised separation
    X = np.zeros_like(X_a)

    for i, (X_ai, X_di, Ci) in enumerate(zip(X_a, X_d, C)):
        wgt = np.linalg.inv(np.mat([[1, Ci],
                                    [Ci, 1]]))
        Xmat = np.array([X_ai, X_di])
        X[i] = sqrt(reduce(np.dot, (Xmat, wgt, Xmat)))

    return ang_sep, X_a, X_d, X


def vlbi_gaia_sep(RA1, DC1, RA1_err, DC1_err, Cor1,
                  RA2, DC2, RA2_err, DC2_err, Cor2,
                  arccof=None):
    '''Calculate the normalized seperation between VLBI and Gaia positions.


    Parameters
    ----------
    RA / DC : Right Ascension / Declination, degress
    e_RA / e_DC : formal uncertainty of RA / DC, micro-as
    Cor : correlation coeffient between RA and DC.
    arccof : cos(Dec.)

    Note: suffix 'G' stands for GaiaDR1 and I for VLBI catalog.

    Returns
    ----------
    ang_sep : angular seperation in micro-as
    X_a / X_d : normalized seperation in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    if arccof is None:
        arccof = np.cos(np.deg2rad(DC1))

    # deg -> uas
    dRA = (RA1 - RA2) * 3.6e9 * arccof
    dRA_err = sqrt(RA1_err**2 + RA2_err**2)
    dDC = (DC1 - DC2) * 3.6e9
    dDC_err = sqrt(DC1_err**2 + DC2_err**2)

    # Correlation coefficient of combined errors
    C = (RA1_err * DC1_err * Cor1 +
         RA2_err * DC2_err * Cor2) / (dRA_err * dDC_err)

    # Normalised separation
    ang_sep, X_a, X_d, X = nor_sep_calc(dRA, dDC, dRA_err, dDC_err, C)

    return ang_sep, X_a, X_d, X


def test_code():
    '''Verify the codes using ICRF2 and GaiaDR1.
    '''

    datfil = "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/gaiadr1_icrf2_1.dat"

    RA1, DC1, RA1_err, DC1_err, Cor1 = np.genfromtxt(
        datfil, usecols=range(3, 8), delimiter="|", unpack=True)
    RA2, DC2, RA2_err, DC2_err, Cor2 = np.genfromtxt(
        datfil, usecols=range(8, 13), delimiter="|", unpack=True)
    sourcename = np.genfromtxt(datfil, usecols=(0,), dtype=str,
                               delimiter="|")

    ang_sep, X_a, X_d, X = vlbi_gaia_sep(RA1, DC1, RA1_err, DC1_err, Cor1,
                                         RA2, DC2, RA2_err, DC2_err, Cor2)

    for (soui, ang_sepi, X_ai, X_di, Xi) in zip(
            sourcename, ang_sep, X_a, X_d, X):

        print(soui, ang_sepi, X_ai, X_di, Xi)


# --------------------------------- END --------------------------------
# test_code()
# The results are the same with that in Mignard et al. 2017.
# --------------------------------- END --------------------------------
