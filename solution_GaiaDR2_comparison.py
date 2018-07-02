#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: solution_GaiaDR2_comparison.py
"""
Created on Fri May  4 16:20:21 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import cos
from os import path
from functools import reduce
import matplotlib.pyplot as plt

from read_sou import read_sou_pos
from VSH_analysis import vsh_analysis
from cross_match import pos_Xmatch, postional_difference_calc
from read_GaiaDR2 import read_gaiadr2_iers_position

# deg -> as
deg2as = 3.6e3


# -----------------------------  FUNCTIONS -----------------------------
def sol_Gaia_diff_calc(cat, cat_dif, sollab):
    '''Calculate the quasar position difference between two vlbi-solutions.

    Paramters
    ---------
    cat : filename with path of the .cat file of a VLBI solution
    cat_dif : filename with path into whom the psotional differences will
            be written.
    sollab : solution label.

    Returns
    ---------
    soucom : source names of common sources, IVS designation
    RA1n, Dec1n : source positions in the 1st catalog, degree
    dRAc, dDec : source position differences, micro-as
    dRAc_err, dDec_err : uncertainty of dRA*cos(Dec) and dDec, micro-as
    cov : covariance between dRAc and dDec, micro-as^2
    pos_sep : positional differences in micro-as
    X_a, X_d : normalized coordinate differences in RA / DC, unit-less
    X : Normalized separations, unit-less
    flgcom : source type flags
    '''

    # Read .cat file
    iers_name1 = np.genfromtxt(cat, dtype=str, usecols=(1,))
    ra1, dec1, ra_error1, dec_error1, ra_dec_corr1 = np.genfromtxt(
        cat, usecols=range(2, 7), unpack=True)

    # # mas -> uas
    # ra_error1 *= 1.e3
    # dec_error1 *= 1.e3

    # Read Gaia DR2 data
    [iers_name2, ra2, ra_error2,
     dec2, dec_error2, ra_dec_corr2] = read_gaiadr2_iers_position(
        "/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/gaiadr2_iers.fits")

    # # mas -> uas
    # ra_error2 *= 1.e3
    # dec_error2 *= 1.e3

    [soucom,
     RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
     RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n] = pos_Xmatch(
        iers_name1, ra1, ra_error1, dec1, dec_error1, ra_dec_corr1,
        iers_name2, ra2, ra_error2, dec2, dec_error2, ra_dec_corr2)

    [dRAc, dDec, dRAc_err, dDec_err, cov,
     pos_sep, X_a, X_d, X, X2] = postional_difference_calc(
        RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
        RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n)

    fdif = open(cat_dif, "w")

    print("# Position difference between two catalogs:\n"
          "# %s \n# GaiaDR2" % cat, file=fdif)
    print("# Sou  RA  Dec  dRAc  dRAcerr  dDec  dDecerr  COV  "
          "pos_sep  X_a  X_d  X  X2\n"
          "#      deg deg  uas   uas      uas   uas     uas^2"
          "uas\n"
          "# X_a/X_d are normalized \n"
          "# X/X2 are Mignard's and normal normalized seperations.\n"
          "# For flag, 'D', 'V', and 'N' are same as in ICRF2, "
          "'O' for non-icrf2 source", file=fdif)

    # for i, soui in enumerate(soucom):
    for (soui, RA1ni, Dec1ni, dRAci, dRAc_erri, dDeci, dDec_erri,
         covi, pos_sepi, X_ai, X_di, Xi) in zip(
            soucom, RA1n, Dec1n, dRAc, dRAc_err, dDec, dDec_err, cov,
            pos_sep, X_a, X_d, X):

        print("%9s  %14.10f  %+14.10f  %+10.4f  %10.4f  %+10.4f  %10.4f  "
              "%+14.6f  %10.4f  %+8.3f  %+8.3f  %8.3f" %
              (soui, RA1ni, Dec1ni,
               dRAci, dRAc_erri, dDeci, dDec_erri, covi,
               pos_sepi, X_ai, X_di, Xi),
              file=fdif)

    fdif.close()

    flg = np.empty_like(soucom)

    return [soucom, RA1n, Dec1n,
            dRAc, dRAc_err, dDec, dDec_err, cov,
            pos_sep, X_a, X_d, X, flg]


def solution_GaiaDR2_analysis(datadir, datafile, sollab):
    '''Comparison between the VLBI solution and Gaia DR1
    '''

    DiffData = sol_Gaia_diff_calc(
        "%s/%s" % (datadir, datafile),
        "%s/%s_Gaiadr2_dif.sou" % (datadir, datafile[:-4]),
        "%s_GaiaDR2" % sollab)

    # VSH function parameters
    print("# VSH analysis:")
    # vsh_analysis(cat_dif, datadir, "%s_GaiaDR1" % sollab)
    vsh_analysis(DiffData,
                 "%s/%s_Gaiadr2_dif.sou" % (datadir, datafile[:-4]),
                 label="%s_GaiaDR2" % sollab)


# --------------------------------- END --------------------------------
