#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_icrf2.py
"""
Created on Mon Apr  9 11:45:31 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import genfromtxt, cos


# -----------------------------  FUNCTIONS -----------------------------
def read_icrf2(icrf2_file="/Users/Neo/Astronomy/Data/catalogs/"
               "icrf/icrf2.dat", unit_pos_as=False, arc_err_flag=True):
    '''Read ICRF2 catalog.

    Parameters
    ----------
    icrf2_file : string
        ICRF2 data file
    unit_pos_as : Boolean
        flag to determine if the unit of RA./Dec. is arc-sec.
    arc_err_flag : Boolean
        flag to determine if returning the formal uncertainty of RA or
        RA*cos(Dec.). True for RA*cos(Dec.) while False for RA.


    Returns
    ----------
    icrfn : array of string
        ICRF designation of source name;
    ivsn : array of string
        IVS designation of source name;
    iersn : array of string
        IERS designation of source name;
    RA/Dec : array of float
        Right Ascension/Declination in degree;
    e_RA/e_DE : array of float
        formal error of RA/Dec in milliarcsecond;
    cor : array of float between [-1, +1]
        correlation coefficient between RA and Dec in ICRF2 solution;
    mepo : array of float
        mean epoch of source position in Modified Julian day;
    Flag : array of character
        flag of source classification in ICRF2 catalog.
    '''

    icrfn, ivsn, iersn, Flag = genfromtxt(
        icrf2_file, usecols=(0, 1, 2, 3), dtype=str, unpack=True)
    RAh, RAm, RAs, Decd, Decm, Decs, e_RAs, e_DEas, cor = genfromtxt(
        icrf2_file, usecols=range(4, 13), unpack=True)

    mepo = genfromtxt(icrf2_file, usecols=(13,))

# determine the sign of Declination
    strDecs = genfromtxt(icrf2_file, usecols=(7,), dtype=str)
    signstr = [x[0] for x in strDecs]
    Dec_sign = np.where(np.array(signstr) == '-', -1, 1)

# calculate the position
    RA = (RAh + RAm / 60.0 + RAs / 3600) * 15  # degree
    Dec = Decd + Dec_sign * (Decm / 60.0 + Decs / 3600)  # degree

# unit: as -> mas
    if arc_err_flag:
        e_RA = e_RAs * 15e3 * cos(np.deg2rad(Dec))
    else:
        e_RA = e_RAs * 15e3
    e_DE = e_DEas * 1.0e3

    if unit_pos_as:
        deg2as = 3.6e3
        RA = RA * deg2as
        Dec = Dec * deg2as

    return icrfn, ivsn, iersn, RA, Dec, e_RA, e_DE, cor, mepo, Flag

# --------------------------------- END --------------------------------
