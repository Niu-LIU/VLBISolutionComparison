#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: mv_icrf2.py
"""
Created on Sun Feb  4 18:52:31 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Since ICRF2 source positions are referred to different mean epoch,
and there is a systematic effect which can induce an apparent proper
motion field, it is, therefore, neccessary to propagate the ICRF2 positions
 to the same epoch.

"""

import numpy as np
from write_solvesrc import write_solvesrc, write_NNRS
from comparison_with_GaiaDR1 import comparison_with_Gaia
from declination_bias_plot import dec_bias_wrtICRF2
cos = np.cos
sin = np.sin
mjd_j2000 = 51544.5
uas2deg = 1.0 / 3.6e9


# -----------------------------  FUNCTIONS -----------------------------


def GA_pars_init():
    '''A priori value for Galactic Aberration Effect
    '''

    A = 5.0
    # A = 5.6  # ICRF3 working group recommendation value
    # alpha_G = 267.0
    # delta_G = -29.0
    alpha_G = 266.4
    delta_G = -28.9

    return A, alpha_G, delta_G


# def GA_apriori(A, alpha_G, delta_G):
def GA_apriori():
    '''Calculate the glide vector.

    Parameters
    ----------
    A : float
        amplitude of GA effect, micro-arcsec/yr
    alpha_G, delta_G : float
        direction of Galactic center, degree

    Return
    ----------
    g : array of (3,), float
        glide vector (g1, g2, g3), micro-arcsec/yr
    '''

    A, alpha_G, delta_G = GA_pars_init()

    alpha_G = np.deg2rad(alpha_G)
    delta_G = np.deg2rad(delta_G)

    g1 = A * cos(alpha_G) * cos(delta_G)
    g2 = A * sin(alpha_G) * cos(delta_G)
    g3 = A * sin(delta_G)

    return np.array([g1, g2, g3])


# def GA_effect_calc(RA, DC, g, flag="arc"):
def GA_effect_calc(RA, DC, ts, flag="arc"):
    '''Calculate the positional offset caused by GA effect


    Parameters
    ----------
    RA, DC : array, float
        Right ascension / Declination of radio sources, degree
    g : array of (3,0), float
        glide vector
    flag : string
        determine whether return the dRA or dRA * cos(DC)
    ts : array, float
        time span between Mean epoch and J2015.0/J2000.0

    Returns
    ----------
    dRA, dDC : array, float
        positional offsets caused by GA effect
    '''

    RAr = np.deg2rad(RA)
    DCr = np.deg2rad(DC)

    g1, g2, g3 = GA_apriori()

    if flag is "arc":
        dRA = (-g1 * sin(RAr) + g2 * cos(RAr)) * ts * uas2deg
    else:
        dRA = (-g1 * sin(RAr) + g2 * cos(RAr)) / cos(DCr) * ts * uas2deg

    dDC = (-g1 * cos(RAr) * sin(DCr) - g2 * sin(RAr) * sin(DCr)
           + g3 * cos(DCr)) * ts * uas2deg

    return dRA + RA, dDC + DC


def read_icrf2():
    '''Read ICRF2 file.

    Parameters
    ----------
    None

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

    icrf2_fil = "/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat"

    icrfn, ivsn, iersn, Flag = np.genfromtxt(
        icrf2_fil, usecols=(0, 1, 2, 3), dtype=str, unpack=True)
    RAh, RAm, RAs, Decd, Decm, Decs, e_RAs, e_DEas, cor = np.genfromtxt(
        icrf2_fil, usecols=range(4, 13), unpack=True)

    mepo = np.genfromtxt(icrf2_fil, usecols=(13,))

# determine the sign of Declination
    strDecs = np.genfromtxt(icrf2_fil, usecols=(7,), dtype=str)
    signstr = [x[0] for x in strDecs]
    Dec_sign = np.where(np.array(signstr) == '-', -1, 1)

# calculate the position
    RA = (RAh + RAm / 60.0 + RAs / 3600) * 15  # degree
    Dec = Decd + Dec_sign * (Decm / 60.0 + Decs / 3600)  # degree

# unit: as -> mas
    e_RA = e_RAs * 15e3
    e_DE = e_DEas * 1.0e3

    return icrfn, ivsn, iersn, RA, Dec, e_RA, e_DE, cor, mepo, Flag


def move_icrf2(year):
    ''' 'Move' ICRF2 mean-epoch to Year(MJD).

    Parametes
    ---------
    year : float
        Target epoch

    Returns
    ---------
    None
    '''

    icrfn, ivsn, iersn, RA, DC, e_RA, e_DC, cor, mepo, Flag = read_icrf2()

    # Time span between mean eooch and year.
    ts = year - 2000 - (mepo - mjd_j2000) / 365.25
    # print(ts)

    RAn, DCn = GA_effect_calc(RA, DC, ts, flag="nonarc")

    # Write data to SOLVE format.
    fsrc = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
                "icrf2_%.0f.src" % year, "w")
    # fsrc = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
    # "icrf2_%.0f_1.src" % year, "w")  # set A=5mas/yr

    # Write header
    print("$$ This file consists of SOLVE-format positions for 3414 "
          "ICRF2 sources with ref epoch J%.1f." % year, file=fsrc)
    # Write positions
    write_solvesrc(ivsn, RAn, DCn, e_DC,
                   "icrf2 mean epoch to J2015.0", fsrc)
    # Close file.
    fsrc.close()

    # Write data for further comparison
    fout = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
                "icrf2_%.0f.dat" % year, "w")
    # fout = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
    # "icrf2_%.0f_1.dat" % year, "w")  # set A=5.6 muas/yr
    print("# Move ICRF2 mean epoch to J%.1f\n"
          "# ICRF IVS IERS RA DC errRA errDC cor flg"
          "#               (deg)    (mas)        D/N/V " % year,
          file=fout)
    lfmt = "%16s %8s %8s %12.3f %+13.3f %8.4f %8.4f %+5.3f %s"
    for (icrfni, ivsni, iersni, RAni, DCni, e_RAi, e_DCi, cori, flgi
         ) in zip(icrfn, ivsn, iersn, RAn, DCn, e_RA, e_DC, cor, Flag):
        print(lfmt % (icrfni, ivsni, iersni, RAni, DCni,
                      e_RAi, e_DCi, cori, flgi), file=fout)
    fout.close()

    # print(icrfn, RAn, e_RA * 1.e3,
    #       DCn, e_DC * 1.e3, cor)

    # Declination wrt ICRF2
    dec_bias_wrtICRF2(ivsn, DCn, e_DC * 1.e3, 'icrf2%.0f' % year)

    # Comparison with GaiaDR1
    comparison_with_Gaia(icrfn, RAn, e_RA * 1.e3,
                         DCn, e_DC * 1.e3, cor,
                         "icrf%.0f" % year,
                         "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3"
                         "/data/icrf2_%.0f.dat" % year)
    # comparison_with_Gaia(icrfn, RAn, e_RA * 1.e3,
    #                      DCn, e_DC * 1.e3, cor,
    #                      "icrf%.0f_1" % year,
    #                      "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3"
    #                      "/data/icrf2_%.0f.dat" % year)


move_icrf2(2000.0)
move_icrf2(2015.0)
# --------------------------------- END -------------------------------
