#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: solution_GaiaDR1_comparison.py
"""
Created on Thu Jan 11 17:14:07 2018

@author: Neo(liuniu@smail.nju.edu.cn)

History
N. Liu, 25 Jan 2018: add functions sub_posplot and post_diff_plot
                     to plot the source position differences;
                     add a parameter 'sollab' to function
                     'sol_icrf2_diff_calc';
N. Liu, 06 Feb 2018: add a new function 'nor_sep_calc';
                     minor changes in functions
                     'solution_Gaia_diff_analysis' and
                     'sol_Gaia_diff_calc'
N. Liu, 09 Apr 2018: minor changes in functions
                     'solution_Gaia_diff_analysis' and
                     'sol_Gaia_diff_calc'

"""

import numpy as np
from os import path
from functools import reduce
import matplotlib.pyplot as plt
from read_sou import read_sou_pos
from nor_sep import nor_sep_calc
from Transformation_cat import catalog_transfor
from VSH_analysis import vsh_analysis
cos = np.cos
deg2as = 3.6e3

__all__ = ["read_gaidadr1", "position_taken", "Xmatch",
           "position_diff_calc", "nor_sep_calc",
           "sub_posplot", "post_diff_plot",
           "sol_Gaia_diff_calc", "solution_Gaia_diff_analysis"]


# -----------------------------  FUNCTIONS -----------------------------
def read_gaiadr1(datafile):
    '''Read Gaia DR1 position

    Note:   1) RA_err = err(RA*cos(DC))
            2) unit for RA/DC : deg
                    for RA_err/DC_err : uas

    Parameters
    ----------
    datafile :
        file name and full path of Gaia DR1 data

    Returns
    ----------
    sourcename :
        source name (ICRF designation) of quasars
    RAdeg / DCdeg :
        Right Ascension / Declination, degreees
    RA_erruas / DC_erruas :
        formal uncertainty of RA*cos(Dec) / DC, micro-as.
    cor :
        correlation coeffient between RA and DC.
    Flag :
        ICRF2 classification of quasars.
    '''

    sourcename, Flag = np.genfromtxt(
        datafile, usecols=(1, 13), dtype=str,
        delimiter="|",
        unpack=True)
    RAdeg, DCdeg, RA_erruas, DC_erruas, cor = np.genfromtxt(
        datafile, usecols=range(3, 8), delimiter="|", unpack=True)

    return sourcename, RAdeg, DCdeg, RA_erruas, DC_erruas, cor, Flag


def position_taken(index, RA, RAc_err, Dec, Dec_err, cor):
    '''Extract the elements from array at specific index.

    Parameters
    ----------
    index :
        the indice corresponding to common sources
    RA / Dec :
        Right Ascension / Declination for all sources, degreees
    RAc_err / Dec_err :
        formal uncertainty of RA*cos(Dec) / Dec for all sources, micro-as.
    cor :
        correlation coeffient between RA and Dec for all sources.

    Returns
    ----------
    RAn / Decn :
        Right Ascension / Declination for common sources, degreees
    RAc_errn / Dec_err :
        formal uncertainty of RA*cos(Dec) / Dec for common sources, micro-as.
    corn :
        correlation coeffient between RA and Dec for common sources.
    '''

    RAn = np.take(RA, index)
    RAc_errn = np.take(RAc_err, index)
    Decn = np.take(Dec, index)
    Dec_errn = np.take(Dec_err, index)
    corn = np.take(cor, index)

    return RAn, RAc_errn, Decn, Dec_errn, corn


def Xmatch(sou1, RA1, RAc_err1, Dec1, Dec_err1, cor1,
           sou2, RA2, RAc_err2, Dec2, Dec_err2, cor2, flg):
    '''Crossmatch between Gaia and VLBI catalogs.

    Parameters
    ----------
    sou :
        source name (ICRF designation)
    RA / Dec :
        Right Ascension / Declination, degreees
    RAc_err / DC_err :
        formal uncertainty of RA*cos(Dec) / Dec, micro-as.
    cor :
        correlation coeffient between RA and Dec.

    Returns
    ----------
    soucom :
        name (ICRF designation) of common sources
    RAn / Decn :
        Right Ascension / Declination for common sources, degreees
    RAc_errn / Dec_err :
        formal uncertainty of RA*cos(Dec) / Dec for common sources, micro-as.
    cor :
        correlation coeffient between RA and Dec for common sources.
    '''

    soucom = []
    index1 = []
    index2 = []
    flgcom = []

    for i, soui in enumerate(sou1):
        indarr = np.where(sou2 == soui)[0]

        if indarr:
            soucom.append(soui)
            index1.append(i)
            j = indarr[0]
            index2.append(j)
            flgcom.append(flg[j])

    RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n = position_taken(
        index1, RA1, RAc_err1, Dec1, Dec_err1, cor1)
    RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n = position_taken(
        index2, RA2, RAc_err2, Dec2, Dec_err2, cor2)

    flgcom = np.array(flgcom)

    return [soucom, flgcom,
            RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
            RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n]


def position_diff_calc(dat1, dat2):
    '''Calculate the position difference.

    Parameters
    ----------
    dat1, dat2 : list, containing
            RA : Right ascension, deg
            RA_err : formal uncertainty of RA, uas
            DC : declination, deg
            DC_err : formal uncertainty of DC, uas
            cor : correlation coefficient between RA and DC

    Returns
    ----------
    dif : list, containing:
            dRAc : difference of \Detal_RA*cos(Dec), uas
            dRAc_err : formal uncertainty of dRAc=sqrt(RAc_err1^2 + RAc_err2^2), uas
            dDec : difference of DC, uas
            dDec_err : formal uncertainty of dDec=sqrt(Dec_err1^2 + Dec_err2^2), uas
            cov : covriance between dRAc and dDec, uas ^2
            cof : correlation coeffient between dRAc and dDec
                see Appendix B in Mignard et al 2016
    '''

    RA1, RAc_err1, Dec1, Dec_err1, cor1 = dat1
    RA2, RAc_err2, Dec2, Dec_err2, cor2 = dat2

    # deg -> uas
    dRAc = (RA1 - RA2) * 3.6e9 * cos(np.deg2rad(Dec1))
    dDec = (Dec1 - Dec2) * 3.6e9

    dRAc_err = np.sqrt(RAc_err1**2 + RAc_err2**2)
    dDec_err = np.sqrt(Dec_err1**2 + Dec_err2**2)
    cov = (cor1 * RAc_err1 * Dec_err1 + cor2 * RAc_err2 * Dec_err2)
    cof = cov / dRAc_err / dDec_err

    return [dRAc, dRAc_err, dDec, dDec_err, cov, cof]


def sub_posplot(ax0, ax1, dRA, dRA_err, dDC, dDC_err, DC, soutp, mk):
    '''
    '''
    ax0.plot(DC, dRA, mk, markersize=1, label=soutp)
    ax1.plot(DC, dDC, mk, markersize=1, label=soutp)


def post_diff_plot(dRA, dRA_err, dDC, dDC_err, DC, tp, main_dir, lab):
    '''
    '''

    print("%d source for plot.\n" % DC.size)

    unit = "mas"
    lim = 5
    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    flgs = ['D', 'V', 'N']
    soutps = ['Defining', 'VCS', 'Non-VCS']
    mks = ['r.', 'b*', 'g<']
    for flg, soutp, mk in zip(flgs, soutps, mks):
        cond = tp == flg
        sub_posplot(ax0, ax1,
                    dRA[cond], dRA_err[cond],
                    dDC[cond], dDC_err[cond],
                    DC[cond], soutp, mk)

    ax0.set_ylabel("$\Delta \\alpha^* (%s)$" % unit)
    ax0.set_ylim([-lim, lim])

    # ax1.errorbar(DC, dDC, yerr=dDC_err, fmt='.', elinewidth=0.01,
    #              markersize=1)
    ax1.set_ylabel("$\Delta \delta (%s)$" % unit)
    ax1.set_ylim([-lim, lim])
    ax1.set_xlabel("Dec. (deg)")
    ax1.set_xlim([-90, 90])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    plt.legend(fontsize='xx-small')
    plt.savefig("%s/plots/%s_sou_dif.eps" % (main_dir, lab))
    plt.close()


def sol_Gaia_diff_calc(cat, catdif, sollab):
    '''Calculate the quasar position difference between two vlbi-solutions.

    Paramters
    ---------
    cat : filename with path of the .sou file of a VLBI solution
    catdif : filename with path into whom the psotional differences will
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

    sou1, RA1, RAc_err1, Dec1, Dec_err1, cor1 = read_sou_pos(
        cat, unit_deg=True, arcerr=True)
    # mas -> uas
    RAc_err1 = RAc_err1 * 1.e3
    Dec_err1 = Dec_err1 * 1.e3

    sou2, RA2, Dec2, RAc_err2, Dec_err2, cor2, flg = read_gaiadr1(
        # "/home/nliu/Data/gaiadr1_icrf2_1.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/GaiaDR1_cds/"
        "gaiadr1_icrf2_1.dat")  # My MacOS

    [soucom, flgcom,
     RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
     RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n] = Xmatch(
        sou1, RA1, RAc_err1, Dec1, Dec_err1, cor1,
        sou2, RA2, RAc_err2, Dec2, Dec_err2, cor2, flg)

    dRAc, dRAc_err, dDec, dDec_err, cov, cof = position_diff_calc(
        [RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n],
        [RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n])

    pos_sep, X_a, X_d, X = nor_sep_calc(dRAc, dRAc_err, dDec, dDec_err, cof)

    fdif = open(catdif, "w")

    print("# Position difference between two catalogs:\n"
          "# %s \n# GaiaDR1" % cat, file=fdif)
    print("# Sou  RA  Dec  dRAc  dRAcerr  dDec  dDecerr  COV  "
          "pos_sep  X_a  X_d  X  flag\n"
          "#      deg deg  uas   uas      uas   uas     uas^2"
          "uas\n"
          "# X_a/X_d are normalized \n"
          "# For flag, 'D', 'V', and 'N' are same as in ICRF2, "
          "'O' for non-icrf2 source", file=fdif)

    # for i, soui in enumerate(soucom):
    for (soui, RA1ni, Dec1ni, dRAci, dRAc_erri, dDeci, dDec_erri,
         covi, pos_sepi, X_ai, X_di, Xi, flgcomi) in zip(
            soucom, RA1n, Dec1n, dRAc, dRAc_err, dDec, dDec_err, cov,
            pos_sep, X_a, X_d, X, flgcom):

        print("%9s  %14.10f  %14.10f  %+10.1f  %10.1f  %+10.1f  %10.1f  "
              "%14.1f  %+10.1f  %+8.3f  %+8.3f  %+8.3f  %s" %
              (soui, RA1ni, Dec1ni,
               dRAci, dRAc_erri, dDeci, dDec_erri, covi,
               pos_sepi, X_ai, X_di, Xi,
               flgcomi),
              file=fdif)

    fdif.close()

    # main_dir = "/home/nliu/solutions/GalacticAberration" # vlbi2
    # main_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3"  # My MacOS
    main_dir = path.dirname(cat)
    post_diff_plot(dRAc / 1.e3, dRAc_err / 1.e3,
                   dDec / 1.e3, dDec_err / 1.e3, Dec1n,
                   flgcom, main_dir, sollab)

    return [soucom, RA1n, Dec1n,
            dRAc, dRAc_err, dDec, dDec_err, cov,
            pos_sep, X_a, X_d, X, flgcom]


def solution_Gaia_diff_analysis(datadir, datafile, sollab):
    '''Comparison between the VLBI solution and Gaia DR1
    '''

    DiffData = sol_Gaia_diff_calc(
        "%s/%s" % (datadir, datafile),
        "%s/%s_Gaiadr1_dif.sou" % (datadir, datafile[:-4]),
        "%s_GaiaDR1" % sollab)
    # DiffData = [soucom, RAdeg, DCdeg, dRA, dRA_err, dDC, dDC_err, cov,
    #             pos_sep, X_a, X_d, X, flgcom]

    # IERS transformation parameter
    print("# IERS transformation:")
    # catalog_transfor(catdif, datadir, "%s_GaiaDR1" % sollab)
    catalog_transfor(DiffData,
                     "%s/%s_Gaiadr1_dif.sou" % (datadir, datafile[:-4]),
                     datadir, "%s_GaiaDR1" % sollab)

    # VSH function parameters
    print("# VSH analysis:")
    # vsh_analysis(catdif, datadir, "%s_GaiaDR1" % sollab)
    vsh_analysis(DiffData,
                 "%s/%s_Gaiadr1_dif.sou" % (datadir, datafile[:-4]),
                 datadir, "%s_GaiaDR1" % sollab)


def test_code():
    solution_Gaia_diff_analysis("/Users/Neo/Astronomy/Works"
                                "/201711_GDR2_ICRF3/solutions/"
                                "GalacticAberration/opa2018a_ref00s",
                                "opa2018a_ref00s.sou", "x")

# test_code()
# --------------------------------- END --------------------------------
