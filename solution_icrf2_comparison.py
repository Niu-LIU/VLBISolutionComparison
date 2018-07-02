#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: sol_icrf2_diff.py
"""
Created on Thu Jan  4 10:37:01 2018

@author: Neo(liuniu@smail.nju.edu.cn)

History
N. Liu, 25 Jan 2018: add functions 'sub_posplot' and 'post_diff_plot'
                     to plot the source position differences;
                     add a parameter 'sollab' to function
                     'sol_icrf2_diff_calc';
N. Liu, 06 Feb 2018: add a new function 'nor_sep_calc';
                     minor changes in functions
                     'solution_Gaia_diff_analysis' and
                     'sol_Gaia_diff_calc'

"""

import numpy as np
import matplotlib.pyplot as plt
from os import path
from functools import reduce
from read_sou import read_sou_pos
from nor_sep import nor_sep_calc
from Transformation_cat import catalog_transfor
from VSH_analysis import vsh_analysis
cos = np.cos

__all__ = ["read_icrf2", "position_taken", "Xmatch",
           "position_diff_calc", "nor_sep_calc",
           "sub_posplot", "post_diff_plot",
           "solution_icrf2_diff_analysis", "sol_icrf2_diff_calc"]


# -----------------------------  FUNCTIONS -----------------------------
def read_icrf2(datafile):
    '''Read icrf2 position

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
        datafile, usecols=(1, 3), dtype=str, unpack=True)
    RAh, RAm, RAs, Decd, Decm, Decs, e_RAs, e_DEas, cor = np.genfromtxt(
        datafile, usecols=range(4, 13), unpack=True)

# determine the sign of Declination
    strDecs = np.genfromtxt(datafile, usecols=(7,), dtype=str)
    signstr = [x[0] for x in strDecs]
    Dec_sign = np.where(np.array(signstr) == '-', -1, 1)

# calculate the position
    RA = (RAh + RAm / 60.0 + RAs / 3600) * 15
    Dec = Decd + Dec_sign * (Decm / 60.0 + Decs / 3600)  # degree
# unit: s/as -> uas
    e_RA = e_RAs * 15e6 * cos(np.deg2rad(Dec))
    e_DE = e_DEas * 1.0e6
    return sourcename, RA, Dec, e_RA, e_DE, cor, Flag


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
            dRAc : difference of Detal_RA*cos(Dec), uas
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

    return [dRAc, dRAc_err, dDec,  dDec_err, cov, cof]


def sub_posplot(ax0, ax1, dRA, dRA_err, dDC, dDC_err, DC, soutp, mk):
    '''
    '''
    ax0.plot(DC, dRA, mk, markersize=1, label=soutp)
    ax1.plot(DC, dDC, mk, markersize=1, label=soutp)


def post_diff_plot(dRA, dRA_err, dDC, dDC_err, DC, tp, main_dir, sollab):
    '''
    '''

    print("%d source for plot.\n" % DC.size)

    unit = "mas"
    lim = 2.5
    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    # ax0.errorbar(DC, dRA, yerr=dRA_err, fmt='.', elinewidth=0.01,
    #              markersize=1)
    # ax0.plot(DC, dRA, '.', markersize=1)
    # ax1.plot(DC, dDC, '.', markersize=1)

    flgs = ['D', 'V', 'N']
    soutps = ['Defining', 'VCS', 'Non-VCS']
    mks = ['r.', 'b*', 'g<']
    for flg, soutp, mk in zip(flgs, soutps, mks):
        cond = tp == flg
        sub_posplot(ax0, ax1,
                    dRA[cond], dRA_err[cond],
                    dDC[cond], dDC_err[cond],
                    DC[cond], soutp, mk)

    ax0.set_ylabel("$\\Delta \\alpha^* (%s)$" % unit)
    ax0.set_ylim([-lim, lim])

    # ax1.errorbar(DC, dDC, yerr=dDC_err, fmt='.', elinewidth=0.01,
    #              markersize=1)
    ax1.set_ylabel("$\\Delta \\delta (%s)$" % unit)
    ax1.set_ylim([-lim, lim])
    ax1.set_xlabel("Dec. (deg)")
    ax1.set_xlim([-90, 90])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

    plt.legend(fontsize='xx-small')
    plt.savefig("%s/plots/%s_sou_dif.eps" % (main_dir, sollab))
    plt.close()


def sol_icrf2_diff_calc(cat, catdif, sollab):
    '''Calculate the quasar position difference between two vlbi-solutions.

    Paramters
    ---------
    cat : filename with path of the .cat file of a VLBI solution
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

    # sou1, RA1, RAc_err1, Dec1, Dec_err1, cor1 = read_sou_pos(
    #     cat, unit_deg=True, arcerr=True)
    sou1 = np.genfromtxt(cat, dtype=str, usecols=(0,))
    RA1, Dec1, RAc_err1, Dec_err1, cor1 = np.genfromtxt(
        cat, usecols=range(2, 7), unpack=True)

    # # mas -> uas
    # RAc_err1 = RAc_err1 * 1.e3
    # Dec_err1 = Dec_err1 * 1.e3

    #
    # ra_err = RAc_err1[RAc_err1 == 0.0]
    # sou = sou1[RAc_err1 == 0.0]
    # print(ra, sou)
    # exit()

    sou2, RA2, Dec2, RAc_err2, Dec_err2, cor2, flg = read_icrf2(
        # "/home/nliu/Data/icrf2.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat")  # My MacOS

    [soucom, flgcom,
     RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
     RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n] = Xmatch(
        sou1, RA1, RAc_err1, Dec1, Dec_err1, cor1,
        sou2, RA2, RAc_err2, Dec2, Dec_err2, cor2, flg)

    dRAc, dRAc_err, dDec, dDec_err, cov, cof = position_diff_calc(
        [RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n],
        [RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n])

    pos_sep, X_a, X_d, X = nor_sep_calc(
        dRAc, dRAc_err, dDec, dDec_err, cof)
    # pos_sep, X_a, X_d, X, X2 = nor_sep_calc(
    #     dRAc, dRAc_err, dDec, dDec_err, cof)

    fdif = open(catdif, "w")

    print("# Position difference between two catalogs:\n"
          "# %s \n# ICRF2" % cat, file=fdif)
    print("# Sou  RA  DC  dRA  dRAerr  dDec  dDCerr  COV  "
          "pos_sep  X_a  X_d  X  flag\n"
          "#      deg deg uas  uas     uas  uas     uas^2"
          "uas\n"
          "# X_a/X_d are normalized \n"
          "# For flag, 'D', 'V', and 'N' are same as in ICRF2, "
          "'O' for non-icrf2 source", file=fdif)

    # for i, soui in enumerate(soucom):
    # for (soui, RA1ni, Dec1ni, dRAci, dRAc_erri, dDeci, dDec_erri,
    #      covi, pos_sepi, X_ai, X_di, Xi, X2i, flgcomi) in zip(
    #         soucom, RA1n, Dec1n, dRAc, dRAc_err, dDec, dDec_err, cov,
    #         pos_sep, X_a, X_d, X, X2, flgcom):
    #     print("%9s  %14.10f  %14.10f  %+8.1f  %8.1f  %+8.1f  %8.1f  "
    #           "%14.1f  %+10.1f  %+8.3f  %+8.3f  %8.3f  %8.3f  %s" %
    #           (soui, RA1ni, Dec1ni, dRAci, dRAc_erri, dDeci, dDec_erri, covi,
    #            pos_sepi, X_ai, X_di, Xi, X2i, flgcomi), file=fdif)

    for (soui, RA1ni, Dec1ni, dRAci, dRAc_erri, dDeci, dDec_erri,
         covi, pos_sepi, X_ai, X_di, Xi, flgcomi) in zip(
            soucom, RA1n, Dec1n, dRAc, dRAc_err, dDec, dDec_err, cov,
            pos_sep, X_a, X_d, X, flgcom):
        print("%9s  %14.10f  %14.10f  %+8.1f  %8.1f  %+8.1f  %8.1f  "
              "%14.1f  %+10.1f  %+8.3f  %+8.3f  %8.3f  %s" %
              (soui, RA1ni, Dec1ni, dRAci, dRAc_erri, dDeci, dDec_erri, covi,
               pos_sepi, X_ai, X_di, Xi, flgcomi), file=fdif)

    fdif.close()

    # main_dir = "/home/nliu/solutions/GalacticAberration" # vlbi2
    # main_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3"  # My MacOS
    # main_dir = path.dirname(cat)
    # post_diff_plot(dRAc / 1.e3, dRAc_err / 1.e3,
    #                dDec / 1.e3, dDec_err / 1.e3, Dec1n,
    #                flgcom, main_dir, sollab)

    return [soucom, RA1n, Dec1n,
            dRAc, dRAc_err, dDec, dDec_err, cov,
            pos_sep, X_a, X_d, X, flgcom]


def solution_icrf2_analysis(datadir, datafile, sollab):
    '''Compare the VLBI solution with ICRF2 solution.

    The comparison is done by two methods:
        1) Use the transformation equation recommended by IERS group;
        2) Use degree-2 VSH analysis.

    Parameters
    ----------
    datadir : string begining with '/'
        directory where the data are stored.


    Returns
    ----------
    None
    '''

    # catdif = "%s/%s_icrf2_dif.sou" % (datadir, datafile[:-4])

    DiffData = sol_icrf2_diff_calc("%s/%s" % (datadir, datafile),
                                   "%s/%s_icrf2_dif.sou" %
                                   (datadir, datafile[:-4]),
                                   "%s_icrf2" % sollab)

    # # IERS transformation parameters.
    # print("# IERS transformation:")
    # # catalog_transfor(catdif, datadir, "%s_icrf2" % sollab)
    # catalog_transfor(DiffData,
    #                  "%s/%s_icrf2_dif.sou" % (datadir, datafile[:-4]),
    #                  label="%s_icrf2" % sollab)
    # # datadir, "%s_icrf2" % sollab)

    # VSH function parameters.
    print("# VSH analysis:")
    # vsh_analysis(catdif, datadir, "%s_icrf2" % sollab)
    vsh_analysis(DiffData,
                 "%s/%s_icrf2_dif.sou" % (datadir, datafile[:-4]),
                 label="%s_icrf2" % sollab)
    # datadir, "%s_icrf2" % sollab)

# Test
# cat = "/home/nliu/solutions/test/a1/result_a1.sou"
# cat = "/home/nliu/solutions/test/GA/opa2017a_aprx.sou"


# #  ------------ Just change this block ----------------
#   | catdir = "/home/nliu/solutions/test/opa2017a_SL"   |
#   | catfil = "opa2017a.sou"                            |
#   | cat = "%s/%s" % (catdir, catfil)                   |
# #  ----------------------------------------------------

# catdif = "%s_icrf2_diff.sou" % cat[:-4]

# sol_icrf2_diff_calc(cat, catdif)

# catalog_transfor(catdif)

# vsh_analysis(catdif)
# --------------------------------- END --------------------------------
