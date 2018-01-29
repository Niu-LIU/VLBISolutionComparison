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

"""

import numpy as np
import matplotlib.pyplot as plt
from os import path
from read_sou import read_sou_pos
from Transformation_cat import cat_transfor
from VSH_analysis import vsh_analysis
cos = np.cos
deg2as = 3.6e3


# -----------------------------  FUNCTIONS -----------------------------
def read_icrf2(datafile):
    sourcename, Flag = np.genfromtxt(
        datafile, usecols=(1, 3), dtype=str, unpack=True)
    RAh, RAm, RAs, Decd, Decm, Decs, e_RAs, e_DEas, cor = np.genfromtxt(
        datafile, usecols=range(4, 13), unpack=True)

# determine the sign of Declination
    strDecs = np.genfromtxt(datafile, usecols=(7,), dtype=str)
    # Method 01
    # Dec_sign = np.ones_like(Decd)
    # for i, strDec in enumerate(strDecs):
    #     if strDec[0] == '-':
    #         Dec_sign[i] = -1
    # Method 02
    signstr = [x[0] for x in strDecs]
    Dec_sign = np.where(np.array(signstr) == '-', -1, 1)

# calculate the position
    RA = (RAh + RAm / 60.0 + RAs / 3600) * 15 * deg2as  # degree -> as
    Dec = Decd + Dec_sign * (Decm / 60.0 + Decs / 3600)  # degree
    Dec = Dec * deg2as  # degree ->as
# unit: s/as -> mas
    e_RA = e_RAs * 15e3
    e_DE = e_DEas * 1.0e3
    return sourcename, RA, Dec, e_RA, e_DE, cor, Flag


def position_taken(index, RA, RA_err, DC, DC_err, cor):
    '''Extract the elements from array at specific index.
    '''

    RAn = np.take(RA, index)
    RA_errn = np.take(RA_err, index)
    DCn = np.take(DC, index)
    DC_errn = np.take(DC_err, index)
    corn = np.take(cor, index)

    return RAn, RA_errn, DCn, DC_errn, corn


def Xmatch(sou1, RA1, RA_err1, DC1, DC_err1, cor1,
           sou2, RA2, RA_err2, DC2, DC_err2, cor2, flg):
    '''Crossmatch
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

    # print(len(index1), len(index2))
    # print(index1)
    # print(index2)

    RA1n, RA_err1n, DC1n, DC_err1n, cor1n = position_taken(
        index1, RA1, RA_err1, DC1, DC_err1, cor1)
    RA2n, RA_err2n, DC2n, DC_err2n, cor2n = position_taken(
        index2, RA2, RA_err2, DC2, DC_err2, cor2)

    # print(RA1n.size, RA2n.size)

    return [soucom, flgcom,
            RA1n, RA_err1n, DC1n, DC_err1n, cor1n,
            RA2n, RA_err2n, DC2n, DC_err2n, cor2n]


def position_diff_calc(dat1, dat2):
    '''Calculate the position difference.

    Parameters
    ----------
    dat1, dat2 : list, containing
            RA : Right ascension, arc-sec
            RA_err : formal uncertainty of RA, mas
            DC : declination, arc-sec
            DC_err : formal uncertainty of DC, mas
            cor : correlation coefficient between RA and DC

    Returns
    ----------
    dif : list, containing:
            dRA : difference of RA, uas
            dRA_err : formal uncertainty of dRA sqrt(RA_err1^2 + RA_err2^2), uas
            dDC : difference of DC, uas
            dDC_err : formal uncertainty of dDC sqrt(DC_err1^2 + DC_err2^2), uas
            cov : covriance between dRA and dDC, uas ^2
                see Appendix B in Mignard et al 2016
    '''

    RA1, RA_err1, DC1, DC_err1, cor1 = dat1
    RA2, RA_err2, DC2, DC_err2, cor2 = dat2

    arccof = np.cos(np.deg2rad(DC1 / 3600.))

    dRA = (RA1 - RA2) * 1.e6 * arccof
    dRA_err = np.sqrt(RA_err1**2 + RA_err2**2) * 1.e3 * np.fabs(arccof)
    dDC = (DC1 - DC2) * 1.e6
    dDC_err = np.sqrt(DC_err1**2 + DC_err2**2) * 1.e3
    cov = (cor1 * RA_err1 * DC_err1 + cor2 * RA_err2 * DC_err2) * 1.e6

    return [dRA, dRA_err, dDC, dDC_err, cov]


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

    ax0.set_ylabel("$\Delta \\alpha^* (%s)$" % unit)
    ax0.set_ylim([-lim, lim])

    # ax1.errorbar(DC, dDC, yerr=dDC_err, fmt='.', elinewidth=0.01,
    #              markersize=1)
    ax1.set_ylabel("$\Delta \delta (%s)$" % unit)
    ax1.set_ylim([-lim, lim])
    ax1.set_xlabel("Dec. (deg)")
    ax1.set_xlim([-90, 90])

    plt.legend(fontsize='xx-small')
    plt.savefig("%s/plots/%s_sou_dif.eps" % (main_dir, sollab))
    plt.close()


def sol_icrf2_diff_calc(cat, catdif, sollab):
    '''Calculate the quasar position difference between two vlbi-solutions.
    '''

    sou1, RA1, RA_err1, DC1, DC_err1, cor1 = read_sou_pos(cat)
    sou2, RA2, DC2, RA_err2, DC_err2, cor2, flg = read_icrf2(
        # "/home/nliu/Data/icrf2.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat")  # My MacOS

    [soucom, flgcom,
     RA1n, RA_err1n, DC1n, DC_err1n, cor1n,
     RA2n, RA_err2n, DC2n, DC_err2n, cor2n] = Xmatch(
        sou1, RA1, RA_err1, DC1, DC_err1, cor1,
        sou2, RA2, RA_err2, DC2, DC_err2, cor2, flg)

    dRA, dRA_err, dDC, dDC_err, cov = position_diff_calc(
        [RA1n, RA_err1n, DC1n, DC_err1n, cor1n],
        [RA2n, RA_err2n, DC2n, DC_err2n, cor2n])

    fdif = open(catdif, "w")

    print("# Position difference between two catalogs:\n"
          "# %s \n# ICRF2" % cat, file=fdif)
    print("# Sou  RA  DC  dRA  dRAerr  dDC  dDCerr  COV  flag\n"
          "#      deg deg uas  uas     uas  uas     uas^2\n"
          "# For flag, 'D', 'V', and 'N' are same as in ICRF2, "
          "'O' for non-icrf2 source", file=fdif)

    for i, soui in enumerate(soucom):
        print("%9s  %14.10f  %14.10f  %+8.1f  %8.1f  %+8.1f  %8.1f  "
              "%14.1f  %s" %
              (soui, RA1n[i]/3.6e3, DC1n[i]/3.6e3,
               dRA[i], dRA_err[i],
               dDC[i], dDC_err[i],
               cov[i], flgcom[i]),
              file=fdif)

    fdif.close()

    # main_dir = "/home/nliu/solutions/GalacticAberration" # vlbi2
    # main_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3"  # My MacOS
    main_dir = path.dirname(cat)
    post_diff_plot(dRA / 1.e3, dRA_err / 1.e3,
                   dDC / 1.e3, dDC_err / 1.e3, DC1n / 3.6e3,
                   np.array(flgcom), main_dir, sollab)


def solution_icrf2_diff_analysis(datadir, datafile, sollab):
    '''
    '''

    catdif = "%s/%s_icrf2_dif.sou" % (datadir, datafile[:-4])

    sol_icrf2_diff_calc("%s/%s" % (datadir, datafile), catdif,
                        "%s_icrf2" % sollab)

    # IERS transformation parameter
    print("# IERS transformation:")
    cat_transfor(catdif, datadir, "%s_icrf2" % sollab)

    print("# VSH analysis:")
    # VSH function parameters
    vsh_analysis(catdif, datadir, "%s_icrf2" % sollab)

# Test
# cat = "/home/nliu/solutions/test/a1/result_a1.sou"
# cat = "/home/nliu/solutions/test/GA/opa2017a_aprx.sou"


# # ---------- Just change this block ----------------
# catdir = "/home/nliu/solutions/test/opa2017a_SL"  # |
# catfil = "opa2017a.sou"                           # |
# cat = "%s/%s" % (catdir, catfil)                  # |
# # --------------------------------------------------

# catdif = "%s_icrf2_diff.sou" % cat[:-4]

# sol_icrf2_diff_calc(cat, catdif)

# cat_transfor(catdif)

# vsh_analysis(catdif)
# --------------------------------- END --------------------------------
