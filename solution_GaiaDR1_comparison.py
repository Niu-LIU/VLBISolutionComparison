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

"""

import numpy as np
from os import path
import matplotlib.pyplot as plt
from read_sou import read_sou_pos
from Transformation_cat import cat_transfor
from VSH_analysis import vsh_analysis
cos = np.cos
deg2as = 3.6e3


# -----------------------------  FUNCTIONS -----------------------------
def read_gaiadr1(datafile):
    '''Read Gaia DR1 position

    Note:   1) e_RA = err(RA*cos(DC))
            2) unit for RA/DC : deg
                    for e_RA/e_DC : uas
    '''

    sourcename, Flag = np.genfromtxt(
        datafile, usecols=(1, 13), dtype=str,
        delimiter="|",
        unpack=True)
    RAdeg, DCdeg, e_RAuas, e_DCuas, cor =\
        np.genfromtxt(datafile, usecols=range(3, 8),
                      delimiter="|", unpack=True)

    return sourcename, RAdeg, DCdeg, e_RAuas, e_DCuas, cor, Flag


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

    # as/mas -> uas
    dRA = (RA1 - RA2) * 1.e6 * arccof
    dRA_err = np.sqrt(RA_err1**2 + RA_err2**2)
    dDC = (DC1 - DC2) * 1.e6
    dDC_err = np.sqrt(DC_err1**2 + DC_err2**2)
    cov = (cor1 * RA_err1 * DC_err1 + cor2 * RA_err2 * DC_err2)

    return [dRA, dRA_err, dDC, dDC_err, cov]


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

    plt.legend(fontsize='xx-small')
    plt.savefig("%s/plots/%s_sou_dif.eps" % (main_dir, lab))
    plt.close()


def sol_Gaia_diff_calc(cat, catdif, sollab):
    '''Calculate the quasar position difference between two vlbi-solutions.
    '''

    sou1, RA1, RA_err1, DC1, DC_err1, cor1 = read_sou_pos(cat)
    # err(RA) -> err(RA * cos(DC)), mas -> uas
    arccof = np.cos(np.deg2rad(DC1 / 3600.))
    RA_err1 = RA_err1 * arccof * 1.e3
    DC_err1 = DC_err1 * 1.e3

    sou2, RA2, DC2, RA_err2, DC_err2, cor2, flg = read_gaiadr1(
        # "/home/nliu/Data/gaiadr1_icrf2_1.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/"
        "gaiadr1_icrf2_1.dat")  # My MacOS

    # deg -> as
    deg2as = 3.6e3
    RA2 = RA2 * deg2as
    DC2 = DC2 * deg2as

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
          "# %s \n# GaiaDR1" % cat, file=fdif)
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


def solution_Gaia_diff_analysis(datadir, datafile, sollab):
    '''
    '''

    catdif = "%s/%s_Gaiadr1_dif.sou" % (datadir, datafile[:-4])

    sol_Gaia_diff_calc("%s/%s" % (datadir, datafile), catdif,
                       "%s_GaiaDR1" % sollab)

    # IERS transformation parameter
    print("# IERS transformation:")
    cat_transfor(catdif, datadir, "%s_GaiaDR1" % sollab)

    print("# VSH analysis:")
    # VSH function parameters
    vsh_analysis(catdif, datadir, "%s_GaiaDR1" % sollab)


# solution_Gaia_diff_analysis(
#     "/home/nliu/solutions/GalacticAberration/opa2018a_ga/opa2018a_ga.sou")
# --------------------------------- END --------------------------------
