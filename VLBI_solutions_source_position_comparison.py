#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: VLBI_solutions_comparison.py
"""
Created on Fri Dec 22 14:09:36 2017

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import matplotlib.pyplot as plt
from read_sou import read_sou_pos
from Transformation_cat import cat_transfor
from VSH_analysis import vsh_analysis


# -----------------------------  FUNCTIONS -----------------------------
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
           sou2, RA2, RA_err2, DC2, DC_err2, cor2):
    '''Crossmatch
    '''

    soucom = []
    index1 = []
    index2 = []

    for i, soui in enumerate(sou1):
        indarr = np.where(sou2 == soui)[0]

        if indarr is not np.array([]):
            soucom.append(soui)
            index1.append(i)
            index2.append(indarr[0])

    RA1n, RA_err1n, DC1n, DC_err1n, cor1n = position_taken(
        index1, RA1, RA_err1, DC1, DC_err1, cor1)
    RA2n, RA_err2n, DC2n, DC_err2n, cor2n = position_taken(
        index2, RA2, RA_err2, DC2, DC_err2, cor2)

    return [soucom,
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

    # cond = arccof == 0
    # print(arccof[cond])

    dRA = (RA1 - RA2) * 1.e6 * arccof
    dRA_err = np.sqrt(RA_err1**2 + RA_err2**2) * 1.e3 * np.fabs(arccof)
    dDC = (DC1 - DC2) * 1.e6
    dDC_err = np.sqrt(DC_err1**2 + DC_err2**2) * 1.e3
    cov = (cor1 * RA_err1 * DC_err1 + cor2 * RA_err2 * DC_err2) * 1.e6

    return [dRA, dRA_err, dDC, dDC_err, cov]


def sub_posplot(ax0, ax1, dRA, dRA_err, dDC, dDC_err, DC, soutp):
    '''
    '''
    ax0.plot(DC, dRA, '.', markersize=1, label=soutp)
    ax1.plot(DC, dDC, '.', markersize=1, label=soutp)


def post_diff_plot(dRA, dRA_err, dDC, dDC_err, DC, tp):
    '''
    '''

    unit = "\mu as"
    lim = 50
    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    # ax0.errorbar(DC, dRA, yerr=dRA_err, fmt='.', elinewidth=0.01,
    #              markersize=1)
    ax0.plot(DC, dRA, '.', markersize=1)
    ax0.plot(DC, dDC, '.', markersize=1)

    flgs = ['D', 'V', 'N', 'O']
    soutps = ['Defining', 'VCS', 'Non-VCS', 'Non-icrf2']
    for flg, soutp in zip(flgs, soutps):
        cond = tp == flg
        sub_posplot(ax0, ax1,
                    dRA[cond], dRA_err[cond],
                    dDC[cond], dDC_err[cond],
                    DC[cond], soutp)

    ax0.set_ylabel("$\Delta R.A. (%s)$" % unit)
    ax0.set_ylim([-lim, lim])
    # ax1.errorbar(DC, dDC, yerr=dDC_err, fmt='.', elinewidth=0.01,
    #              markersize=1)
    ax1.set_ylabel("$\Delta Dec. (%s)$" % unit)
    ax1.set_ylim([-lim, lim])
    ax1.set_xlabel("Dec. (deg)")
    ax1.set_xlim([-90, 90])
    plt.legend(fontsize='xx-small')
    plt.savefig("/home/nliu/solutions/test/GA/sou_post_diff.eps")  # vlbi2
    # plt.savefig("/home/nliu/solutions/test/GA/sou_post_diff.eps")
    plt.close()


def catlg_diff_calc(cat1, cat2, catdif):
    '''Calculate the quasar position difference between two vlbi-solutions.
    '''

    [sou1, RA1, RA_err1, DC1, DC_err1, cor1] = read_sou_pos(cat1)
    [sou2, RA2, RA_err2, DC2, DC_err2, cor2] = read_sou_pos(cat2)

    [soucom,
     RA1n, RA_err1n, DC1n, DC_err1n, cor1n,
     RA2n, RA_err2n, DC2n, DC_err2n, cor2n] = Xmatch(
        sou1, RA1, RA_err1, DC1, DC_err1, cor1,
        sou2, RA2, RA_err2, DC2, DC_err2, cor2)

    dRA, dRA_err, dDC, dDC_err, cov = position_diff_calc(
        [RA1n, RA_err1n, DC1n, DC_err1n, cor1n],
        [RA2n, RA_err2n, DC2n, DC_err2n, cor2n])

    souicrf2, flg = np.genfromtxt(
        "/home/nliu/Data/icrf2.dat",
        usecols=(1, 3), dtype=str, unpack=True)

    fdif = open(catdif, "w")

    print("# Position difference between two catalogs:\n"
          "# %s\n# %s" % (cat1, cat2), file=fdif)

    RA1n = RA1n / 3.6e3
    DC1n = DC1n / 3.6e3

    tp = []

    for i, soui in enumerate(soucom):

        indarr = np.where(souicrf2 == soui)[0]

        if indarr:
            j = indarr[0]

            print("%9s  %14.10f  %14.10f  %+8.1f  %8.1f  %+8.1f  %8.1f  "
                  "%14.1f  %s" %
                  (soui, RA1n[i], DC1n[i],
                   dRA[i], dRA_err[i],
                   dDC[i], dDC_err[i],
                   cov[i], flg[j]),
                  file=fdif)

            tp.append(flg[j])

        else:

            print("%9s  %14.10f  %14.10f  %+8.1f  %8.1f  %+8.1f  %8.1f"
                  "  %14.1f  O" %
                  (soui, RA1n[i], DC1n[i],
                   dRA[i], dRA_err[i], dDC[i], dDC_err[i], cov[i]),
                  file=fdif)
            # O for non-icrf2 source

            tp.append('O')

    post_diff_plot(dRA, dRA_err, dDC, dDC_err, DC1n, np.array(tp))

    fdif.close()

# -------------------- MAIN ----------------------------------

# # Test
# cat1 = "/home/nliu/solutions/test/a1/result_a1.sou"
# cat2 = "/home/nliu/solutions/test/a2/result_a2.sou"
# catdif = "/home/nliu/solutions/test/a1_a2_dif.sou"
# catlg_diff_calc(cat1, cat2, catdif)
# # cat_transfor(catdif)
# vsh_analysis(catdif)

# GA - Non-GA
cat1 = "/home/nliu/solutions/test/GA/opa2017a_aprx.sou"
cat2 = "/home/nliu/solutions/test/GA/opa2017a_ga.sou"
catdif = "/home/nliu/solutions/test/GA/opa2017a_ga_dif.sou"

catlg_diff_calc(cat1, cat2, catdif)

cat_transfor(catdif)

vsh_analysis(catdif)
# --------------------------------- END -------------------------------
