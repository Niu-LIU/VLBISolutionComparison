#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: comparison_with_GaiaDR1.py
"""
Created on Mon Feb  5 11:31:11 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import time
from os import path
from functools import reduce
import matplotlib.pyplot as plt
from VSH_analysis import write_residual, catalog_comparison_VSH


# -----------------------------  FUNCTIONS -----------------------------
def read_gaiadr1(datafile):
    '''Read Gaia DR1 position.

    NOTICE!!!
    The file being read is not the original data file from CDS.

    Note:   1) e_RA = err(RA*cos(DC))
            2) unit for RA/DC : deg
                    for e_RA/e_DC : uas
    '''

    icrfn, Flag = np.genfromtxt(
        datafile, usecols=(0, 13), dtype=str,
        delimiter="|",
        unpack=True)
    RAdeg, DCdeg, e_RAuas, e_DCuas, cor =\
        np.genfromtxt(datafile, usecols=range(3, 8),
                      delimiter="|", unpack=True)

    return icrfn, RAdeg, DCdeg, e_RAuas, e_DCuas, cor, Flag


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
    '''Crossmatch between catalog and GaiaDR1.
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

    return [soucom, np.array(flgcom),
            RA1n, RA_err1n, DC1n, DC_err1n, cor1n,
            RA2n, RA_err2n, DC2n, DC_err2n, cor2n]


def nor_sep_calc(RA1, DC1, RA1_err, DC1_err, Cor1,
                 RA2, DC2, RA2_err, DC2_err, Cor2,
                 arccof=None):
    '''Calculate the normalized seperation.


    Parameters
    ----------
    RA / DC : Right Ascension / Declination, degress
    e_RA / e_DC : formal uncertainty of RA / DC, micro-as.
    Cor : correlation coeffient between RA and DC.

    Suffix 'G' stands for GaiaDR1 and I for VLBI catalog.

    Returns
    ----------
    X_a / X_d : normalized coordinate differences in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    if arccof is None:
        arccof = np.cos(np.deg2rad(DC1 / 3600.))

    # deg -> uas
    dRA = (RA1 - RA2) * 3.6e9 * arccof
    dRA_err = np.sqrt(RA1_err**2 + RA2_err**2)
    dDC = (DC1 - DC2) * 3.6e9
    dDC_err = np.sqrt(DC1_err**2 + DC2_err**2)

    pos_seps = np.sqrt(dRA**2 + dDC**2)

    # Normalised coordinate difference
    X_a = dRA / np.sqrt(RA1_err**2 + RA2_err**2)
    X_d = dDC / np.sqrt(DC1_err**2 + DC2_err**2)

    # Correlation coefficient of combined errors
    C = (RA1_err * DC1_err * Cor1 +
         RA2_err * DC2_err * Cor2) / (dRA_err * dDC_err)

    # Normalised separation
    X = np.zeros_like(X_a)

    for i, (X_ai, X_di, Ci) in enumerate(zip(X_a, X_d, C)):

        wgt = np.linalg.inv(np.mat([[1, Ci],
                                    [Ci, 1]]))
        Xmat = np.array([X_ai, X_di])
        # Xi = np.sqrt(reduce(np.dot, (Xmat, wgt, Xmat)))
        X[i] = np.sqrt(reduce(np.dot, (Xmat, wgt, Xmat)))

    return X_a, X_d, X


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
            dRA : difference of RA, uas
            dRA_err : formal uncertainty of dRA sqrt(RA_err1^2 + RA_err2^2), uas
            dDC : difference of DC, uas
            dDC_err : formal uncertainty of dDC sqrt(DC_err1^2 + DC_err2^2), uas
            cov : covriance between dRA and dDC, uas ^2
                see Appendix B in Mignard et al 2016
    '''

    RA1, RA_err1, DC1, DC_err1, cor1 = dat1
    RA2, RA_err2, DC2, DC_err2, cor2 = dat2

    arccof = np.cos(np.deg2rad(DC1))

    # deg -> uas
    deg2uas = 3.6e9
    dRA = (RA1 - RA2) * deg2uas * arccof
    dRA_err = np.sqrt(RA_err1**2 + RA_err2**2)
    dDC = (DC1 - DC2) * deg2uas
    dDC_err = np.sqrt(DC_err1**2 + DC_err2**2)
    cov = (cor1 * RA_err1 * DC_err1 + cor2 * RA_err2 * DC_err2)

    return [dRA, dRA_err, dDC, dDC_err, cov]


def sub_posplot(ax0, ax1, dRA, dRA_err, dDC, dDC_err, DC,
                soutp, mk):
    '''
    '''
    ax0.plot(DC, dRA, mk, markersize=1, label=soutp)
    ax1.plot(DC, dDC, mk, markersize=1, label=soutp)


def post_diff_plot(dRA, dRA_err, dDC, dDC_err, DC, tp,
                   main_dir, lab):
    '''Source position difference plot
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


def comparison_with_Gaia(sou, RAd, RA_err, DCd, DC_err, cor,
                         catlab, datafile):
    '''Comparison between some catalog with Gaia DR1 qso solution.

    Parameters
    ----------
    sou : array, string
        source name, ICRF designation
    RA / DC : array, float
        Right ascension / Declination of sources given in degrees.
    RA_err / DC_err : array, float
        formal uncertainty of RA / DC, given in micro-as.
    cor : array, float
        corelation coffecients between RA and DC, unitless.
    catlab : string
        identifier
    datafile : string
        filename with path of the data file


    Returns
    ---------
    '''

    # Gaia DR1 source positions
    soug, RAgd, DCgd, RA_errg, DC_errg, corg, flg = read_gaiadr1(
        # "/home/nliu/Data/gaiadr1_icrf2_1.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/"
        "gaiadr1_icrf2_1.dat")  # My MacOS

    [soucom, flgcom,
     RAdn, RA_errn, DCdn, DC_errn, corn,
     RAgdn, RA_errgn, DCgdn, DC_errgn, corgn] = Xmatch(
        sou, RAd, RA_err, DCd, DC_err, cor,
        soug, RAgd, RA_errg, DCgd, DC_errg, corg, flg)

    dRAu, dRA_err, dDCu, dDC_err, cov = position_diff_calc(
        [RAdn, RA_errn, DCdn, DC_errn, corn],
        [RAgdn, RA_errgn, DCgdn, DC_errgn, corgn])

    X_a, X_d, X = nor_sep_calc(
        RAdn, DCdn, RA_errn, DC_errn, corn,
        RAgdn, DCgdn, RA_errgn, DC_errgn, corgn)

    # VSH analysis
    main_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3"
    label = "%sVGaiaDR1" % catlab
    # Log file.
    FLOG = open("%s/logs/%s_vsh_param.log" % (main_dir, label), "w")
    print('## LOG FILE\n'
          '## Data: %s \n%s' %
          (datafile,
           time.strftime('##%Y-%m-%d %H:%M:%S Begins!',
                         time.localtime(time.time()))),
          file=FLOG)
    # Log file of tex table.
    FTEX = open("%s/logs/%s_vsh_param.tex" % (main_dir, label), "w")
    print('## LOG FILE\n'
          '## The result of different kinds of transformation\n',
          '## Data: %s \n%s' %
          (datafile,
           time.strftime('## %Y-%m-%d %H:%M:%S Begins!',
                         time.localtime(time.time()))),
          file=FTEX)

    ###########################################
    RArn, DCrn = np.deg2rad(RAdn), np.deg2rad(DCdn)  # Unit: rad

    print('###########################################', file=FLOG)
    print('###########################################', file=FTEX)

    # print("# Fitting")
    RdRA1ig, RdDC1ig, RdRA2ig, RdDC2ig = catalog_comparison_VSH(
        datafile,
        soucom, dRAu, dDCu, dRA_err, dDC_err, cov,
        RArn, DCrn, flgcom, FLOG, FTEX, X_a, X_d, X)

    # Write the post-fitting residual.
    FLG = (flgcom == 'D')
    # Boolean -> int, 1 for defining source while 0 for non-def
    FLG = FLG.astype(int)
    # the 1st degree
    write_residual(RAdn, DCdn, dRA_err, dDC_err, RdRA1ig, RdDC1ig, FLG,
                   '%s/logs/%s_vsh.res1' % (main_dir, label))
    # the first two degrees
    write_residual(RAdn, DCdn, dRA_err, dDC_err, RdRA2ig, RdDC2ig, FLG,
                   '%s/logs/%s_vsh.res2' % (main_dir, label))

    # Plot the position differences
    post_diff_plot(dRAu / 1.e3, dRA_err / 1.e3,
                   dDCu / 1.e3, dDC_err / 1.e3, DCdn,
                   flgcom, main_dir, label)

    print('Done!')

# --------------------------------- END --------------------------------
