#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: comparison_with_GaiaDR1.py
"""
Created on Mon Feb  5 11:31:11 2018

@author: Neo(liuniu@smail.nju.edu.cn)

N. Liu, 09 Apr 2018: minor changes in functions
                     'solution_Gaia_diff_analysis' and
                     'sol_Gaia_diff_calc'
"""

import numpy as np
import time
from os import path
from functools import reduce
import matplotlib.pyplot as plt
from nor_sep import nor_sep_calc
from VSH_analysis import write_residual, catalog_comparison_VSH


# -----------------------------  FUNCTIONS -----------------------------
def read_gaiadr1(datafile):
    '''Read Gaia DR1 position

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
    dRAc = (RA1 - RA2) * 3.6e9
    dDec = (Dec1 - Dec2) * 3.6e9

    dRAc_err = np.sqrt(RAc_err1**2 + RAc_err2**2)
    dDec_err = np.sqrt(Dec_err1**2 + Dec_err2**2)
    cov = (cor1 * RAc_err1 * Dec_err1 + cor2 * RAc_err2 * Dec_err2)
    cof = cov / dRAc_err / dDec_err

    return [dRAc, dRAc_err, dDec  dDec_err, cov, cof]


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


def comparison_with_Gaia(sou, RA, RAc_err, Dec, Dec_err, cor,
                         catlab, datafile):
    '''Comparison between some catalog with Gaia DR1 qso solution.

    Parameters
    ----------
    sou : array, string
        source name, ICRF designation
    RA / Dec : array, float
        Right ascension / Declination of sources given in degrees.
    RAc_err / Dec_err : array, float
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
    soug, RAg, Decg, RAc_errg, Dec_errg, corg, flg = read_gaiadr1(
        # "/home/nliu/Data/gaiadr1_icrf2_1.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/"
        "gaiadr1_icrf2_1.dat")  # My MacOS

    [soucom, flgcom,
     RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
     RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n] = Xmatch(
        sou, RA, RAc_err, Dec, Dec_err, cor,
        soug, RAg, RAc_errg, Decg, Dec_errg, corg, flg)

    dRAc, dRAc_err, dDec, dDec_err, cov, cof = position_diff_calc(
        [RAn, RAc_errn, Decn, Dec_errn, corn],
        [RAgn, RAc_errgn, Decgn, Dec_errgn, corgn])

    pos_sep, X_a, X_d, X = nor_sep_calc(dRAc, dRAc_err, dDec, dDec_err, cof)

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

    print('###########################################', file=FLOG)
    print('###########################################', file=FTEX)

    # print("# Fitting")
    RdRA1ig, RdDC1ig, RdRA2ig, RdDC2ig = catalog_comparison_VSH(
        datafile,
        soucom, dRAc, dDec, dRAc_err, dDec_err, cov,
        np.deg2rad(RAn), np.deg2rad(Decn), flgcom, FLOG, FTEX,
        pos_sep, X_a, X_d, X)

    # Write the post-fitting residual.
    FLG = (flgcom == 'D')
    # Boolean -> int, 1 for defining source while 0 for non-def
    FLG = FLG.astype(int)
    # the 1st degree
    write_residual(RAdn, DCdn, dRAc_err, dDec_err, RdRA1ig, RdDC1ig, FLG,
                   '%s/logs/%s_vsh.res1' % (main_dir, label))
    # the first two degrees
    write_residual(RAdn, DCdn, dRAc_err, dDec_err, RdRA2ig, RdDC2ig, FLG,
                   '%s/logs/%s_vsh.res2' % (main_dir, label))

    # Plot the position differences
    post_diff_plot(dRAc / 1.e3, dRAc_err / 1.e3,
                   dDec / 1.e3, dDec_err / 1.e3, DCdn,
                   flgcom, main_dir, label)

    print('Done!')

# --------------------------------- END --------------------------------
