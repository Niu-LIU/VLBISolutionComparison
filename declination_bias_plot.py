#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: declination_bias_plot.py
"""
Created on Tue Mar 13 11:01:41 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Plot the declination differences for ICRF2 defining sources
to evince the declination bias or regional error in VLBI, or more
specifically ICRF, system.

"""

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------  FUNCTIONS -----------------------------
def dec_bias_plot(fname, Dec, DDec, EDDec, ymin=-500, ymax=500):
    '''Plot the declination differences.

    Parameters
    ----------
        fname : string
            output figure name, including the full path;
        Dec : array, float
            Declinations in degree;
        DDec : array of float
            Declination differences in micro-arcsecond;
        EDDec : array of float
            Formal uncertainties of declination difference in micro-as;
        ymin/ymax : array of float
            minimum/maximum of the y-axis.

    Returns
    ----------
        None
    '''

    plt.figure(figsize=(10, 6))
    # plt.text(-90, 0.51, '$\\times 10^3$', fontsize=12)
    plt.hlines(y=0, xmin=-90, xmax=90, linestyles='dashed', lw=0.5)
    plt.errorbar(Dec, DDec, yerr=EDDec, fmt='b.', elinewidth=0.025,
                 markersize=3)
    plt.xlabel('Dec.(deg)', fontsize=18)
    plt.ylabel('$\Delta\delta$ ($\mu$as)', fontsize=18)
    plt.xticks(np.arange(-90, 91, 30), fontsize=18)
    # plt.yticks(np.arange(-0.4, 0.5, 0.2),
    # ('$-400$', '$-200$', '0', '200', '400'), fontsize=20)
    # ('$-400$', '$-200$', '0', '200', '400')
    plt.xlim([-90, 90])
    # plt.ylim([-500, 500])
    plt.ylim([ymin, ymax])
    plt.savefig(fname)
    plt.close()


# ------------------------------------------------------
def read_icrf2():
    '''Read ICRF2 file.

    Parameters
    ----------
    None

    Returns
    ----------
    ivsn : array of string
        IVS designation of source name;
    Dec : array of float
        Declination in degree;
    e_Dec : array of float
        formal error of RA/Dec in micro-arcsecond;
    Flag : array of character
        flag of source classification in ICRF2 catalog.
    '''

    icrf2_fil = "/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat"

    ivsn, Flag = np.genfromtxt(icrf2_fil,
                               usecols=(1, 3), dtype=str, unpack=True)
    Decd, Decm, Decs, e_Dec_as = np.genfromtxt(icrf2_fil,
                                               usecols=(7, 8, 9, 11),
                                               unpack=True)

# determine the sign of Declination
    strDecs = np.genfromtxt(icrf2_fil, usecols=(7,), dtype=str)
    signstr = [x[0] for x in strDecs]
    Dec_sign = np.where(np.array(signstr) == '-', -1, 1)

# calculate the position
    Dec = Decd + Dec_sign * (Decm / 60.0 + Decs / 3600)  # degree

# unit: as -> uas
    e_Dec = e_Dec_as * 1.0e6

    return ivsn, Dec, e_Dec, Flag


# ------------------------------------------------------
def get_icrf2_def():
    '''Get ICRF2 defining source position.

    Parameters
    ----------
    None

    Returns
    ----------
    ivsnD : array of string
        IVS designation of source name;
    DecD : array of float
        Declination in degree;
    e_DecD : array of float
        formal error of RA/Dec in micro-arcsecond;
    '''

    # fetch the whole ICRF2 data
    ivsn, Dec, e_Dec, Flag = read_icrf2()

    # find ICRF2 defining sources
    index = (Flag == 'D')

    # Extract the data for subset of ICRF2 defining.
    ivsnD = ivsn[index]
    DecD = Dec[index]
    e_DecD = e_Dec[index]

    # To verify if we get the correct data
    if ivsnD.size == 295:
        return ivsnD, DecD, e_DecD
    else:
        print("Error in program get_icrf2_def")
        exit()


# ------------------------------------------------------
def Xmatch(sou1, DC1, DC_err1, sou2, DC2, DC_err2):
    '''Crossmatch between two catalogs.

    Parameters
    ----------
    sou1/sou2 : array of string
        Source name
    DC1/DC2 : array of float
        Declination given in degree
    DC_err1/DC_err2 : array of float
        Formal uncertainty of declination in micro-arcsecond

    Returns
    ----------
    soucom : array of string
        Common source name
    DC1com/DC2com : array of float
        Declinations for common sources in degree
    DC_err1com/DC_err2com : array of float
        Formal uncertainty of declination in micro-arcsecond

    '''

    soucom = []
    index1 = []
    index2 = []

    for i, soui in enumerate(sou1):
        indarr = np.where(sou2 == soui)[0]

        if indarr:
            soucom.append(soui)
            index1.append(i)
            j = indarr[0]
            index2.append(j)

    # Retieve the data for common sources
    soucom = np.array(soucom)
    DC1com = np.take(DC1, index1)
    DC2com = np.take(DC2, index2)
    DC_err1com = np.take(DC_err1, index1)
    DC_err2com = np.take(DC_err2, index2)

    return [soucom,
            DC1com, DC_err1com,
            DC2com, DC_err2com]


# ------------------------------------------------------
def position_diff_calc(dat1, dat2):
    '''Calculate the declination difference.

    Parameters
    ----------
    dat1, dat2 : list, containing
            sou : source name
            DC : declination, degree
            DC_err : formal uncertainty of DC, mirco-as

    Returns
    ----------
    dif : list, containing:
            soucom : common source names
            DCcom : declination of common sources
            dDC : difference of DC, uas
            dDC_err : formal uncertainty of dDC sqrt(DC_err1^2 + DC_err2^2), uas
    '''

    sou1, DC1, DC_err1 = dat1
    sou2, DC2, DC_err2 = dat2

    # Cross-match
    soucom, DC1com, DC_err1com, DC2com, DC_err2com = Xmatch(
        sou1, DC1, DC_err1, sou2, DC2, DC_err2)

    # Degree -> micro-as
    dDC = (DC1com - DC2com) * 3.6e9
    dDC_err = np.sqrt(DC_err1com**2 + DC_err2com**2)

    return [soucom, DC2com, dDC, dDC_err]


# ------------------------------------------------------
def dec_bias_wrtICRF2(sou, DC, DC_err, label):
    '''Plot the declination difference wrt ICRF2 catalogs via defining sources.

    Take 'sol - ICRF2' in sense.

    Parameters
    ----------
    sou : array of string
        source name
    DC : array of float
        declination in degree
    DC_err : array of float
        formal error of DC

    Returns
    ----------
    None
    '''

    # Get the declination of 295 ICRF2 defining sources.
    ivsnD, DecD, e_DecD = get_icrf2_def()

    # Calculate the difference
    soucom, DCcom, dDC, dDC_err = position_diff_calc(
        [sou, DC, DC_err], [ivsnD, DecD, e_DecD])

    # Plot the difference
    dec_bias_plot('/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/'
                  '%s_icrf2_DecDif.eps' % label, DCcom, dDC, dDC_err,
                  ymin=-200, ymax=200)


# --------------------------------- END --------------------------------
