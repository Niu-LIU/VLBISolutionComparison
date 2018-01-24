#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: write_solvesrc.py
"""
Created on Mon Jan 22 15:04:40 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np


# -----------------------------  FUNCTIONS -----------------------------
def deg2dms(angle):
    '''RA(deg) -> D, M, S
    '''

    sgn = np.sign(angle)
    angle = np.fabs(angle)

    D = angle // 1
    M = angle % 1 * 60 // 1
    S = angle % 1 * 60 % 1 * 60

    return sgn, D, M, S

# Test codes
# print(deg2dms(-32.3503429220))
# print(deg2dms(0.0849978345 / 15))


# def write_solvesrc(ivs, RA, DC, errDC):
def write_solvesrc(ivs, RA, DC, errDC, comment, fout):
    '''Write qso position used for SOLVE.

    Parameters
    ----------
    ivs : array, string
        IVS designation name of radio source
    RA : array, float
        Right ascension, degree
    DC : array, float
        Declination, degree
    errDC : array, float
        Formal uncertainty of DC, micro-arcsec
    comment : string
        comment
    fout :
        Output file

    Returns
    ----------
    None
    '''

    _, RAh, RAm, RAs = deg2dms(RA / 15.0)     # arc-sec -> second
    DCsgn, DCd, DCm, DCs = deg2dms(DC)

    sgn = np.where(DCsgn == -1, '-', ' ')

    # # Output file
    # outfile = "/Users/Neo/Astronomy/Data/catalogues/Gaia_cds/qso_solve.dat"

    # Loop for writing data
    linefmt = "    %8s  %2d %2d %11.8f   %s%02d %2d  %10.7f   %7.3f   %s"

    # for i, souname in enumerate(ivs):
    # if DCsgn[i] == -1:
    #     sgn = '-'
    # else:
    #     sgn = ' '
    # print(linefmt % (souname,
    #                  RAh[i], RAm[i], RAs[i],
    #                  sgn[i], DCd[i], DCm[i], DCs[i], errDC[i]),
    #       file=fout)

    for (souni, RAhi, RAmi, RAsi, sgni, DCdi, DCmi, DCsi,
         errDCi) in zip(ivs, RAh, RAm, RAs, sgn, DCd, DCm, DCs, errDC):

        print(linefmt % (souni, RAhi, RAmi, RAsi,
                         sgni, DCdi, DCmi, DCsi, errDCi, comment),
              file=fout)


def write_NNRS(soun, flst):
    '''Write the No-Net-Rotation source list.

    Parameters
    ----------
    soun: array, string
        IVS designation name of radio source
    flst: Output file

    Returns
    ----------
    None
    '''

    N = 8  # 8 source pre line
    # m = soun.size % N

    # linefmt = "     %s \\" % ("%-8s " * N)
    # for i in range(m + 1):
    #     [s0, s1, s2, s3, s4, s5, s6, s7] = soun[i * N: (i+1) * N]
    #     print(linefmt % (s0, s1, s2, s3, s4, s5, s6, s7))

    for i, souni in enumerate(soun):
        if not i % N:
            if i:
                print("\\\n     ", end="", file=flst)
            else:
                print("     ", end="", file=flst)

        print("%-9s" % souni, end="", file=flst)

    # remnfmt = "     %s" % ("%-8s " * (soun.size % N))
    # print(remnfmt % ((soun[i * m:])))

# --------------------------------- END -------------------------------
