#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:25:02 2017

@author: Neo
"""

from ensure_dir import ensure_dir


# ------------------------------  FUNCTIONS  ---------------------------
def write_sta_pv(staname, ident, eposta,
                 Xsta, Ysta, Zsta, Xsta_err, Ysta_err, Zsta_err, pv='P'):
    '''write the Position/Velocity of station into a text file.

    Parameters
    ----------
    staname : string
        name of station
    ident : string
        'XYZ' or 'UEN'
    eposta : array, float
        epoch of time series
    Xsta : array, float
        X- / U- component
    Ysta : array, float
        Y- / E- component
    Zsta : array, float
        Z- / N- component
    Xsta_err : array, float
        formal uncertainty of X- / U- component
    Ysta_err : array, float
        formal uncertainty of Y- / E- component
    Zsta_err : array, float
        formal uncertainty of Z- / N- component
    pv : string
        'P' or 'V'

    Returns
    ----------
    None.
    '''

    if pv == 'P' or pv == 'p':
        lab = 'sta'
    elif pv == 'V' or pv == 'v':
        lab = 'vel'
    else:
        print("pv can only be set one of 'PVpv'")
        exit()

    datname = 'ts_%s/%s_%s.dat' % (lab, staname, ident)
    ensure_dir(datname)
    fdat = open(datname, 'w')
    opfmt = "%%10.5f%s%s" % ("|%15.2f" * 3, "|%10.3f" * 3)

    for i, epoi in enumerate(eposta):
        print(opfmt % (epoi,
                       Xsta[i], Ysta[i], Zsta[i],
                       Xsta_err[i], Ysta_err[i], Zsta_err[i]), file=fdat)

    fdat.close()


def write_sou_pv(souname, eposou,
                 RAsou, DCsou, RAsou_err, DCsou_err, corsou, pv='P'):
    '''Write the Position/PM of sources into text files.

    Parameters
    ----------
    souname : string
        IVS source name
    eposou : array, float
        epoch of time series
    RAsou : array, float
        RA component
    DCsou : array, float
        Dec. component
    RAsou_err : array, float
        formal uncertainty of RA component
    DCsou_err : array, float
        formal uncertainty of Dec. component
    pv : string
        'P' or 'V'

    Returns
    ----------
    None.
    '''

    if pv == 'P' or pv == 'p':
        lab = 'sou'
    elif pv == 'V' or pv == 'v':
        lab = 'pmt'
    else:
        print("pv can only be set one of 'PVpv'")
        exit()

    datname = 'ts_%s/%s.dat' % (lab, souname)
    ensure_dir(datname)
    fdat = open(datname, 'w')
    opfmt = "%%10.5f%s%s%s" % ("|%15.2f" * 2, "|%10.3f" * 2, "|%7.4f")

    for i, epoi in enumerate(eposou):
        print(opfmt % (epoi,
                       RAsou[i], DCsou[i],
                       RAsou_err[i], DCsou_err[i], corsou[i]), file=fdat)

    fdat.close()

# ------------------------------ END -----------------------------------
