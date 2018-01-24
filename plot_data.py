#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:26:18 2017

@author: Neo
"""


from ensure_dir import ensure_dir


# ------------------------------  FUNCTIONS  ---------------------------
def plot_sta_pv(staname, ident, eposta,
                Xsta, Ysta, Zsta, Xsta_err, Ysta_err, Zsta_err, pv='P'):
    '''Plot the Position/Velocity of stations.

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
        unit = 'mm'
    elif pv == 'V' or pv == 'v':
        unit = 'mm/yr'
    else:
        print("pv can only be set one of 'PVpv'")
        exit()

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)
    ax0.errorbar(eposta, Xsta, yerr=Xsta_err, fmt='.')
    ax0.set_title("%s (%s)" % (ident[0], unit))
    # ax0.set_ylim([-500, 500])
    ax1.errorbar(eposta, Ysta, yerr=Ysta_err, fmt='.')
    ax1.set_title("%s (%s)" % (ident[1], unit))
    # ax1.set_ylim([-200, 800])
    ax2.errorbar(eposta, Zsta, yerr=Zsta_err, fmt='.')
    ax2.set_title("%s (%s)" % (ident[2], unit))
    ax2.set_xlabel("Epoch (year)")
    ax2.set_xlim([1979.0, 2018.0])

    figname = "figures/ts_sta/%s_%s.eps" % (staname, ident)
    ensure_dir(figname)
    plt.savefig(figname)
    plt.close()


def plot_sou_pv(souname, eposou,
                RAsou, DCsou, RAsou_err, DCsou_err, pv='P'):
    '''Plot the Position/PM of sources.

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
        unit = 'mas'
        lab = 'sou'
    elif pv == 'V' or pv == 'v':
        unit = 'mas/yr'
        lab = 'pmt'
    else:
        print("pv can only be set one of 'PVpv'")
        exit()

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    ax0.errorbar(eposou, RAsou * 1000, yerr=RAsou_err, fmt='.')
    ax0.set_title("R.A. (%s)" % unit)
    ax1.errorbar(eposou, DCsou * 1000, yerr=DCsou_err, fmt='.')
    ax1.set_title("Dec. (%s)" % unit)
    ax1.set_xlabel("Epoch (year)")
    ax1.set_xlim([1979.0, 2018.0])

    figname = "figures/ts_%s/%s.eps" % (lab, souname)
    ensure_dir(figname)
    plt.savefig(figname)
    plt.close()
# ------------------------------ END -----------------------------------
