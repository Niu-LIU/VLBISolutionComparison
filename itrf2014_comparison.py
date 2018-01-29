#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: itrf2014_comparison.py
"""
Created on Thu Jan 25 17:11:19 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from read_itrf import read_itrf, find_stable_sta
from itrf_trans import itrf_trans, trans_fitting
from read_sta import read_sta
from read_vel import read_vel
import time


# -----------------------------  FUNCTIONS -----------------------------
def itrf_sit(datafile):
    '''Read irtf position/velocity of sites.
    '''
    # if label is "p":
    #     datafile = "/home/oper/traitement/itrf2014.sit"
    # else:
    #     datafile = "/home/oper/traitement/itrf2014.vel"

    # datafile = "/home/oper/traitement/itrf2014.sit"
    sta, pvX, pvY, pvZ = read_itrf(datafile)
    stan, pvXn, pvYn, pvZn = find_stable_sta(sta, pvX, pvY, pvZ)

    return stan, pvXn, pvYn, pvZn


def position_taken(index, X, X_err, Y, Y_err, Z, Z_err):
    '''Extract the elements from array at specific index.
    '''

    Xn = np.take(X, index)
    X_errn = np.take(X_err, index)
    Yn = np.take(Y, index)
    Y_errn = np.take(Y_err, index)
    Zn = np.take(Z, index)
    Z_errn = np.take(Z_err, index)

    return Xn, X_errn, Yn, Y_errn, Zn, Z_errn


def Xmatch(sta1, X1, X_err1, Y1, Y_err1, Z1, Z_err1,
           sta2, X2, X_err2, Y2, Y_err2, Z2, Z_err2):
    '''Crossmatch
    '''

    stacom = []
    index1 = []
    index2 = []

    for i, stai in enumerate(sta1):
        indarr = np.where(sta2 == stai)[0]

        # print(indarr)

        if indarr:
            stacom.append(stai)
            index1.append(i)
            j = indarr[0]
            index2.append(j)

    X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n = position_taken(
        index1, X1, X_err1, Y1, Y_err1, Z1, Z_err1)
    X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n = position_taken(
        index2, X2, X_err2, Y2, Y_err2, Z2, Z_err2)

    return [np.array(stacom),
            X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n,
            X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n]


def pv_diff_calc(dat1, dat2):
    '''Calculate the position difference.

    Parameters
    ----------
    dat1, dat2 : list, containing
            X / Y / Z : position (m) / velocity (mm/yr) for X-/Y-/Z-component
            X_err / Y_err / Z_err : formal uncertainty of X / Y / Z,
                                    units are same as X / Y / Z.

    Returns
    ----------
    dif : list, containing:
            dX / dY / dZ : position (m) / velocity (mm/yr) difference
                            for X-/Y-/Z-component
            dX_err / dY_err / dZ_err : formal uncertainty of dX / dY / dZ,
                                    units are same as dX / dY / dZ.
    '''

    X1, X_err1, Y1, Y_err1, Z1, Z_err1 = dat1
    X2, X_err2, Y2, Y_err2, Z2, Z_err2 = dat2

    dX = (X1 - X2)
    dX_err = np.sqrt(X_err1**2 + X_err2**2)
    dY = (Y1 - Y2)
    dY_err = np.sqrt(Y_err1**2 + Y_err2**2)
    dZ = (Z1 - Z2)
    dZ_err = np.sqrt(Z_err1**2 + Z_err2**2)

    return [dX, dX_err, dY, dY_err, dZ, dZ_err]


def itrf14_diff_calc(dat, datdif, label):
    '''Calculate the quasar position difference between two vlbi-solutions.
    '''

    if label is "p" or label is "P":
        [sta, X, X_err, Y, Y_err, Z, Z_err,
         U, U_err, E, E_err, N, N_err,
         ObsUse, ObsTot, SesUse, SesTot, DateBeg, DateEnd,
         XpYp, XpZp, YpZp, XpXv, YpXv, ZpXv,
         XpYv, YpYv, ZpYv, XvYv, XpZv, YpZv,
         ZpZv, XvZv, YvZv] = read_sta(dat)
        datatp = "Position"
        itrffile = "/home/oper/traitement/itrf2014.sit"

    elif label is "v" or label is "V":
        [sta, X, X_err, Y, Y_err, Z, Z_err,
            U, U_err, E, E_err, N, N_err] = read_vel(dat)
        itrffile = "/home/oper/traitement/itrf2014.vel"
        datatp = "Velocity"

    else:
        print(" ERROR is input parameter of itrf14_diff_calc")
        exit()

    sta14, pvX14, pvY14, pvZ14 = read_itrf(itrffile)
    zeroerr = np.zeros_like(pvX14)

    [stacom,
     X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n,
     X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n] = Xmatch(
        sta, X, X_err, Y, Y_err, Z, Z_err,
        sta14, pvX14, zeroerr, pvY14, zeroerr, pvZ14, zeroerr)

    dX, dX_err, dY, dY_err, dZ, dZ_err = pv_diff_calc(
        [X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n],
        [X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n])

    fdif = open(datdif, "w")

    print("# %s difference between two solutions:\n"
          "# %s \n# ITRF2014" % (datatp, dat), file=fdif)
    print("# Sta  X  Y  Z  dX  dXerr  dY  dYerr  dZ  dZerr\n"
          "# Unit: m     for position\n"
          "#       mm/yr for velocity", file=fdif)

    for i, soui in enumerate(soucom):
        print("%9s  %14.10f  %14.10f  %+8.1f  %8.1f  %+8.1f  %8.1f  "
              "%14.1f  %s" %
              (soui, RA1n[i]/3.6e3, DC1n[i]/3.6e3,
               dRA[i], dRA_err[i],
               dDC[i], dDC_err[i],
               cov[i], flg[i]),
              file=fdif)

    fdif.close()

    return stacom, dX, dX_err, dY, dY_err, dZ, dZ_err


def itrf2014_comparison(datadir, stafile, sollab, datatp):
    '''Comparison of site position or velocity between itrf2014 and solution.

    '''

    if datatp is "P" or datatp is "p":
        prefix = "sta"

    elif datatp is "V" or datatp is "v":
        prefix = "vel"

    else:
        print(" ERROR is input parameter of itrf2014_comparison")
        exit()

    stadif = "%s/%s_itrf14_dif.%s" % (datadir, stafile[:-4], prefix)

    stacom, dX, dX_err, dY, dY_err, dZ, dZ_err = itrf14_diff_calc(
        "%s/%s" % (datadir, stafile), stadif, datatp)

    # Log file.
    flog = open("%s/logs/%s_itrf14_param.log" % (datadir, sollab), "w")
    print('## LOG FILE\n'
          '## Data: %s \n%s' %
          (stafile,
           time.strftime('##%Y-%m-%d %H:%M:%S Begins!',
                         time.localtime(time.time()))),
          file=flog)

    # Log file of tex table.
    ftex = open("%s/logs/%s_itrf14_param.tex" % (datadir, sollab), "w")
    print('## LOG FILE\n'
          '## The result of  transformation\n',
          '## Data: %s \n%s' %
          (stafile,
           time.strftime('## %Y-%m-%d %H:%M:%S Begins!',
                         time.localtime(time.time()))),
          file=ftex)

    zerocor = np.zeros_like(dX)

    # ITRF transformation parameter
    print("# ITRF transformation:")
    itrf_trans(sta, pvX, pvY, pvZ,
               pvXerr, pvYerr, pvZerr,
               zerocor, zerocor, zerocor,
               flog, ftex)
    print("Done!")


def test_code():
    '''This is a test code.

    I compare the site position and velocity between itrf08 and itrf14,
    to verify the codes.
    '''

    # compare the position
    # itrf2014
    sta14, pX14, pY14, pZ14 = itrf_sit(
        # "/home/oper/traitement/itrf2014.sit")  # vlbi2
        "/Users/Neo/Astronomy/Data/SOLVE/itrf2014_sitmod")  # My Mac
    sta08, pX08, pY08, pZ08 = itrf_sit(
        # "/obs/nliu/traitement/itrf2008.sit")  # vlbi2
        "/Users/Neo/Astronomy/Data/SOLVE/itrf2008.sit")  # My Mac

    one14 = np.ones_like(pX14)
    one08 = np.ones_like(pX08)

    [stacom,
     X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n,
     X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n] = Xmatch(
        sta08, pX08, one08, pY08, one08, pZ08, one08,
        sta14, pX14, one14, pY14, one14, pZ14, one14)

    dX, dX_err, dY, dY_err, dZ, dZ_err = pv_diff_calc(
        [X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n],
        [X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n])

    zerocor = np.zeros_like(X1n)

    x, sig, cof, ind, _, _, _ = trans_fitting(
        dX, dY, dZ,
        dX_err, dY_err, dZ_err,
        zerocor, zerocor, zerocor,
        X1n, Y1n, Z1n)
    [tx,  ty,  tz,  d, rx,  ry,  rz] = x
    [txerr,  tyerr,  tzerr,  derr, rxerr,  ryerr,  rzerr] = sig

    # Unit transformation
    tx, ty, tz = tx * 1.e3, ty * 1.e3, tz * 1.e3
    txerr, tyerr, tzerr = txerr * 1.e3, tyerr * 1.e3, tzerr * 1.e3
    d, derr = d * 1.e9, derr * 1.e9
    # rad -> mas
    rad2mas = 206265.e3
    rx, ry, rz = rx * rad2mas, ry * rad2mas, rz * rad2mas
    rxerr, ryerr, rzerr = rxerr * rad2mas, ryerr * rad2mas, rzerr * rad2mas

    print("===============================================")
    print("Position comparison:")
    print("#### Translation component:\n",
          " %+8.3f +/- %8.3f |" * 3 % (tx, txerr, ty, tyerr, tz, tzerr))
    print("#### Scale factor:\n",
          " %+8.3f +/- %8.3f" % (d, derr))
    print("#### Rotation component:\n",
          " %+8.3f +/- %8.3f |" * 3 % (rx, rxerr, ry, ryerr, rz, rzerr))
    print("##   correlation coefficients are:\n", cof)
    print("===============================================")

    # compare the velocity
    # itrf2014
    sta14, pX14, pY14, pZ14 = itrf_sit(
        # "/home/oper/traitement/itrf2014.vel")  # vlbi2
        "/Users/Neo/Astronomy/Data/SOLVE/itrf2014_velmod")  # My Mac
    sta08, pX08, pY08, pZ08 = itrf_sit(
        # "/obs/nliu/traitement/itrf2008.vel") # vlbi2
        "/Users/Neo/Astronomy/Data/SOLVE/itrf2014_velmod")  # My Mac
    one14 = np.ones_like(pX14)
    one08 = np.ones_like(pX08)

    [stacom,
     X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n,
     X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n] = Xmatch(
        sta08, pX08, one08, pY08, one08, pZ08, one08,
        sta14, pX14, one14, pY14, one14, pZ14, one14)

    dX, dX_err, dY, dY_err, dZ, dZ_err = pv_diff_calc(
        [X1n, X_err1n, Y1n, Y_err1n, Z1n, Z_err1n],
        [X2n, X_err2n, Y2n, Y_err2n, Z2n, Z_err2n])

    zerocor = np.zeros_like(X1n)

    x, sig, cof, ind, _, _, _ = trans_fitting(
        dX, dY, dZ,
        dX_err, dY_err, dZ_err,
        zerocor, zerocor, zerocor,
        X1n, Y1n, Z1n)
    # print(X1n.size)
    [tx,  ty,  tz,  d, rx,  ry,  rz] = x
    [txerr,  tyerr,  tzerr,  derr, rxerr,  ryerr,  rzerr] = sig

    d, derr = d * 1.e9, derr * 1.e9
    # rad -> mas
    rad2mas = 206265.e3
    rx, ry, rz = rx * rad2mas, ry * rad2mas, rz * rad2mas
    rxerr, ryerr, rzerr = rxerr * rad2mas, ryerr * rad2mas, rzerr * rad2mas

    print("===============================================")
    print("Velocity comparison:")
    print("#### Translation component:\n",
          " %+8.3f +/- %8.3f |" * 3 % (tx, txerr, ty, tyerr, tz, tzerr))
    print("#### Scale factor:\n",
          " %+8.3f +/- %8.3f" % (d, derr))
    print("#### Rotation component:\n",
          " %+8.3f +/- %8.3f |" * 3 % (rx, rxerr, ry, ryerr, rz, rzerr))
    print("##   correlation coefficients are:\n", cof)
    print("===============================================")


# --------------------------------- MAIN --------------------------------
test_code()
# --------------------------------- END --------------------------------
