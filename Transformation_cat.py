#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: Transformation_cat.py
"""
Created on Wed Jan  3 18:08:42 2018

@author: Neo(liuniu@smail.nju.edu.cn)

History
N. Liu, 10 Feb 2018: change the input parameters of function
                     'catalog_transfor', replacing variable 'datafile'
                     with 'DiffData'
"""

import numpy as np
import time
from os import path
# from post_trans import tran01_01, tran01_02, \
#     tran01_03, tran02_fitting, tran03_fitting
from post_trans import tran02_fitting, tran03_fitting
from tex_table import write_result_deg1


# -----------------------------  FUNCTIONS -----------------------------
def mod_calc(x):
    return np.sqrt(np.dot(x, x))


#------------------------------------------
def sig_calc(x, s):
    return np.sqrt(np.dot(x**2, s**2) / np.dot(x, x))


#------------------------------------------
# def Orientation_calc(Condition, RA, DE, dRA, dDC, dRA_err, dDC_err, Flog):
def tran_fitting(Condition, RA, DC, dRA, dDC,
                 dRA_err, dDC_err, COV, Flog, Ftex):

    ra = np.extract(Condition, RA)
    dc = np.extract(Condition, DC)
    dra = np.extract(Condition, dRA)
    ddc = np.extract(Condition, dDC)
    dra_err = np.extract(Condition, dRA_err)
    ddc_err = np.extract(Condition, dDC_err)
    cov = np.extract(Condition, COV)

    print('## Number of the Sample: %d' % ra.size, file=Flog)

# Only rigid rotation => tran03_fitting.
    w, sig, cormat, _ = tran03_fitting(dra, ddc, dra_err, ddc_err,
                                       cov, ra, dc)

    [wx,  wy,  wz] = w
    [ewx, ewy, ewz] = sig
    wtol, ewtol = mod_calc(w), sig_calc(w, sig)

    # Write result into log file.
    print('###### Rigid rotation only.\n'
          '## The orientation are(uas):\n',
          ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
          '=> %4d +/- %3d' % (wtol, ewtol), file=Flog)
    print('##    correlation coefficients are:\n', cormat, file=Flog)
    # For tex table.
    print('###### Rigid rotation only.', file=Ftex)
    print('## (r_x, r_y, r_z, r)', file=Ftex)
    print('  &$%+4d \pm %3d$' * 3 % (wx, ewx, wy, ewy, wz, ewz) +
          '  &$%4d \pm %3d$' % (wtol, ewtol), file=Ftex)
    write_result_deg1(["A_1", "A_2", "A_3"], w, sig, cormat, Ftex)

# consider a possible bias in declination => tran02_fitting.
    w, sig, cormat, _ = tran02_fitting(dra, ddc, dra_err, ddc_err,
                                       cov, ra, dc)
    [wx, wy, z, b] = w
    [ewx, ewy, ewz, eb] = sig
    wtol, ewtol = mod_calc(w[:3]), sig_calc(w[:3], sig[:3])
    # Write result into log file.
    print('###### Consider a bias in declination.', file=Flog)
    print('## The orientation are(uas):\n',
          ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
          '=> %4d +/- %3d' % (wtol, ewtol), file=Flog)
    print('## The bias in declination(uas):\n',
          ' %+4d +/- %3d' % (b, eb), file=Flog)
    print('##   correlation coefficients are:\n', cormat, file=Flog)
    # For tex table.
    print('###### Consider a bias in declination.', file=Ftex)
    print('## (r_x, r_y, r_z, b, r)', file=Ftex)
    print('  &$%+4d \pm %3d$' * 4 % (wx, ewx, wy, ewy, wz, ewz, b, eb) +
          '  &$%4d \pm %3d$' % (wtol, ewtol), file=Ftex)
    write_result_deg1(["A_1", "A_2", "A_3", "dz"], w, sig, cormat, Ftex)


def cat_comparison(tag, RA, DC, dRA, dDC, dRA_err, dDC_err, COV,
                   flg, FLOG, FTEX):

    print('##============================================\n'
          '## %s' % tag, file=FLOG)
    print('  # ============================================\n'
          '## %s' % tag, file=FTEX)

    # All sources
    con = (DC <= np.pi / 2)
    print('##--------- For All sources:', file=FLOG)
    print('##--------- For All sources:', file=FTEX)
    tran_fitting(con, RA, DC, dRA, dDC, dRA_err, dDC_err, COV,
                 FLOG, FTEX)

    # Divide data into 2 sets: Nouthern and Southern declinations.
    # For Nouthern hemisphere
    # con = (DC > 0)  # Nouthern
    # print('##--------- For Northern hemisphere:', file=FLOG)
    # print('##--------- For Northern hemisphere:', file=FTEX)
    # tran_fitting(con, RA, DC, dRA, dDC, dRA_err, dDC_err, COV, FLOG, FTEX)

    # # For Southern hemisphere
    # con = (DC < 0)  # Southern
    # print('##--------- For Southern hemisphere:', file=FLOG)
    # print('##--------- For Southern hemisphere:', file=FTEX)
    # tran_fitting(con, RA, DC, dRA, dDC, dRA_err, dDC_err, COV, FLOG, FTEX)

    # For Defining sources.
    con = (flg == 'D')
    print('##--------- For Defining sources:', file=FLOG)
    print('##--------- For Defining sources:', file=FTEX)
    tran_fitting(con, RA, DC, dRA, dDC, dRA_err, dDC_err, COV, FLOG, FTEX)

    # For Non-defining sources.
    con = (flg != 'D')
    print('##--------- For Non-defining sources:', file=FLOG)
    print('##--------- For Non-defining sources:', file=FTEX)
    tran_fitting(con, RA, DC, dRA, dDC, dRA_err, dDC_err, COV, FLOG, FTEX)


def catalog_transfor(DiffData, datafile,
                     # vlbi2
                     # main_dir="/home/nliu/solutions/GalacticAberration",
                     # My MacOS
                     main_dir="/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3",
                     label=''):
    '''Estiamte the coordinate transformation parameters from positional offset.
    '''

    # main_dir = "/home/nliu/solutions/GalacticAberration"

    # Log file.
    FLOG = open("%s/logs/%s_trans_param.log" % (main_dir, label), "w")
    print('## LOG FILE\n'
          '## Data: %s \n'
          '## The result of different kinds of transformation\n'
          '%s' % (datafile,
                  time.strftime('## %Y-%m-%d %H:%M:%S Begins!',
                                time.localtime(time.time()))),
          file=FLOG)
    # Log file of tex table.
    FTEX = open("%s/logs/%s_trans_param.tex" % (main_dir, label), "w")
    print('## TEX FILE\n'
          '## Data: %s \n'
          '## The result of different kinds of transformation\n'
          '%s' % (datafile,
                  time.strftime('## %Y-%m-%d %H:%M:%S Begins!',
                                time.localtime(time.time()))),
          file=FTEX)

    # Load data
    # datafile = '/home/nliu/solutions/test/a1_a2_dif.sou'
    # (RAdeg, DEdeg, dRA, dRA_err, dDC, dDC_err, COV) = np.genfromtxt(
    #     datafile, usecols=range(1, 8), unpack=True)
    # flg = np.genfromtxt(datafile, usecols=(8,), dtype=str)
    [souIVS, RAdeg, DCdeg, dRA, dRA_err, dDC, dDC_err, COV,
     ang_sep, X_a, X_d, X, flg] = DiffData

    RArad = np.deg2rad(RAdeg)  # deg -> rad
    DCrad = np.deg2rad(DCdeg)  # deg -> rad
    # ## TEST
    # COV = np.zeros_like(RArad)
    # GA - non-GA
    print('###########################################', file=FLOG)
    cat_comparison('##  GA - non-GA', RArad, DCrad, dRA, dDC,
                   dRA_err, dDC_err, COV, flg, FLOG, FTEX)

    FLOG.close()
    FTEX.close()
    print('Done!')

# --------------------------------- END --------------------------------
