# !/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: VSH_analysis.py
"""
Created on Thu Jan  4 12:17:39 2018

@author: Neo(liuniu@smail.nju.edu.cn)

History
N. Liu, 10 Feb 2018: change the input parameters of function
                     'vsh_analysis', replacing variable 'datafile'
                     with 'DiffData'.
N. Liu, 12 Feb 2018: add three parameters 'X_a', 'X_d', 'X' to the input
                     variables of functions 'catalog_comparison_VSH',
                     'VSH_analysis', 'apply_condition';
                     modified function 'print_outlier' to print the
                     normalized seperation information.

"""

import numpy as np
import time
from os import path
from VSH_deg2_cor import VSHdeg01_fitting, VSHdeg02_fitting
from tex_table import write_result_deg1, write_result_deg2
from vector_direction import vec6_calc
sin = np.sin
cos = np.cos


# -----------------------------  FUNCTIONS -----------------------------
def write_residual(
        RA, DE, ERRdRA, ERRdDE, RdRA, RdDE, flg, fname):

    DATA = np.transpose((RA, DE, ERRdRA, ERRdDE, RdRA, RdDE, flg))
    np.savetxt(
        fname, DATA, fmt="%5.1f " * 2 + "%9.1f " * 4 + "%d",
        delimiter=',',
        header="Residual after the VSH fitting\n"
        "RdRAig  RdDcig\n"
        "suffix: 'ig' for ICRF2 - gsf2016g\n"
        "%s" % time.strftime('##%Y-%m-%d %H:%M:%S Begins!',
                             time.localtime(time.time())))


def print_outlier(source, X_a, X_d, X, ind_outl, flog):
    '''
    '''

    outliers = np.extract(ind_outl, source)
    X_a1 = np.extract(ind_outl, X_a)
    X_d1 = np.extract(ind_outl, X_d)
    X1 = np.extract(ind_outl, X)

    print('## %d   Outliers: \n'
          '##  Source     X_a    X_d    X' % outliers.size,
          file=flog)

    for (outlier, X_ai, X_di, Xi) in zip(outliers, X_a1, X_d1, X1):
        print('## %10s  %+8.3f  %+8.3f  %7.3f' %
              (outlier, X_ai, X_di, Xi), file=flog)


def VSH_analysis(sou, d_RA, d_DE, e_dRA, e_dDE, cor, RArad, DErad,
                 flog, ftex, X_a, X_d, X):
    '''
    '''

    # Name of estimated parameters.
    x1name = ['G_1', 'G_2', 'G_3', 'R_1', 'R_2', 'R_3']
    x2name = ['E_{2,2}^{\\rm Re}', 'E_{2,2}^{\\rm Im}',
              'E_{2,1}^{\\rm Re}', 'E_{2,1}^{\\rm Im}', 'E_{2,0}',
              'M_{2,2}^{\\rm Re}', 'M_{2,2}^{\\rm Im}',
              'M_{2,1}^{\\rm Re}', 'M_{2,1}^{\\rm Im}', 'M_{2,0}']

    x1, sig1, corr1, ind_outl1, RdRA1, RdDE1 = VSHdeg01_fitting(
        d_RA, d_DE, e_dRA, e_dDE, cor, RArad, DErad)
    [gx,  gy,  gz,  wx,  wy,  wz] = x1
    [egx, egy, egz, ewx, ewy, ewz] = sig1
    (r1, alr1, der1, errr1, erralr1, errder1,
        g1, alg1, deg1, errg1, erralg1, errdeg1) = vec6_calc(x1, sig1)

# Print the result
    # write_result_deg1(x1name, x1, sig1, corr1, flog)
    # For log file.
    print('#### for degree 1:\n',
          '## Rotation component:\n',
          ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
          '=> %4d +/- %3d' % (r1, errr1), file=flog)
    print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
          (alr1, erralr1, der1, errder1), file=flog)
    print('## Glide component:\n',
          ' %+4d +/- %3d |' * 3 % (gx, egx, gy, egy, gz, egz),
          '=> %4d +/- %3d' % (g1, errg1), file=flog)
    print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
          (alg1, erralg1, deg1, errdeg1), file=flog)
    print('##   correlation coefficients are:\n', corr1, file=flog)

    # Print the outliers
    print_outlier(sou, X_a, X_d, X, ind_outl1, flog)
    # For tex file.
    print('## for degree 1:\n',
          '## Rotation component:\n',
          ' $%+4d \pm$ %3d &' * 3 % (wx, ewx, wy, ewy, wz, ewz),
          ' $%4d \pm$ %3d &(%.0f $\pm$ %.0f, $%+.0f \pm$ %.0f)' %
          (r1, errr1, alr1, erralr1, der1, errder1), file=ftex)
    print('## Glide component:\n',
          ' $%+4d \pm$ %3d &' * 3 % (gx, egx, gy, egy, gz, egz),
          ' $%4d \pm$ %3d &(%.0f $\pm$ %.0f, $%+.0f \pm$ %.0f)' %
          (g1, errg1, alg1, erralg1, deg1, errdeg1), file=ftex)
    write_result_deg1(x1name, x1, sig1, corr1, ftex)

    x2, sig2, corr2, ind_outl2, RdRA2, RdDE2 = VSHdeg02_fitting(
        d_RA, d_DE, e_dRA, e_dDE, cor, RArad, DErad)
    [gx,  gy,  gz,  wx,  wy,  wz] = x2[:6]
    [egx, egy, egz, ewx, ewy, ewz] = sig2[:6]
    (r2, alr2, der2, errr2, erralr2, errder2,
        g2, alg2, deg2, errg2, erralg2, errdeg2) = vec6_calc(
        x2[:6], sig2[:6])

# Print the result
    # For log file.
    print('#### for degree 2:\n',
          '## Rotation component:\n',
          ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
          '=> %4d +/- %3d' % (r2, errr2), file=flog)
    print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
          (alr2, erralr2, der2, errder2), file=flog)
    print('## Glide component:\n',
          ' %+4d +/- %3d |' * 3 % (gx, egx, gy, egy, gz, egz),
          '=> %4d +/- %3d' % (g2, errg2), file=flog)
    print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
          (alg2, erralg2, deg2, errdeg2), file=flog)
    print('## quadrupole component:\n',   x2[6:], file=flog)
    print('## formal uncertainties:\n', sig2[6:], file=flog)
    print('##   correlation coefficients are:\n', corr2, file=flog)
    # Print the outliers
    print_outlier(sou, X_a, X_d, X, ind_outl2, flog)

    # For tex file.
    print('## for degree 2:\n',
          '## Rotation component:\n',
          ' $%+4d \pm$ %3d &' * 3 % (wx, ewx, wy, ewy, wz, ewz),
          ' $%4d \pm$ %3d &(%.0f $\pm$ %.0f, $%+.0f \pm$ %.0f)' %
          (r2, errr2, alr2, erralr2, der2, errder2), file=ftex)
    print('## Glide component:\n',
          ' $%+4d \pm$ %3d &' * 3 % (gx, egx, gy, egy, gz, egz),
          ' $%4d \pm$ %3d &(%.0f $\pm$ %.0f, $%+.0f \pm$ %.0f)' %
          (g2, errg2, alg2, erralg2, deg2, errdeg2), file=ftex)
    write_result_deg2(x1name, x2name, x2, sig2, corr2, ftex)

    # Return the residual
    return RdRA1, RdDE1, RdRA2, RdDE2
    # return RdRA1, RdDE1, RdRA1, RdDE1


# -----------------------------------
def apply_condition(sou, d_RA, d_DE, e_dRA, e_dDE, cor,
                    RArad, DErad, flog, ftex, condition,
                    X_a, X_d, X):

    d_RA1,  d_DE1  = np.extract(condition, d_RA),  \
        np.extract(condition, d_DE)
    e_dRA1, e_dDE1 = np.extract(condition, e_dRA), \
        np.extract(condition, e_dDE)
    RArad1, DErad1 = np.extract(condition, RArad), \
        np.extract(condition, DErad)
    cor1 = np.extract(condition, cor)

    # Added on 12 Feb 2018
    X_a1 = np.extract(condition, X_a)
    X_d1 = np.extract(condition, X_d)
    X1 = np.extract(condition, X)

    RdRA1, RdDE1, RdRA2, RdDE2 = VSH_analysis(
        sou, d_RA1, d_DE1, e_dRA1, e_dDE1, cor1,
        RArad1, DErad1, flog, ftex,
        X_a1, X_d1, X1)

    return RdRA1, RdDE1, RdRA2, RdDE2


def catalog_comparison_VSH(tag, sou, d_RA, d_DE, e_dRA, e_dDE, cor,
                           RArad, DErad, flg, FLOG, FTEX,
                           X_a, X_d, X):

    print('## %s' % tag, file=FLOG)
    print('## %s' % tag, file=FTEX)

# For all sources
    print('##--------- For All sources:', file=FLOG)
    print('##--------- For All sources:', file=FTEX)
    RdRA1a, RdDE1a, RdRA2a, RdDE2a = apply_condition(
        sou, d_RA, d_DE, e_dRA, e_dDE, cor,
        RArad, DErad, FLOG, FTEX, np.ones_like(d_RA),
        X_a, X_d, X)

# For defining sources
    print('##--------- For Defining sources:', file=FLOG)
    print('##--------- For Defining sources:', file=FTEX)
    RdRA1d, RdDE1d, RdRA2d, RdDE2d = apply_condition(
        sou, d_RA, d_DE, e_dRA, e_dDE, cor,
        RArad, DErad, FLOG, FTEX, flg == 'D',
        X_a, X_d, X)

# For non-defining sources
    print('##--------- For Non-defining sources:', file=FLOG)
    print('##--------- For Non-defining sources:', file=FTEX)
    RdRA1n, RdDE1n, RdRA2n, RdDE2n = apply_condition(
        sou, d_RA, d_DE, e_dRA, e_dDE, cor,
        RArad, DErad, FLOG, FTEX, flg != 'D',
        X_a, X_d, X)

    # Return resiudal of all sources
    return RdRA1a, RdDE1a, RdRA2a, RdDE2a


def vsh_analysis(DiffData, datafile,
                 # vlbi2
                 # main_dir="/home/nliu/solutions/GalacticAberration",
                 # My MacOS
                 main_dir="/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3",
                 label=''):
    '''
    '''

    # main_dir = "/home/nliu/solutions/GalacticAberration"

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
    # # Data with covariance.
    # (RAdeg, DEdeg, D_RA, ERR_RA, D_DE, ERR_DE, COR) = np.genfromtxt(
    #     datafile, usecols=range(1, 8), unpack=True)
    # sou, flg = np.genfromtxt(
    #     datafile, usecols=(0, 8), dtype=str, unpack=True)

    [sou, RAdeg, DEdeg, D_RA, ERR_RA, D_DE, ERR_DE, COR,
     X_a, X_d, X, flg] = DiffData

    ###########################################
    RArad, DErad = np.deg2rad(RAdeg), np.deg2rad(DEdeg)  # Unit: rad

    print('###########################################', file=FLOG)
    print('###########################################', file=FTEX)

    # print("# Fitting")
    RdRA1ig, RdDE1ig, RdRA2ig, RdDE2ig = catalog_comparison_VSH(
        '##  %s ' % datafile,
        sou, D_RA, D_DE, ERR_RA, ERR_DE, COR,
        RArad, DErad, flg, FLOG, FTEX,
        X_a, X_d, X)

    # Write the post-fitting residual.
    FLG = (flg == 'D')
    # Boolean -> int, 1 for defining source while 0 for non-def
    FLG = FLG.astype(int)
    # the 1st degree
    write_residual(RAdeg, DEdeg, ERR_RA, ERR_DE, RdRA1ig, RdDE1ig, FLG,
                   '%s/logs/%s_vsh.res1' % (main_dir, label))
    # the first two degrees
    write_residual(RAdeg, DEdeg, ERR_RA, ERR_DE, RdRA2ig, RdDE2ig, FLG,
                   '%s/logs/%s_vsh.res2' % (main_dir, label))

    print('Done!')

# --------------------------------- END --------------------------------
