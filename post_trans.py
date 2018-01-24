# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 16:00:08 2017

Positional transformation.

@author: Neo
"""

import numpy as np
sin = np.sin
cos = np.cos
pi = np.pi
# Notice !!!
# unit for RA and DE are rad.
# ------------------ FUNCTION --------------------------------


def elimination(y1r, y2r, n=3.0):
    # std1 = np.sqrt(np.sum(y1r**2 / y1_err**2) / np.sum(y1_err**-2))
    # std2 = np.sqrt(np.sum(y2r**2 / y2_err**2) / np.sum(y2_err**-2))
    std1 = np.sqrt(np.sum(y1r**2) / (y1r.size-1))
    std2 = np.sqrt(np.sum(y2r**2) / (y2r.size-1))
    # print(std1, std2)
    indice1 = np.where(np.fabs(y1r) <= n*std1)
    indice2 = np.where(np.fabs(y2r) <= n*std2)
    indice = np.intersect1d(indice1, indice2)
    # return indice, std1, std2
    return indice
# ----------------------------------------------------


def wgt_mat(e_dRA, e_dDE, Cor):
    err = np.hstack((e_dRA, e_dDE))
# Covariance.
    cov = np.diag(err**2)
    # print(cov.shape)
# Take the correlation into consideration.
    num = e_dRA.size
    # print(num)
    for i, C in enumerate(Cor):
        cov[i, i+num] = C
        cov[i+num, i] = C
# Inverse it.
    wgt = np.linalg.inv(cov)
    # print(wgt[num-1, 2*num-1])
# Return the matrix.
    return wgt
# ---------------------------------------------------
#######################  Tran-01-01 ###################################


def Jac_mat_tran01_01(RA, DE):
    # Partial array dRA and dDE, respectively.
    par11 = -sin(DE)*cos(RA)
    par12 = sin(RA)
    par21 = -sin(DE)*sin(RA)
    par22 = -cos(RA)
    par31 = cos(DE)
    par32 = np.zeros_like(DE)
    par41 = np.rad2deg(DE)*cos(DE)  # unit: deg
    par42 = np.zeros_like(DE)
    par51 = np.zeros_like(DE)
    par52 = np.rad2deg(DE)  # unit: deg
    par61 = np.zeros_like(DE)
    par62 = np.ones_like(DE)
# (dRA, dDE).
    par1 = np.hstack((par11, par12))
    par2 = np.hstack((par21, par22))
    par3 = np.hstack((par31, par32))
    par4 = np.hstack((par41, par42))
    par5 = np.hstack((par51, par52))
    par6 = np.hstack((par61, par62))
# Jacobian matrix.
    JacMatT = np.vstack((par1, par2, par3, par4, par5, par6))
    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT
# ---------------------------------------------------


def res_arr01_01(dRA, dDE, RA, DE, w):
    # Observables
    dPos = np.hstack((dRA, dDE))
# Jacobian matrix and its transpose.
    JacMat, _ = Jac_mat_tran01_01(RA, DE)
# Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResRA, ResDE = np.resize(ResArr, (2, dRA.size))
    return ResRA, ResDE
# ---------------------------------------------------


def tran01_01(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    # Rigid rotation + slope01 + offset.
    # Equation:
    # d_RA^* = -r_x*sin(DE)*cos(RA) - r_y*sin(DE)*sin(RA) + r_z*cos(DE)\
    ##              + D_1*(DE-DE0)*cos(DE)
    ##     d_DE   = +r_x*sin(RA)         - r_y*cos(RA) \
    ##              + D_2*(DE-DE0) + B_2
    ## DE0 = 0.0
    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_tran01_01(RA, DE)
# Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cor)
# Calculate matrix A and b of matrix equation:
# A * w = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)
# Solve the equations.
##  w = (r_x, r_y, r_z, D_1, D_2, B_2)
    w = np.linalg.solve(A, b)
# Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())
# Correlation coefficient.
    corrcoef = np.array([cov[i, j]/sig[i]/sig[j]
                         for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
# Return the result.
    return w, sig, corrcoef
# ----------------------------------------------------


def tran01_01_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    w, sig, cof = tran01_01(dRA, dDE, e_dRA, e_dDE, cor, RA, DE)
# Iteration.
    num1 = 1
    num2 = 0
    while(num1 != num2):
        num1 = num2
# Calculate the residual. ( O - C )
        rRA, rDE = res_arr01_01(dRA, dDE, RA, DE, w)
        indice = elimination(rRA, rDE)
        num2 = dRA.size - indice.size
        dRAn, dDEn, e_dRAn, e_dDEn = \
            np.take(  dRA, indice), np.take(  dDE, indice),\
            np.take(e_dRA, indice), np.take(e_dDE, indice)
        corn = np.take(cor, indice)
        RAn, DEn = np.take(RA, indice), np.take(DE, indice)
        wn, sign, cofn = tran01_01(dRAn, dDEn, e_dRAn, e_dDEn, corn, RAn, DEn)
        w = wn
        # print('# Number of sample: %d  %d' % (dRA.size-num1, dRA.size-num2),\
        #     file = flog)
    ind_outl = np.setxor1d(np.arange(dRA.size), indice)
    return wn, sign, cofn, ind_outl
#######################  Tran-01-01 ###################################
# ---------------------------------------------------
#######################  Tran-01-02 ###################################


def Jac_mat_tran01_02(RA, DE):
    # Partial array dRA and dDE, respectively.
    par11 = -sin(DE)*cos(RA)
    par12 = sin(RA)
    par21 = -sin(DE)*sin(RA)
    par22 = -cos(RA)
    par31 = cos(DE)
    par32 = np.zeros_like(DE)
    par41 = np.rad2deg(DE)  # unit: deg
    par42 = np.zeros_like(DE)
    par51 = np.zeros_like(DE)
    par52 = np.rad2deg(DE)  # unit: deg
    par61 = np.zeros_like(DE)
    par62 = np.ones_like(DE)
# (dRA, dDE).
    par1 = np.hstack((par11, par12))
    par2 = np.hstack((par21, par22))
    par3 = np.hstack((par31, par32))
    par4 = np.hstack((par41, par42))
    par5 = np.hstack((par51, par52))
    par6 = np.hstack((par61, par62))
# Jacobian matrix.
    JacMatT = np.vstack((par1, par2, par3, par4, par5, par6))
    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT
# ---------------------------------------------------


def res_arr01_02(dRA, dDE, RA, DE, w):
    # Observables
    dPos = np.hstack((dRA, dDE))
# Jacobian matrix and its transpose.
    JacMat, _ = Jac_mat_tran01_02(RA, DE)
# Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResRA, ResDE = np.resize(ResArr, (2, dRA.size))
    return ResRA, ResDE
# ---------------------------------------------------


def tran01_02(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    # Rigid rotation + slope02 + offset.
    # Equation:
    # d_RA^* = -r_x*sin(DE)*cos(RA) - r_y*sin(DE)*sin(RA) + r_z*cos(DE)\
    ##              + D_1*(DE-DE0)
    ##     d_DE   = +r_x*sin(RA)         - r_y*cos(RA) \
    ##              + D_2*(DE-DE0) + B_2
    ## DE0 = 0.0
    # unit for D_1 / D_2: mas/deg
    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_tran01_02(RA, DE)
# Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cor)
# Calculate matrix A and b of matrix equation:
# A * w = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)
# Solve the equations.
##  w = (r_x, r_y, r_z, D_1, D_2, B_2)
    w = np.linalg.solve(A, b)
# Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())
# Correlation coefficient.
    corrcoef = np.array([cov[i, j]/sig[i]/sig[j]
                         for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
# Return the result.
    return w, sig, corrcoef
# ---------------------------------------------------


def tran01_02_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    w, sig, cof = tran01_01(dRA, dDE, e_dRA, e_dDE, cor, RA, DE)
# Iteration.
    num1 = 1
    num2 = 0
    while(num1 != num2):
        num1 = num2
# Calculate the residual. ( O - C )
        rRA, rDE = res_arr01_02(dRA, dDE, RA, DE, w)
        indice = elimination(rRA, rDE)
        num2 = dRA.size - indice.size
        dRAn, dDEn, e_dRAn, e_dDEn = \
            np.take(  dRA, indice), np.take(  dDE, indice),\
            np.take(e_dRA, indice), np.take(e_dDE, indice)
        corn = np.take(cor, indice)
        RAn, DEn = np.take(RA, indice), np.take(DE, indice)
        wn, sign, cofn = tran01_02(dRAn, dDEn, e_dRAn, e_dDEn, corn, RAn, DEn)
        w = wn
        # print('# Number of sample: %d  %d' % (dRA.size-num1, dRA.size-num2),\
        #     file = flog)
    ind_outl = np.setxor1d(np.arange(dRA.size), indice)
    return wn, sign, cofn, ind_outl
#######################  Tran-01-02 ###################################
# ---------------------------------------------------
#######################  Tran-01-03 ###################################


def Jac_mat_tran01_03(RA, DE):
    # Partial array dRA and dDE, respectively.
    par11 = -sin(DE)*cos(RA)
    par12 = sin(RA)
    par21 = -sin(DE)*sin(RA)
    par22 = -cos(RA)
    par31 = cos(DE)
    par32 = np.zeros_like(DE)
    par41 = np.rad2deg(DE)*cos(DE)  # unit: deg
    par42 = np.zeros_like(DE)
    par51 = np.zeros_like(DE)
    par52 = np.rad2deg(DE)  # unit: deg
# (dRA, dDE).
    par1 = np.hstack((par11, par12))
    par2 = np.hstack((par21, par22))
    par3 = np.hstack((par31, par32))
    par4 = np.hstack((par41, par42))
    par5 = np.hstack((par51, par52))
# Jacobian matrix.
    JacMatT = np.vstack((par1, par2, par3, par4, par5))
    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT
# ---------------------------------------------------


def res_arr01_03(dRA, dDE, RA, DE, w):
    # Observables
    dPos = np.hstack((dRA, dDE))
# Jacobian matrix and its transpose.
    JacMat, _ = Jac_mat_tran01_03(RA, DE)
# Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResRA, ResDE = np.resize(ResArr, (2, dRA.size))
    return ResRA, ResDE
# ---------------------------------------------------


def tran01_03(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    # Rigid rotation + slope01 + offset.
    # Equation:
    # d_RA^* = -r_x*sin(DE)*cos(RA) - r_y*sin(DE)*sin(RA) + r_z*cos(DE)\
    ##              + D_1*(DE-DE0)*cos(DE)
    ##     d_DE   = +r_x*sin(RA)         - r_y*cos(RA) \
    ##              + D_2*(DE-DE0) + B_2
    ## DE0 = 0.0
    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_tran01_03(RA, DE)
# Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cor)
# Calculate matrix A and b of matrix equation:
# A * w = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)
# Solve the equations.
##  w = (r_x, r_y, r_z, D_1, D_2, B_2)
    w = np.linalg.solve(A, b)
# Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())
# Correlation coefficient.
    corrcoef = np.array([cov[i, j]/sig[i]/sig[j]
                         for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
# Return the result.
    return w, sig, corrcoef
# ----------------------------------------------------


def tran01_03_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    w, sig, cof = tran01_03(dRA, dDE, e_dRA, e_dDE, cor, RA, DE)
# Iteration.
    num1 = 1
    num2 = 0
    while(num1 != num2):
        num1 = num2
# Calculate the residual. ( O - C )
        rRA, rDE = res_arr01_03(dRA, dDE, RA, DE, w)
        indice = elimination(rRA, rDE)
        num2 = dRA.size - indice.size
        dRAn, dDEn, e_dRAn, e_dDEn = \
            np.take(  dRA, indice), np.take(  dDE, indice),\
            np.take(e_dRA, indice), np.take(e_dDE, indice)
        corn = np.take(cor, indice)
        RAn, DEn = np.take(RA, indice), np.take(DE, indice)
        wn, sign, cofn = tran01_03(dRAn, dDEn, e_dRAn, e_dDEn, corn, RAn, DEn)
        w = wn
        # print('# Number of sample: %d  %d' % (dRA.size-num1, dRA.size-num2),\
        #     file = flog)
    ind_outl = np.setxor1d(np.arange(dRA.size), indice)
    return wn, sign, cofn, ind_outl
#######################  Tran-01-03 ###################################
# ---------------------------------------------------
#######################  Tran-02    ###################################


def Jac_mat_tran02(RA, DE):
    # Partial array dRA and dDE, respectively.
    par11 = -sin(DE)*cos(RA)
    par12 = sin(RA)
    par21 = -sin(DE)*sin(RA)
    par22 = -cos(RA)
    par31 = cos(DE)
    par32 = np.zeros_like(DE)
    par41 = np.zeros_like(DE)
    par42 = np.ones_like(DE)
# (dRA, dDE).
    par1 = np.hstack((par11, par12))
    par2 = np.hstack((par21, par22))
    par3 = np.hstack((par31, par32))
    par4 = np.hstack((par41, par42))
# Jacobian matrix.
    JacMatT = np.vstack((par1, par2, par3, par4))
    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT
# ---------------------------------------------------


def res_arr02(dRA, dDE, RA, DE, w):
    # Observables
    dPos = np.hstack((dRA, dDE))
# Jacobian matrix and its transpose.
    JacMat, _ = Jac_mat_tran02(RA, DE)
# Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResRA, ResDE = np.resize(ResArr, (2, dRA.size))
    return ResRA, ResDE
# ---------------------------------------------------


def tran02(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    # Rigid rotation + equatorial tilt.
    # Equation:
    # d_RA^* = -r_x*sin(DE)*cos(RA) - r_y*sin(DE)*sin(RA) + r_z*cos(DE)
    ##     d_DE   = +r_x*sin(RA)         - r_y*cos(RA) + dz
    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_tran02(RA, DE)
# Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cor)
# Calculate matrix A and b of matrix equation:
# A * w = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)
# Solve the equations.
##  w = (r_x, r_y, r_z, dz)
    w = np.linalg.solve(A, b)
# Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())
# Correlation coefficient.
    corrcoef = np.array([cov[i, j]/sig[i]/sig[j]
                         for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
# Return the result.
    return w, sig, corrcoef
# ---------------------------------------------------


def tran02_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    w, sig, cof = tran02(dRA, dDE, e_dRA, e_dDE, cor, RA, DE)
# Iteration.
    num1 = 1
    num2 = 0
    while(num1 != num2):
        num1 = num2
# Calculate the residual. ( O - C )
        rRA, rDE = res_arr02(dRA, dDE, RA, DE, w)
        indice = elimination(rRA, rDE)
        num2 = dRA.size - indice.size
        dRAn, dDEn, e_dRAn, e_dDEn = \
            np.take(  dRA, indice), np.take(  dDE, indice),\
            np.take(e_dRA, indice), np.take(e_dDE, indice)
        corn = np.take(cor, indice)
        RAn, DEn = np.take(RA, indice), np.take(DE, indice)
        wn, sign, cofn = tran02(dRAn, dDEn, e_dRAn, e_dDEn, corn, RAn, DEn)
        w = wn
        # print('# Number of sample: %d  %d' % (dRA.size-num1, dRA.size-num2),\
        #     file = flog)
    ind_outl = np.setxor1d(np.arange(dRA.size), indice)
    return wn, sign, cofn, ind_outl
#######################  Tran-02    ###################################
# ---------------------------------------------------
#######################  Tran-03    ###################################


def Jac_mat_tran03(RA, DE):
    # Partial array dRA and dDE, respectively.
    parx1 = -sin(DE)*cos(RA)
    parx2 = sin(RA)
    pary1 = -sin(DE)*sin(RA)
    pary2 = -cos(RA)
    parz1 = cos(DE)
    parz2 = np.zeros_like(DE)
# (dRA, dDE).
    parx = np.hstack((parx1, parx2))
    pary = np.hstack((pary1, pary2))
    parz = np.hstack((parz1, parz2))
# Jacobian matrix.
    JacMatT = np.vstack((parx, pary, parz))
    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT
# ---------------------------------------------------


def res_arr03(dRA, dDE, RA, DE, w):
    # Observables
    dPos = np.hstack((dRA, dDE))
# Jacobian matrix and its transpose.
    JacMat, _ = Jac_mat_tran03(RA, DE)
# Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResRA, ResDE = np.resize(ResArr, (2, dRA.size))
    return ResRA, ResDE
# ---------------------------------------------------


def tran03(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    # Rigid rotation.
    # Equation:
    # d_RA^* = -r_x*sin(DE)*cos(RA) - r_y*sin(DE)*sin(RA) + r_z*cos(DE)
    ##     d_DE   = +r_x*sin(RA)         - r_y*cos(RA)
    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_tran03(RA, DE)
# Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cor)
    # print(WgtMat)
# Calculate matrix A and b of matrix equation:
# A * w = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)
# Solve the equations.
##  w = (r_x, r_y, r_z)
    w = np.linalg.solve(A, b)
# Covariance.
    cov = np.linalg.inv(A)
    # print(cov)
    sig = np.sqrt(cov.diagonal())
# Correlation coefficient.
    corrcoef = np.array([cov[i, j]/sig[i]/sig[j]
                         for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
# Return the result.
    return w, sig, corrcoef
# ---------------------------------------------------


def tran03_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    w, sig, cof = tran03(dRA, dDE, e_dRA, e_dDE, cor, RA, DE)
# Iteration.
    num1 = 1
    num2 = 0
    while(num1 != num2):
        num1 = num2
# Calculate the residual. ( O - C )
        rRA, rDE = res_arr03(dRA, dDE, RA, DE, w)
        # print(rRA, rDE)
        indice = elimination(rRA, rDE)
        num2 = dRA.size - indice.size
        dRAn, dDEn, e_dRAn, e_dDEn = \
            np.take(  dRA, indice), np.take(  dDE, indice),\
            np.take(e_dRA, indice), np.take(e_dDE, indice)
        corn = np.take(cor, indice)
        RAn, DEn = np.take(RA, indice), np.take(DE, indice)
        wn, sign, cofn = tran03(dRAn, dDEn, e_dRAn, e_dDEn, corn, RAn, DEn)
        w = wn
        # print('# Number of sample: %d  %d' % (dRA.size-num1, dRA.size-num2),\
        #     file = flog)
        # print('# Number of sample: %d  %d' % (dRA.size-num1, indice.size))
    ind_outl = np.setxor1d(np.arange(dRA.size), indice)
    return wn, sign, cofn, ind_outl
#######################  Tran-02    ###################################
# -------------------- MAIN ----------------------------------
# sampleNum = 1000
# mu = 0
# sigma = 1
# np.random.seed(0)
# RA = np.random.normal(mu, sigma, sampleNum) * 2 * pi
# DE = np.random.normal(mu, sigma, sampleNum) * pi / 2
# R1, R2, R3 = 1.3, 3.2, 5.6
# dRA = +R1 * cos(RA) * sin(DE) + R2 * sin(RA) * sin(DE) - R3 * cos(DE) \
#     + np.random.normal(mu, sigma, sampleNum) * 0.5
# dDE = -R1 * sin(RA) + R2 * cos(RA) \
#     + np.random.normal(mu, sigma, sampleNum) * 0.8
# err1, err2 = np.arange(1, 1001, 1) * 0.3, np.arange(103, 1103, 1) * 0.4
# # err1, err2 = np.ones_like(RA), np.ones_like(RA)
# cor1 = np.random.normal(mu, sigma, sampleNum) * 0.9
# # cor = np.zeros_like(err1)
# ## Using new subroutine.
# print('Trans01_01:')
# # w, sig, corrcoef = tran01_01(dRA, dDE, err1, err2, cor1, RA, DE)
# w, sig, corrcoef, _ = tran01_01_fitting(dRA, dDE, err1, err2, cor1, RA, DE)
# print('w = ', w)
# print('sigma: ', sig)
# ## Using new subroutine.
# print('Trans01_02:')
# # w, sig, corrcoef = tran01_02(dRA, dDE, err1, err2, cor1, RA, DE)
# w, sig, corrcoef, _ = tran01_02_fitting(dRA, dDE, err1, err2, cor1, RA, DE)
# print('w = ', w)
# print('sigma: ', sig)
# ## Using new subroutine.
# print('Trans02:')
# # w, sig, corrcoef = tran02(dRA, dDE, err1, err2, cor1, RA, DE)
# w, sig, corrcoef, _ = tran02_fitting(dRA, dDE, err1, err2, cor1, RA, DE)
# print('w = ', w)
# print('sigma: ', sig)
# ## Using new subroutine with covariance.
# print('Trans03:')
# # w, sig, corrcoef = tran03(dRA, dDE, err1, err2, cor1, RA, DE)
# w, sig, corrcoef, _ = tran03_fitting(dRA, dDE, err1, err2, cor1, RA, DE)
# print('w = ', w)
# print('sigma: ', sig)
# print('Done!')
# -------------------- END -----------------------------------
