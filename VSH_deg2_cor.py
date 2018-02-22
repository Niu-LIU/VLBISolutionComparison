# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 22:54:38 2017

@author: Neo

VSH function.
The full covariance matrix is used.

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
    std1 = np.sqrt(np.sum(y1r**2) / (y1r.size - 1))
    std2 = np.sqrt(np.sum(y2r**2) / (y2r.size - 1))
    indice1 = np.where(np.fabs(y1r) - n * std1 <= 0)
    indice2 = np.where(np.fabs(y2r) - n * std2 <= 0)
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
        cov[i, i + num] = C
        cov[i + num, i] = C
# Inverse it.
    wgt = np.linalg.inv(cov)
    # print(wgt[num-1, 2*num-1])
# Return the matrix.
    return wgt


# ---------------------------------------------------
def Jac_mat_deg01(RA, DE):
    # Partial array dRA and dDE, respectively.
    # For RA
    # dipole glide
    par1_d1 = -sin(RA)
    par1_d2 = cos(RA)
    par1_d3 = np.zeros_like(RA)
    # dipole rotation
    par1_r1 = -cos(RA) * sin(DE)
    par1_r2 = -sin(RA) * sin(DE)
    par1_r3 = cos(DE)
    # For DE
    # dipole glide
    par2_d1 = -cos(RA) * sin(DE)
    par2_d2 = -sin(RA) * sin(DE)
    par2_d3 = cos(DE)
    # dipole rotation
    par2_r1 = sin(RA)
    par2_r2 = -cos(RA)
    par2_r3 = np.zeros_like(RA)
# (dRA, dDE).
    pard1 = np.hstack((par1_d1, par2_d1))
    pard2 = np.hstack((par1_d2, par2_d2))
    pard3 = np.hstack((par1_d3, par2_d3))
    parr1 = np.hstack((par1_r1, par2_r1))
    parr2 = np.hstack((par1_r2, par2_r2))
    parr3 = np.hstack((par1_r3, par2_r3))
# Jacobian matrix.
    JacMatT = np.vstack((pard1, pard2, pard3, parr1, parr2, parr3))
    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT


# ---------------------------------------------------
def res_arr01(dRA, dDE, RA, DE, w):
    # Observables
    dPos = np.hstack((dRA, dDE))
# Jacobian matrix and its transpose.
    JacMat, _ = Jac_mat_deg01(RA, DE)
# Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResRA, ResDE = np.resize(ResArr, (2, dRA.size))
    return ResRA, ResDE


# ---------------------------------------------------
def VSH_deg01(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    # The 1st degree of VSH function: glide and rotation.
    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_deg01(RA, DE)
# Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cor)
# Calculate matrix A and b of matrix equation:
# A * w = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)
# Solve the equations.
##  w = (d1, d2, d3, r1, r2, r3)
    w = np.linalg.solve(A, b)
# Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())
# Correlation coefficient.
    corrcoef = np.array([cov[i, j] / sig[i] / sig[j]
                         for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
# Return the result.
    return w, sig, corrcoef


# ----------------------------------------------------
# def VSHdeg01_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE, flog):
def VSHdeg01_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    w, sig, cof = VSH_deg01(dRA, dDE, e_dRA, e_dDE, cor, RA, DE)
# Iteration.
    num1 = 1
    num2 = 0
    while(num1 != num2):
        num1 = num2
# Calculate the residual. ( O - C )
        rRA, rDE = res_arr01(dRA, dDE, RA, DE, w)
        indice = elimination(rRA, rDE)
        num2 = dRA.size - indice.size
        dRAn, dDEn, e_dRAn, e_dDEn = \
            np.take(  dRA, indice), np.take(  dDE, indice),\
            np.take(e_dRA, indice), np.take(e_dDE, indice)
        corn = np.take(cor, indice)
        RAn, DEn = np.take(RA, indice), np.take(DE, indice)
        wn, sign, cofn = VSH_deg01(dRAn, dDEn, e_dRAn, e_dDEn, corn, RAn, DEn)
        w = wn
        # print('# Number of sample: %d  %d' % (dRA.size-num1, dRA.size-num2),\
        #     file = flog)
    dRAres, dDEres = res_arr01(dRA, dDE, RA, DE, w)
    ind_outl = np.setxor1d(np.arange(dRA.size), indice)
    return wn, sign, cofn, ind_outl, dRAres, dDEres


# ---------------------------------------------------
def Jac_mat_deg02(RA, DE):
    # Partial array dRA and dDE, respectively.
    # For RA
    # dipole glide
    par1_11ER = -sin(RA)
    par1_11EI = cos(RA)  # a_{1,-1}^E
    par1_10E = np.zeros_like(RA)
    # dipole rotation
    par1_11MR = -cos(RA) * sin(DE)
    par1_11MI = -sin(RA) * sin(DE)  # a_{1,-1}^M
    par1_10M = cos(DE)
    # quadrupole
    par1_22ER = -2 * sin(2 * RA) * cos(DE)
    par1_22EI = -2 * cos(2 * RA) * cos(DE)
    par1_21ER = sin(RA) * sin(DE)
    par1_21EI = cos(RA) * sin(DE)
    par1_20E = np.zeros_like(RA)
    par1_22MR = -sin(2 * DE) * cos(2 * RA)
    par1_22MI = sin(2 * DE) * sin(2 * RA)
    par1_21MR = -cos(RA) * cos(2 * DE)
    par1_21MI = sin(RA) * cos(2 * DE)
    par1_20M = sin(2 * DE)
    # For DE
    # dipole glide
    par2_11ER = par1_11MR
    par2_11EI = par1_11MI
    par2_10E = par1_10M
    # dipole rotation
    par2_11MR = -par1_11ER
    par2_11MI = -par1_11EI
    par2_10M = -par1_10E
    # quadrupole
    par2_22ER = par1_22MI
    par2_22EI = par1_22MR
    par2_21ER = par1_21MI
    par2_21EI = par1_21MR
    par2_20E = par1_20M
    par2_22MR = -par1_22EI
    par2_22MI = -par1_22ER
    par2_21MR = -par1_21EI
    par2_21MI = -par1_21ER
    par2_20M = -par1_20E
# (dRA, dDE).
    # dipole glide
    par11ER = np.hstack((par1_11ER, par2_11ER))
    par11EI = np.hstack((par1_11EI, par2_11EI))
    par10E = np.hstack((par1_10E, par2_10E))
    # dipole rotation
    par11MR = np.hstack((par1_11MR, par2_11MR))
    par11MI = np.hstack((par1_11MI, par2_11MI))
    par10M = np.hstack((par1_10M, par2_10M))
    # quadrupole
    par22ER = np.hstack((par1_22ER, par2_22ER))
    par22EI = np.hstack((par1_22EI, par2_22EI))
    par21ER = np.hstack((par1_21ER, par2_21ER))
    par21EI = np.hstack((par1_21EI, par2_21EI))
    par20E = np.hstack((par1_20E, par2_20E))
    par22MR = np.hstack((par1_22MR, par2_22MR))
    par22MI = np.hstack((par1_22MI, par2_22MI))
    par21MR = np.hstack((par1_21MR, par2_21MR))
    par21MI = np.hstack((par1_21MI, par2_21MI))
    par20M = np.hstack((par1_20M, par2_20M))
# Jacobian matrix.
    JacMatT = np.vstack((par11ER, par11EI, par10E, \
                         # dipole glide
                         par11MR, par11MI, par10M,\
                         # dipole rotation
                         par22ER, par22EI, par21ER, par21EI, par20E, \
                         par22MR, par22MI, par21MR, par21MI, par20M))
    # quadrupole
    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT


# ---------------------------------------------------
def res_arr02(dRA, dDE, RA, DE, w):
    # Observables
    dPos = np.hstack((dRA, dDE))
# Jacobian matrix and its transpose.
    JacMat, _ = Jac_mat_deg02(RA, DE)
# Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResRA, ResDE = np.resize(ResArr, (2, dRA.size))
    return ResRA, ResDE


# ---------------------------------------------------
def VSH_deg02(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    # The 2nd degree of VSH function: glide and rotation + quadrupole.
    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_deg02(RA, DE)
# Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cor)
# Calculate matrix A and b of matrix equation:
# A * w = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)
# Solve the equations.
##  w = (d1, d2, d3, r1, r2, r3)
    w = np.linalg.solve(A, b)
# Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())
# Correlation coefficient.
    corrcoef = np.array([cov[i, j] / sig[i] / sig[j]
                         for j in range(len(w)) for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))
# Return the result.
    return w, sig, corrcoef


# ----------------------------------------------------
# def VSHdeg02_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE, flog):
def VSHdeg02_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
    w, sig, cof = VSH_deg02(dRA, dDE, e_dRA, e_dDE, cor, RA, DE)
# Iteration.
    num1 = 100
    num2 = 0

    i = 0
    # while(num2 != num1):
    while(np.fabs(num2 - num1) > 1):

        # For debugging
        i += 1
        print("Iteration times: ", i, num1, num2)

        # Calculate the residual. ( O - C )
        rRA, rDE = res_arr02(dRA, dDE, RA, DE, w)
        indice = elimination(rRA, rDE)

        num1, num2 = num2, dRA.size - indice.size

        dRAn, dDEn, e_dRAn, e_dDEn = \
            np.take(  dRA, indice), np.take(  dDE, indice),\
            np.take(e_dRA, indice), np.take(e_dDE, indice)
        corn = np.take(cor, indice)
        RAn, DEn = np.take(RA, indice), np.take(DE, indice)
        wn, sign, cofn = VSH_deg02(
            dRAn, dDEn, e_dRAn, e_dDEn, corn, RAn, DEn)
        w = wn
        # print('# Number of sample: %d  %d' % (
        # dRA.size-num1, dRA.size-num2),\
        #     file = flog)
    ind_outl = np.setxor1d(np.arange(dRA.size), indice)
    dRAres, dDEres = res_arr02(dRA, dDE, RA, DE, w)
    return wn, sign, cofn, ind_outl, dRAres, dDEres

# -------------------- MAIN ----------------------------------
# Test codes
# ## Log file.
# flog = open('../logs/test.log', 'w')
# ## Generate a data sample.
# sampleNum = 1000
# mu = 0
# sigma = 1
# np.random.seed(0)
# RA = np.random.normal(mu, sigma, sampleNum) * 2 * pi

# R1, R2, R3 = 1.3, 3.2, 5.6
# dRA = R1 * cos(RA) * sin(DE) + R2 * sin(RA) * sin(DE) - R3 * cos(DE) \
#     + np.random.normal(mu, sigma, sampleNum) * 0.5
# dDE = -R1 * sin(RA) + R2 * cos(DE) \
#     + np.random.normal(mu, sigma, sampleNum) * 0.8
# err1, err2 = np.arange(1, 1001, 1) * 0.3, np.arange(103, 1103, 1) * 0.4
# cor1 = np.random.normal(mu, sigma, sampleNum) * 0.1
# cor = np.zeros_like(err1)
# ## Using VSH degree 01.
# print('VSH deg01:')
# w, sig, corrcoef, _ = VSHdeg01_fitting(dRA, dDE, err1, err2, cor, RA, DE, flog)
# print('w = ', w)
# print('sigma: ', sig)
# ## Using VSH degree 01.
# print('VSH deg01: full covariance')
# w, sig, corrcoef, _ = VSHdeg01_fitting(dRA, dDE, err1, err2, cor1, RA, DE, flog)
# print('w = ', w)
# print('sigma: ', sig)
# ## Using VSH degree 02.
# print('VSH deg02:')
# w, sig, corrcoef, _ = VSHdeg02_fitting(dRA, dDE, err1, err2, cor, RA, DE, flog)
# print('w = ', w)
# print('sigma: ', sig)
# ## Using VSH degree 02.
# print('VSH deg02: full covariance')
# w, sig, corrcoef, _ = VSHdeg02_fitting(dRA, dDE, err1, err2, cor1, RA, DE, flog)
# print('w = ', w)
# print('sigma: ', sig)
# print('Done!')
# -------------------- END -----------------------------------
