#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: itrf_trans.py
"""
Created on Sun Dec 17 06:21:29 2017

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from functools import reduce


__all__ = [
    'array_flatten', 'elim', 'elim3d',
    'wgt_mat', 'jacmat_calc', 'res_mat',
    'trans_solve', 'trans_fitting',
]


# -----------------------------  FUNCTIONS -----------------------------
def array_flatten(x1, x2, x3):
    '''
    '''
    x = np.vstack((x1, x2, x3))
    return x.flatten('F')


def elim(yr, n=3.0):
    # std1 = np.sqrt(np.sum(y1r**2 / y1_err**2) / np.sum(y1_err**-2))
    std = np.sqrt(np.sum(yr**2) / (yr.size - 1))
    indice = np.where(np.fabs(yr) - n * std <= 0)

    return indice


def elim3d(y1r, y2r, y3r, n=3.0):

    indice1 = elim(y1r, n)
    indice2 = elim(y2r, n)
    indice3 = elim(y3r, n)

    indice = reduce(np.intersect1d,
                    (indice1, indice2, indice3))

    return indice


def wgt_mat(errX, errY, errZ, corXY, corXZ, corYZ):

    err = np.hstack((errX, errY, errZ))
    # Covariance.
    cov = np.diag(err**2)
    # print(cov.shape)

    # Take the correlation into consideration.
    num = errX.size
    # print(num)
    for i in range(num):
        # X-Y
        cov[i, i+num] = corXY[i]
        cov[i+num, i] = corXY[i]
        # X-Z
        cov[i, i+num*2] = corXZ[i]
        cov[i+num*2, i] = corXZ[i]
        # Y-Z
        cov[i+num, i+num*2] = corYZ[i]
        cov[i+num*2, i+num] = corYZ[i]

    # print(cov)
    # Inverse it.
    wgt = np.linalg.inv(cov)

    return wgt

# errX = np.array([0.3, 0.2, 0.5])
# errY = np.array([1.2, 2.3, 0.3])
# errZ = np.array([0.23, 0.76, .87])
# corXY = np.array([0.2, 0.1, 0.9])
# corXZ = np.array([0, 0.5, 0.0])
# corYZ = np.array([0.1, 1.2, 0.9])
# wgt_mat(errX, errY, errZ, corXY, corXZ, corYZ)
# It's OK!


def jacmat_calc(x, y, z):
    # Jacobian matrix.

    # For \Delta_x
    # translation components
    parx_t1 = np.ones_like(x)
    parx_t2 = np.zeros_like(x)
    parx_t3 = np.zeros_like(x)
    # scale factor
    parx_d = x
    # rotation angles
    parx_r1 = np.zeros_like(x)
    parx_r2 = z
    parx_r3 = -y

    # For \Delta_y
    # translation components
    pary_t1 = np.zeros_like(y)
    pary_t2 = np.ones_like(y)
    pary_t3 = np.zeros_like(y)
    # scale factor
    pary_d = y
    # rotation angles
    pary_r1 = -z
    pary_r2 = np.zeros_like(y)
    pary_r3 = x

    # For \Delta_z
    # translation components
    parz_t1 = np.zeros_like(z)
    parz_t2 = np.zeros_like(z)
    parz_t3 = np.ones_like(z)
    # scale factor
    parz_d = z
    # rotation angles
    parz_r1 = y
    parz_r2 = -x
    parz_r3 = np.zeros_like(z)

# (dx1, dy1, dz1, dx2, dy2, dz3, ...).
    part1 = np.hstack((parx_t1, pary_t1, parz_t1))
    part2 = np.hstack((parx_t2, pary_t2, parz_t2))
    part3 = np.hstack((parx_t3, pary_t3, parz_t3))
    pard = np.hstack((parx_d, pary_d, parz_d))
    parr1 = np.hstack((parx_r1, pary_r1, parz_r1))
    parr2 = np.hstack((parx_r2, pary_r2, parz_r2))
    parr3 = np.hstack((parx_r3, pary_r3, parz_r3))

# Jacobian matrix.
    JacMatT = np.vstack((part1, part2, part3,
                         pard,
                         parr1, parr2, parr3))
    JacMat = np.transpose(JacMatT)

    return JacMat, JacMatT
# ---------------------------------------------------


def res_mat(dx, dy, dz, x, y, z, w):

    # Observables
    dPos = np.hstack((dx, dy, dz))

    # Jacobian matrix and its transpose.
    JacMat, _ = jacmat_calc(x, y, z)

    # Calculate the residual. ( O - C )
    ResArr = dPos - np.dot(JacMat, w)
    ResX, ResY, ResZ = np.resize(ResArr, (3, dx.size))

    return ResX, ResY, ResZ
# ---------------------------------------------------


def trans_solve(dx, dy, dz, x, y, z,
                errX, errY, errZ, corXY, corXZ, corYZ):

    # Jacobian matrix and its transpose.
    JacMat, JacMatT = jacmat_calc(x, y, z)
    # Weighted matrix.
    WgtMat = wgt_mat(errX, errY, errZ, corXY, corXZ, corYZ)

    # Calculate matrix A and b of matrix equation:
    # A * w = b.
    tmpMat = np.dot(JacMatT, WgtMat)
    A = np.dot(tmpMat, JacMat)

    dPos = np.hstack((dx, dy, dz))
    b = np.dot(tmpMat,  dPos)

    # Solve the equations.
    # w = (t1, t2, t3, d, r1, r2, r3)
    w = np.linalg.solve(A, b)
    # Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())
    # Correlation coefficient.
    corrcoef = np.array([cov[i, j] / sig[i] / sig[j]
                         for j in range(len(w))
                         for i in range(len(w))])
    corrcoef.resize((len(w), len(w)))

    return w, sig, corrcoef
# ----------------------------------------------------


def trans_fitting(dx, dy, dz, x, y, z,
                  errX, errY, errZ, corXY, corXZ, corYZ):

    w, sig, cof = trans_solve(dx, dy, dz, x, y, z,
                              errX, errY, errZ, corXY, corXZ, corYZ)

    # Iteration.
    num1 = 1
    num2 = 0
    num_it = 0

    while(num1 != num2):

        num_it += 1
        print("# Iterate %d times." % num_it)

        num1 = num2
        # Calculate the residual. ( O - C )
        rdx, rdy, rdz = res_mat(dx, dy, dz, x, y, z, w)
        indice = elim3d(rdx, rdy, rdz)
        num2 = dx.size - indice.size

        dxn, dyn, dzn = [np.take(dx, indice),
                         np.take(dy, indice),
                         np.take(dz, indice)]

        xn, yn, zn = [np.take(x, indice),
                      np.take(y, indice),
                      np.take(z, indice)]

        errXn, errYn, errZn = [np.take(errX, indice),
                               np.take(errY, indice),
                               np.take(errZ, indice)]

        corXYn, corXZn, corYZn = [np.take(corXY, indice),
                                  np.take(corXZ, indice),
                                  np.take(corYZ, indice)]

        wn, sign, cofn = trans_solve(dxn, dyn, dzn,
                                     xn, yn, zn,
                                     errXn, errYn, errZn,
                                     corXYn, corXZn, corYZn)
        w = wn
        # print('# Number of sample: %d  %d' % (dRA.size-num1, dRA.size-num2),\
        #     file = flog)

    rdx, rdy, rdz = res_mat(dx, dy, dz, x, y, z, w)
    ind_outl = np.setxor1d(np.arange(dx.size), indice)

    return wn, sign, cofn, ind_outl, rdx, rdy, rdz

# -------------------------- TEST ------------------------------------


def test_code():
    # For test
    N = 3000
    x1 = np.arange(N) * 1.e-1
    x2 = x1**3
    x3 = x2 + 3.4 * np.sqrt(x1)
    x = np.vstack((x1, x2, x3))
    # print(x.shape)

    xe = np.random.normal(0, 0.1, N)
    ye = np.random.normal(0, 0.2, N)
    ze = np.random.normal(0, 0.1, N)
    xerr = np.vstack((xe, ye, ze))

    R1, R2, R3 = 0.1, 0.3, 0.5
    R = np.array([[0, -R3, R2],
                  [R3, 0, -R1],
                  [-R2, R1, 0]])
    # print(R.shape)

    T = np.array([[1.5, 0.0, 0.0],
                  [0.0, 2.3, 0.0],
                  [0.0, 0.0, 0.8]])
    # print(np.dot(T, np.ones((3, N))))

    y = np.dot(T, np.ones((3, N))) + 0.95 * (x) + np.dot(R, x) + xerr

    dx1, dx2, dx3 = y
    err = np.ones_like(x1)
    # * np.random.normal(0, 1, N)
    cor = np.zeros_like(x1)

    w, sig, cof, ind, _, _, _ = trans_fitting(dx1, dx2, dx3,
                                              x1, x2, x3,
                                              err, err, err,
                                              cor, cor, cor)
    print("w:", w)
    print("sigma:", sig)
    print("coefficient:", cof)
    print("indice: ", ind)

# Time the cost time for 1 run.
from timeit import Timer
t1 = Timer("test_code()",  "from __main__ import test_code")
print(t1.timeit(1))
# --------------------------------- END --------------------------------
