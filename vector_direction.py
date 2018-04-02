#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 19:57:35 2017

@author: Neo

Functions used for rotation and glide computation.

History
N. Liu 27/03/2018: give more simple expression in function 'vecpar_calc'
"""

import numpy as np
cos = np.cos
sin = np.sin


# ------------------------------  FUNCTIONS  ---------------------------
def rad_to_deg(x):
    '''x(rad) -> x(deg).
    If x is not in [0, 2*pi], then x = x + 2*pi.'''
    if not (0 <= x < 2 * np.pi):
        x = x + 2 * np.pi
    return np.rad2deg(x)


# -----------------------------------
def vecmod_calc(x):
    return np.sqrt(np.sum(x ** 2))


# -----------------------------------
def vecerr_calc(par, err):
    return np.sqrt(np.dot(par**2, err**2))


# -----------------------------------
def vecpar_calc(x, sig):
    '''x is a 3-D vector:
    x = (A*cos(alp)*cos(det), A*sin(alp)*cos(det), A*sin(det))^T.
    Calculate the amplitude A and apex (alp, det) and the uncertainty.'''
# Calculate the parameters.
    A = np.sqrt(np.sum(x ** 2))
    alp_r = np.arctan2(x[1], x[0])
#     alp_r = np.arctan(x[1] / x[0])
# ##  Judge if the x-axis is positive.
# ##  If not, add pi to it.
#     if x[0] < 0:
#         alp_r += np.pi
    det_r = np.arcsin(x[2] / A)

# ##  Calculate the uncertainties.
#     M = np.array([\
#         [( cos(alp_r)*cos(det_r) )**2 , \
#          (A*sin(alp_r)*cos(det_r))**2 , \
#          (A*cos(alp_r)*sin(det_r))**2], \
#         [( sin(alp_r)*cos(det_r) )**2 , \
#          (A*cos(alp_r)*cos(det_r))**2 , \
#          (A*sin(alp_r)*sin(det_r))**2], \
#         [(      sin(det_r)       )**2 , \
#                     0                 , \
#          (      A*cos(det_r)     )**2]])
#     M_1 = np.linalg.inv(M)
#     # print(M_1)
#     sig2 = np.dot(M_1, sig**2)
#     # print(sig2)
#     sig1 = np.sqrt(sig2)
#     errA, erralp_r, errdet_r = sig1
# Uncertainties.
    par_Amp = x / A
    # ------------------ 30 Mar 2018 ----------------------------------
    # These expressions are a little bit complicated.
    # par_alp = np.array([-cos(alp_r)**2 * x[1] / x[0]**2,
    #                     cos(alp_r)**2 / x[0], 0])
    # par_det = np.array([-x[0] * x[2],
    #                     -x[1] * x[2],
    #                     x[0]**2 + x[1]**2]) / cos(det_r) / A**3
    #
    # Let's simplifier them.
    par_alp = np.array([-sin(alp_r), cos(alp_r), 0]) / g / cos(det_r)
    par_det = np.array([-cos(alp_r) * sin(det_r),
                        -sin(alp_r) * sin(det_r),
                        cos(det_r)]) / g
    # ------------------ 30 Mar 2018 ----------------------------------

    errA, erralp_r, errdet_r = [vecerr_calc(par_Amp, sig),
                                vecerr_calc(par_alp, sig),
                                vecerr_calc(par_det, sig)]

# Rad -> deg.
    alp_d, det_d = rad_to_deg(alp_r), np.rad2deg(det_r)
    erralp_d, errdet_d = np.rad2deg(erralp_r), np.rad2deg(errdet_r)
    return A, alp_d, det_d, errA, erralp_d, errdet_d


# -----------------------------------
def vec6_calc(x, sig):
    g, err_g = x[:3], sig[:3]
    r, err_r = x[3:], sig[3:]
    R, alR, deR, errR, erralR, errdeR = vecpar_calc(r, err_r)
    G, alG, deG, errG, erralG, errdeG = vecpar_calc(g, err_g)
    return R, alR, deR, errR, erralR, errdeR,\
        G, alG, deG, errG, erralG, errdeG


# -----------------------------------
def prog_test(x, err_x, alp, err_alp, det, err_det):
    g1 = x * cos(alp) * cos(det)
    g2 = x * sin(alp) * cos(det)
    g3 = x * sin(det)
    M = np.array([
        [(cos(alp) * cos(det))**2,
         (x * sin(alp) * cos(det))**2,
         (x * cos(alp) * sin(det))**2],
        [(sin(alp) * cos(det))**2,
         (x * cos(alp) * cos(det))**2,
         (x * sin(alp) * sin(det))**2],
        [(sin(det))**2,
         0,
         (x * cos(det))**2]])
    sig2 = np.array([err_x**2, err_alp**2, err_det**2])
    err_g1, err_g2, err_g3 = np.sqrt(np.dot(M, sig2))
    return g1, g2, g3, err_g1, err_g2, err_g3

# # ------------------------------  MAIN BODY  ---------------------------
# R, errR = 85, 18
# RA, errRA = np.deg2rad(13.4), np.deg2rad(11.4)
# Dc, errDc = np.deg2rad(-11.2), np.deg2rad(12.8)
# print(prog_test(R, errR, RA, errRA, Dc, errDc))
# R, errR = 253, 17
# RA, errRA = np.deg2rad(124.3), np.deg2rad(5.6)
# Dc, errDc = np.deg2rad(42.9), np.deg2rad(4.0)
# print(prog_test(R, errR, RA, errRA, Dc, errDc))

# print('Done!')
# # ------------------------------ END -----------------------------------
