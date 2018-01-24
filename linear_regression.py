##!/usr/bin/env python3
## -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:10:12 2017

Weighted linear regression fitting.

y = k * x + y0.

For data (x_i, y_i), determinate the parameter 'k' and 'y0'.

@author: Neo
"""

import numpy as np
import sys
import matplotlib.pyplot as plt  ## For test.
##------------------ FUNCTION --------------------------------
def rms_calc(x, y, k, y0):
    res = y - (k*x + y0)
    rms = np.sqrt( np.sum(res**2) / (res.size-1) )
    return res, rms
##----------------------------------------------------
def elimination(x, y, k, y0):
    n = 3.0
    res, rms = rms_calc(x, y, k, y0)
    # print(rms)
    # print(res)
    outlier = np.where( np.fabs(res) >= n * rms )
    # print(outlier)
    return outlier
##----------------------------------------------------
def parameter_calc(x, y, y_err):
    ## A * par = b
    ## calculate the elements of matrix.
    ## A
    A00 = np.sum( x**2 / y_err**2 )
    A01 = np.sum( x    / y_err**2 )
    A10 = A01
    A11 = np.sum( 1.0  / y_err**2 )
    A = np.array([[A00, A01], [A10, A11]])
    ## b
    b00 = np.sum( x*y  / y_err**2 )
    b01 = np.sum( y    / y_err**2 )
    b = np.array([b00, b01])
    ## Solve the equation. par = A^{-1} * b
    par = np.linalg.solve(A, b)
    ## formal error for par
    A_inv = np.linalg.inv(A)
    err = np.sqrt(np.diag(A_inv))
    ## correlation coeffiences
    corrcoef = np.array([ A_inv[i,j] / err[i] / err[j] \
                for j in range(par.size) for i in range(par.size)])
    corrcoef.resize((par.size, par.size))
    ## return thr result.
    return par, err, corrcoef
##----------------------------------------------------
def linear_regression(x, y, y_err):
    if len(x) < 3:
        print('There is not enough data.')
        sys.exit()
    ## Iteration.
    ## Number of source used in present and last calculation.
    num_1 = 1
    num_2 = 0
    while (num_1 != num_2):
        # print(num_1, num_2)
        num_1 = num_2
        par, err, cor = parameter_calc(x, y, y_err)
        outlier = elimination(x, y, par[0], par[1])
        num_2 = x.size - len(outlier[0])
    # print(num_1, num_2)
    return par, err, outlier, cor
##----------------------------------------------------
def linear_regression_unitwght(x, y):
    if len(x) < 3:
        print('There is not enough data.')
        sys.exit()
    ## Iteration.
    ## Number of source used in present and last calculation.
    num_1 = 1
    num_2 = 0
    y_err = np.ones_like(y)
    while (num_1 != num_2):
        # print(num_1, num_2)
        num_1 = num_2
        par, err, cor = parameter_calc(x, y, y_err)
        outlier = elimination(x, y, par[0], par[1])
        num_2 = x.size - len(outlier[0])
    # print(num_1, num_2)
    used = np.setxor1d(np.arange(x.size), outlier)
    x_used, y_used = np.extract(used, x), np.extract(used, y)
    _, rms = rms_calc(x_used, y_used, par[0], par[1]) ## post-fitting rms
    err = err * rms
    return par, err, outlier, cor
##-------------------- For test ---------------------------------
# x = np.arange(1, 17, 1)
# y = np.array([4.00, 6.40, 8.00, 8.80, 9.22, 9.50, 9.70, 9.86, 10.00, 10.20, 10.32, 10.42, 10.50, 10.55, 10.58, 10.60])
# err = np.arange(1, 17) * 0.1
# print('Unweighted fitting:')
# z1 = np.polyfit(x, y, 1)
# print('Using the library function:', z1)
# par, _, _, _ = linear_regression_unitwght(x, y)
# print('Using the script I wrote: ', par)
# print('Weighted fitting:')
# z2 = np.polyfit(x, y, deg = 1, w = 1.0/err)
# print('Using the library function:', z1)
# par, _, _, _ = linear_regression(x, y, err)
# print('Using the script I wrote: ', par)
# ## Plot
# p1 = np.poly1d(z1)
# p2 = np.poly1d(z2) # 绘制曲线
# # 原曲线
# plt.plot(x, y, 'b^', label='Origin Line')
# plt.plot(x, p1(x), 'gv-', label='Poly Fitting Line(equal-weighted)')
# plt.plot(x, p2(x), 'r*-',   label='Poly Fitting Line(normalweighted)')
# plt.axis([0, 18, 0, 18])
# plt.legend()
# # Save figure
# plt.savefig('scipy02.png', dpi=96)
# ## The result shows the script I wrote performs well.
##-------------------- END -----------------------------------
