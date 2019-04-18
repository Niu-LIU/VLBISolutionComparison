# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 22:54:38 2017

@author: Neo

VSH function.
The full covariance matrix is used.

# Notice !!!
# unit for RA and DE are rad.

History

N.Liu, 22/02/2018 : add some comments;
                    add new funtion 'test_code';
                    calculate the rms of residuals and reduced chi-squares.
N.Liu, 31/03/2018 : add new elimination criterior of outliers,
                      i.e., functions 'elim_angsep' and 'elim_norsep';
                    add some comments to understand these codes;
N.Liu, 04/04/2018 : add a new input parameter 'wgt' to function
                      'elim_nsigma';
                    add 2 new input parameters to functions
                      'VSHdeg01_fitting' and 'VSHdeg02_fitting';
                    add a new function 'find_good_obs';
                    functions 'VSHdeg01_fitting' and 'VSHdeg02_fitting'
                      now can use different outlier elimination algorithm;
N.Liu, 04/04/2018 : divide this code into two files "vsh_deg1_cor" and
                    "vsh_deg2_cor";

"""

import numpy as np
from numpy import sin, cos, pi, concatenate
from wrms_calc import calc_wrms, calc_2Dchi2
from nor_sep import nor_sep_calc


__all__ = ["elim_nsigma", "elim_angsep", "elim_norsep", "find_good_obs",
           "wgt_mat",
           "Jac_mat_deg02", "residual_calc02", "VSH_deg02", "VSHdeg02_fitting",
           "test_code"]


# ------------------ FUNCTION --------------------------------
def elim_nsigma(y1r, y2r, n=3.0, wgt_flag=False,
                y1_err=None, y2_err=None):
    '''An outlier elimination using n-sigma criteria.

    Parameters
    ----------
    y1r/y2r : array of float
        residuals of RA and DC
    n : float
        the strength of elimination, default value 3.0
    wgt_flag : True or False, default False
        use the rms or wrms as the unit of n

    Returns
    ----------
    ind_go : array of int
        index of good observations
    '''

    if wgt_flag:
        # wrms
        # std1 = np.sqrt(np.sum(y1r**2 / y1_err**2) / np.sum(y1_err**-2))
        # std2 = np.sqrt(np.sum(y2r**2 / y2_err**2) / np.sum(y2_err**-2))
        indice1 = np.where(np.fabs(y1r) - n * y1_err <= 0)
        indice2 = np.where(np.fabs(y2r) - n * y2_err <= 0)
        # y = np.sqrt(y1r**2 + y2r ** 2)
        # y_err = np.sqrt(y1_err**2 + y2_err**2)
        # ind_go = np.where(y - n * y_err <= 0)[0]
    else:
        # rms
        std1 = np.sqrt(np.sum(y1r**2) / (y1r.size - 1))
        std2 = np.sqrt(np.sum(y2r**2) / (y2r.size - 1))

    # indice1 = np.where(np.fabs(y1r) - n * std1 <= 0)
    # indice2 = np.where(np.fabs(y2r) - n * std2 <= 0)
    ind_go = np.intersect1d(indice1, indice2)

    # return ind_go, std1, std2
    return ind_go


# ----------------------------------------------------
def elim_angsep(angsep, pho_max=10.0e3):
    '''An outlier elimiantion based optic-radio angular seperation.

    Parameters
    ----------
    ang_sep : array of float
        angular seperation, in micro-as
    pho_max : float
        accepted maximum angular seperation, default 10.0 mas

    Returns
    ----------
    ind_go : array of int
        index of good observations
    '''

    ind_go = np.where(angsep <= pho_max)

    return ind_go


# ----------------------------------------------------
def elim_norsep(X, X_max=10.0):
    '''A outlier elimiantion based the normalized optic-radio seperation.

    Parameters
    ----------
    X : array of float
        normalized separations, unit-less.
    X_max : float
        accepted maximum X, default 10.0

    Returns
    ----------
    ind_go : array of int
        index of good observations
    '''

    ind_go = np.where(X <= X_max)

    return ind_go


def find_good_obs(dRA, dDE, e_dRA, e_dDE, cov, RA, DE, ind_go):
    '''Find the good observations based on index.

    Parameters
    ----------
    dRA/dDE : array of float
        R.A.(*cos(Dec.))/Dec. differences in uas
    e_dRA/e_dDE : array of float
        formal uncertainty of dRA(*cos(DE))/dDE in uas
    cov : array of float
        covariance between dRA and dDE in uas^2
    RA/DE : array of float
        Right ascension/Declination in radian
    ind_go : array of int
        index of good observations

    Returns
    ----------
    dRAn/dDEn : array of float
        R.A.(*cos(Dec.))/Dec. differences for good obsevations in uas
    e_dRAn/e_dDEn : array of float
        formal uncertainty of dRA(*cos(DE))/dDE good obsevations in uas
    covn : array of float
        covariance between dRA and dDE good obsevations in uas^2
    RAn/DEn : array of float
        Right ascension/Declination good obsevations in radian
    '''

    dRAn, dDEn, e_dRAn, e_dDEn = [np.take(dRA, ind_go),
                                  np.take(dDE, ind_go),
                                  np.take(e_dRA, ind_go),
                                  np.take(e_dDE, ind_go)]
    if cov is not None:
        covn = np.take(cov, ind_go)
    else:
        covn = None
    RAn, DEn = np.take(RA, ind_go), np.take(DE, ind_go)

    return dRAn, dDEn, e_dRAn, e_dDEn, covn, RAn, DEn


# ----------------------------------------------------
def wgt_mat(e_dRA, e_dDE, cov=None):
    '''Generate the weighted matrix.

    Parameters
    ----------
    e_dRA/e_dDE : array of float
        formal uncertainty of dRA(*cos(DE))/dDE in uas
    cov : array of float
        covariance between dRA and dDE in uas^2, default is None

    Returns
    ----------
    wgt : matrix
        weighted matrix used in the least squares fitting.
    '''

    err = concatenate((e_dRA, e_dDE), axis=0)

    # Covariance matrix.
    covmat = np.diag(err**2)
    # print(covmat.shape)

    if cov is not None:
        # Take the correlation into consideration.
        num = e_dRA.size
        # print(num)
        for i, C in enumerate(cov):
            covmat[i, i + num] = C
            covmat[i + num, i] = C

    # Inverse it to obtain weighted matrix.
    wgt = np.linalg.inv(covmat)
    # print(wgt[num-1, 2*num-1])

    # Return the matrix.
    return wgt


# ---------------------------------------------------
def Jac_mat_deg02(RA, DE):
    '''Generate the Jacobian matrix

    Parameters
    ----------
    RA : array of float
        right ascension in radian
    DE : array of float
        declination in radian

    Returns
    ----------
    JacMat/JacMatT : matrix
        Jacobian matrix and its transpose matrix
    '''

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
    par1_22MI = +sin(2 * DE) * sin(2 * RA)
    par1_21MR = -cos(RA) * cos(2 * DE)
    par1_21MI = +sin(RA) * cos(2 * DE)
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

    # # (dRA, dDE).
    # # dipole glide
    # par11ER = np.hstack((par1_11ER, par2_11ER))
    # par11EI = np.hstack((par1_11EI, par2_11EI))
    # par10E = np.hstack((par1_10E, par2_10E))
    # # dipole rotation
    # par11MR = np.hstack((par1_11MR, par2_11MR))
    # par11MI = np.hstack((par1_11MI, par2_11MI))
    # par10M = np.hstack((par1_10M, par2_10M))
    # # quadrupole
    # par22ER = np.hstack((par1_22ER, par2_22ER))
    # par22EI = np.hstack((par1_22EI, par2_22EI))
    # par21ER = np.hstack((par1_21ER, par2_21ER))
    # par21EI = np.hstack((par1_21EI, par2_21EI))
    # par20E = np.hstack((par1_20E, par2_20E))
    # par22MR = np.hstack((par1_22MR, par2_22MR))
    # par22MI = np.hstack((par1_22MI, par2_22MI))
    # par21MR = np.hstack((par1_21MR, par2_21MR))
    # par21MI = np.hstack((par1_21MI, par2_21MI))
    # par20M = np.hstack((par1_20M, par2_20M))

    # (dRA, dDE).
    # dipole glide
    par11ER = concatenate((par1_11ER, par2_11ER), axis=0)
    par11EI = concatenate((par1_11EI, par2_11EI), axis=0)
    par10E = concatenate((par1_10E, par2_10E), axis=0)
    # dipole rotation
    par11MR = concatenate((par1_11MR, par2_11MR), axis=0)
    par11MI = concatenate((par1_11MI, par2_11MI), axis=0)
    par10M = concatenate((par1_10M, par2_10M), axis=0)
    # quadrupole
    par22ER = concatenate((par1_22ER, par2_22ER), axis=0)
    par22EI = concatenate((par1_22EI, par2_22EI), axis=0)
    par21ER = concatenate((par1_21ER, par2_21ER), axis=0)
    par21EI = concatenate((par1_21EI, par2_21EI), axis=0)
    par20E = concatenate((par1_20E, par2_20E), axis=0)
    par22MR = concatenate((par1_22MR, par2_22MR), axis=0)
    par22MI = concatenate((par1_22MI, par2_22MI), axis=0)
    par21MR = concatenate((par1_21MR, par2_21MR), axis=0)
    par21MI = concatenate((par1_21MI, par2_21MI), axis=0)
    par20M = concatenate((par1_20M, par2_20M), axis=0)

    # Jacobian matrix.
    JacMatT = np.vstack((par11ER, par11EI, par10E,
                         # dipole glide
                         par11MR, par11MI, par10M,
                         # dipole rotation
                         par22ER, par22EI, par21ER, par21EI, par20E,
                         par22MR, par22MI, par21MR, par21MI, par20M))
    # JacMatT = concatenate((par11ER, par11EI, par10E,
    #                        # dipole glide
    #                        par11MR, par11MI, par10M,
    #                        # dipole rotation
    #                        par22ER, par22EI, par21ER, par21EI, par20E,
    #                        par22MR, par22MI, par21MR, par21MI, par20M),
    #                       axis=1)
    # quadrupole

    # # Jacobian matrix -- for test.
    # JacMatT = np.vstack((
    #     # par11ER, par11EI, par10E,
    #     # dipole glide
    #     # par11MR, par11MI, par10M,
    #     # dipole rotation
    #     par22ER, par22EI, par21ER, par21EI, par20E,
    #     par22MR, par22MI, par21MR, par21MI, par20M))

    JacMat = np.transpose(JacMatT)
    return JacMat, JacMatT


# # ---------------------------------------------------
# def res_arr02(dRA, dDE, RA, DE, W):
#     '''Calculate the residuals of RA/Dec

#     Parameters
#     ----------
#     dRA/dDE : array of float
#         R.A.(*cos(Dec.))/Dec. differences in uas
#     RA/DE : array of float
#         Right ascension/Declination in radian
#     W : matrix
#         weighted matrix

#     Returns
#     ----------
#     ResRA/ResDE : array of float
#         residual array of dRA(*cos(Dec))/dDec in uas.
#     '''

#     # Observables
#     dPos = np.hstack((dRA, dDE))

#     # Jacobian matrix and its transpose.
#     JacMat, _ = Jac_mat_deg02(RA, DE)

#     # Calculate the residual. ( O - C )
#     ResArr = dPos - np.dot(JacMat, W)
#     ResRA, ResDE = np.resize(ResArr, (2, dRA.size))

#     return ResRA, ResDE


# ---------------------------------------------------
def vsh_func01(ra, dec,
               r1, r2, r3, g1, g2, g3):
    '''VSH function of the first degree.

    Parameters
    ----------
    ra/dec : array of float
        Right ascension/Declination in radian
    r1, r2, r3 : float
        rotation parameters
    g1, g2, g3 : float
        glide parameters

    Returns
    ----------
    dra/ddec : array of float
        R.A.(*cos(Dec.))/Dec. differences in uas
    '''

    dra = [-r1 * cos(ra) * sin(dec) - r2 * sin(ra) * sin(dec)
           + r3 * cos(dec)
           - g1 * sin(ra) + g2 * cos(ra)][0]
    ddec = [+ r1 * sin(ra) - r2 * cos(ra)
            - g1 * cos(ra) * sin(dec) - g2 * sin(ra) * sin(dec)
            + g3 * cos(dec)][0]

    return dra, ddec


# ---------------------------------------------------
def vsh_func02(ra, dec,
               ER_22, EI_22, ER_21, EI_21, E_20,
               MR_22, MI_22, MR_21, MI_21, M_20):
    '''VSH function of the second degree.

    Parameters
    ----------
    ra/dec : array of float
        Right ascension/Declination in radian
    E_20, ER_21, EI_21, ER_22, EI_22 : float
        quadrupolar parameters of electric type
    M_20, MR_21, MI_21, MR_22, MI_22 : float
        quadrupolar parameters of magnetic type

    Returns
    ----------
    dra/ddec : array of float
        R.A.(*cos(Dec.))/Dec. differences in uas
    '''

    dra = [+M_20 * sin(2 * dec)
           - (MR_21 * cos(ra) - MI_21 * sin(ra)) * cos(2*dec)
           + (ER_21 * sin(ra) + EI_21 * cos(ra)) * sin(dec)
           - (MR_22 * cos(2*ra) - MI_22 * sin(2*ra)) * sin(2*dec)
           - 2*(ER_22 * sin(2*ra) + EI_22 * cos(2*ra)) * cos(dec)][0]
    ddec = [+E_20 * sin(2 * dec)
            - (MR_21 * sin(ra) + MI_21 * cos(ra)) * sin(dec)
            - (ER_21 * cos(ra) - EI_21 * sin(ra)) * cos(2*dec)
            + 2*(MR_22*sin(2*ra) + MI_22 * cos(2*ra)) * cos(dec)
            - (ER_22 * cos(2*ra) - EI_22 * sin(2*ra)) * sin(2*dec)][0]

    return dra, ddec


# ---------------------------------------------------
def vsh_func_calc02(ra, dec, param02):
    '''VSH function of the second degree.

    Parameters
    ----------
    ra/dec : array of float
        Right ascension/Declination in radian
    param02 : array of float
        estimation of rotation, glide, and quadrupolar parameters

    Returns
    ----------
    dra/ddec : array of float
        R.A.(*cos(Dec.))/Dec. differences in uas
    '''

    dra1, ddec1 = vsh_func01(ra, dec, *param02[:6])
    dra2, ddec2 = vsh_func02(ra, dec, *param02[6:])

    return dra1 + dra2, ddec1 + ddec2


# ---------------------------------------------------
def residual_calc02(dRA, dDE, RA, DE, param02):
    '''Calculate the residuals of RA/Dec

    Parameters
    ----------
    dRA/dDE : array of float
        R.A.(*cos(Dec.))/Dec. differences in uas
    RA/DE : array of float
        Right ascension/Declination in radian
    param02 : array of float
        estimation of rotation, glide, and quadrupolar parameters

    Returns
    ----------
    ResRA/ResDE : array of float
        residual array of dRA(*cos(Dec))/dDec in uas.
    '''

    # Theoritical value
    dra, ddec = vsh_func_calc02(RA, DE, param02)

    # Calculate the residual. ( O - C )
    ResRA, ResDE = dRA - dra, dDE - ddec

    return ResRA, ResDE


# ---------------------------------------------------
def VSH_deg02(dRA, dDE, e_dRA, e_dDE, cov, RA, DE):
    '''2rd degree of VSH function: glide and rotation.

    The 2nd degree of VSH function: glide and rotation + quadrupole.

    Parameters
    ----------
    dRA/dDE : array of float
        R.A.(*cos(Dec.))/Dec. differences in uas
    e_dRA/e_dDE : array of float
        formal uncertainty of dRA(*cos(DE))/dDE in uas
    cov : array of float
        covariance between dRA and dDE in uas^2
    RA/DE : array of float
        Right ascension/Declination in radian

    Returns
    ----------
    x : array of float
        estimaation of (d1, d2, d3,
                        r1, r2, r3,
                        ER_22, EI_22, ER_21, EI_21, E_20,
                        MR_22, MI_22, MR_21, MI_21, M_20) in uas
    sig : array of float
        uncertainty of x in uas
    corrmat : matrix
        matrix of correlation coefficient.
    '''

    # Jacobian matrix and its transpose.
    JacMat, JacMatT = Jac_mat_deg02(RA, DE)

    # Weighted matrix.
    WgtMat = wgt_mat(e_dRA, e_dDE, cov)

    # Calculate matrix A and b of matrix equation:
    # A * x = b.
    A = np.dot(np.dot(JacMatT, WgtMat), JacMat)
    dPos = np.hstack((dRA, dDE))
    b = np.dot(np.dot(JacMatT, WgtMat),  dPos)

    # Solve the equations.
    # x = (d1, d2, d3,
    # r1, r2, r3,
    # ER_22, EI_22, ER_21, EI_21, E_20,
    # MR_22, MI_22, MR_21, MI_21, M_20)
    x = np.linalg.solve(A, b)

    # Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())

    # Correlation coefficient.
    corrmat = np.array([cov[i, j] / sig[i] / sig[j]
                        for j in range(len(x))
                        for i in range(len(x))])
    corrmat.resize((len(x), len(x)))

    # Return the result.
    return x, sig, corrmat


# ----------------------------------------------------
# def VSHdeg02_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE):
# def VSHdeg02_fitting(dRA, dDE, e_dRA, e_dDE, cor, RA, DE, flog):
def VSHdeg02_fitting(dRA, dDE, e_dRA, e_dDE, cov, RA, DE, flog,
                     elim_flag="sigma", N=3.0):
    '''2rd-degree vsh fitting.

    Parameters
    ----------
    dRA/dDE : array of float
        R.A.(*cos(Dec.))/Dec. differences in uas
    e_dRA/e_dDE : array of float
        formal uncertainty of dRA(*cos(DE))/dDE in uas
    cov : array of float
        covariance between dRA and dDE in uas^2
    RA/DE : array of float
        Right ascension/Declination in radian
    flog :
        handlings of output file.
    elim_flag : string
        "sigma" uses n-sigma principle
        "angsep" uses angular seperation as the criteria
        "norsep" uses normalized seperation as the criteria
        "nor_ang" uses both normalized and angular seperation as the criteria
        "None" or "none" doesn't use any criteria
    N : float
        N-sigma principle for eliminating the outliers
        or
        the maximum seperation (uas)

    Returns
    ----------
    x : array of float
        estimaation of (d1, d2, d3, r1, r2, r3) in uas
    sig : array of float
        uncertainty of x in uas
    cofmat : matrix
        matrix of correlation coefficient.
    ind_outl : array of int
        index of outliers
    dRAres/dDEres: array of float
        residual array of dRA(*cos(Dec))/dDec in uas.
    '''

    # # Calculate the apriori wrms
    # meanRA, wrmsRA, stdRA = calc_wrms(dRA, e_dRA)
    # meanDE, wrmsDE, stdDE = calc_wrms(dDE, e_dDE)
    # print("# apriori statistics (weighted):\n"
    #       "#         mean for RA: %10.3f uas\n"
    #       "#         wrms for RA: %10.3f uas\n"
    #       "#          std for RA: %10.3f uas\n"
    #       "#        mean for Dec: %10.3f uas\n"
    #       "#        wrms for Dec: %10.3f uas\n"
    #       "#         std for Dec: %10.3f uas\n" %
    #       (meanRA, wrmsRA, stdRA, meanDE, wrmsDE, stdDE), file=flog)

    # # Calculate the apriori rms
    # meanRA, wrmsRA, stdRA = calc_wrms(dRA)
    # meanDE, wrmsDE, stdDE = calc_wrms(dDE)
    # print("# apriori statistics (unweighted):\n"
    #       "#         mean for RA: %10.3f uas\n"
    #       "#          rms for RA: %10.3f uas\n"
    #       "#          std for RA: %10.3f uas\n"
    #       "#        mean for Dec: %10.3f uas\n"
    #       "#         rms for Dec: %10.3f uas\n"
    #       "#         std for Dec: %10.3f uas\n" %
    #       (meanRA, wrmsRA, stdRA, meanDE, wrmsDE, stdDE), file=flog)

    # # Calculate the reduced Chi-square
    # # print("# Apriori")  # for debugging
    # print("# apriori reduced Chi-square for: %10.3f" %
    #       calc_2Dchi2(dRA, e_dRA, dDE, e_dDE, cov, reduced=True),
    #       file=flog)

    # Now we can use different criteria of elimination.
    if elim_flag is "None" or elim_flag is "none":
        x, sig, corrmat = VSH_deg02(dRA, dDE, e_dRA, e_dDE, cov, RA, DE)
        ind_go = np.arange(dRA.size)

    elif elim_flag is "sigma":
        x, sig, corrmat = VSH_deg02(dRA, dDE, e_dRA, e_dDE, cov, RA, DE)
        # Iteration.
        num1 = 1
        num2 = 0
        while(num1 != num2):
            num1 = num2

            # Calculate the residual. ( O - C )
            rRA, rDE = residual_calc02(dRA, dDE, RA, DE, x)
            ind_go = elim_nsigma(rRA, rDE, N, wgt_flag=True,
                                 y1_err=e_dRA, y2_err=e_dDE)
            # ind_go = elim_nsigma(rRA, rDE, N)
            num2 = dRA.size - ind_go.size
            # num2 = dRA.size - len(ind_go)

            dRAn, dDEn, e_dRAn, e_dDEn, covn, RAn, DEn = find_good_obs(
                dRA, dDE, e_dRA, e_dDE, cov, RA, DE, ind_go)

            xn, sign, corrmatn = VSH_deg02(dRAn, dDEn, e_dRAn, e_dDEn,
                                           covn, RAn, DEn)

            x, sig, corrmat = xn, sign, corrmatn
            print('# Number of sample: %d' % (dRA.size-num2),
                  file=flog)
    else:
        ang_sep, X_a, X_d, X = nor_sep_calc(dRA, dDE, e_dRA, e_dDE, cov)

        if elim_flag is "angsep":
            ind_go = elim_angsep(ang_sep, N)
        elif elim_flag is "norsep":
            ind_go = elim_norsep(X, N)
        elif elim_flag is "nor_ang":
            ind_go_nor = elim_norsep(X, N)
            ind_go_ang = elim_angsep(ang_sep, N)
            ind_go = np.intersect1d(ind_go_nor, ind_go_ang)
        else:
            print("ERROR: elim_flag can only be sigma, angsep, norsep,"
                  " or nor_ang!")
            exit()

        # Find all good observations
        dRAn, dDEn, e_dRAn, e_dDEn, covn, RAn, DEn = find_good_obs(
            dRA, dDE, e_dRA, e_dDE, cov, RA, DE, ind_go)
        x, sig, corrmat = VSH_deg02(dRAn, dDEn, e_dRAn, e_dDEn, covn,
                                    RAn, DEn)
        print('# Number of sample: %d' % (dRA.size-num2),
              file=flog)

    ind_outl = np.setxor1d(np.arange(dRA.size), ind_go)
    dRAres, dDEres = residual_calc02(dRA, dDE, RA, DE, x)

    # # Calculate the posteriori wrms
    # meanRA, wrmsRA, stdRA = calc_wrms(dRAres, e_dRA)
    # meanDE, wrmsDE, stdDE = calc_wrms(dDEres, e_dDE)
    # print("# posteriori statistics (weighted):\n"
    #       "#        mean for RAc: %10.3f uas\n"
    #       "#        wrms for RAc: %10.3f uas\n"
    #       "#         std for RAc: %10.3f uas\n"
    #       "#        mean for Dec: %10.3f uas\n"
    #       "#        wrms for Dec: %10.3f uas\n"
    #       "#         std for Dec: %10.3f uas\n" %
    #       (meanRA, wrmsRA, stdRA, meanDE, wrmsDE, stdDE), file=flog)

    # # Calculate the posteriori rms
    # meanRA, wrmsRA, stdRA = calc_wrms(dRAres)
    # meanDE, wrmsDE, stdDE = calc_wrms(dDEres)
    # print("# posteriori statistics (unweighted):\n"
    #       "#        mean for RAc: %10.3f uas\n"
    #       "#         rms for RAc: %10.3f uas\n"
    #       "#         std for RAc: %10.3f uas\n"
    #       "#        mean for Dec: %10.3f uas\n"
    #       "#         rms for Dec: %10.3f uas\n"
    #       "#         std for Dec: %10.3f uas\n" %
    #       (meanRA, wrmsRA, stdRA, meanDE, wrmsDE, stdDE), file=flog)

    # # Calculate the reduced Chi-square
    # # print("# posteriori")  # For debugging
    # print("# posteriori reduced Chi-square for: %10.3f" %
    #       calc_2Dchi2(dRAres, e_dRA, dDEres, e_dDE, cov, reduced=True),
    #       file=flog)

    # return x, sig, corrmat
    return x, sig, corrmat, ind_outl, dRAres, dDEres


# ----------------------------------------------------
def test_code():
    '''Code testing
    '''

    # # Log file.
    # flog = open('../logs/test.log', 'w')

    # # Generate a data sample.
    # Num = 1000
    # mu = 0
    # sigma = 1
    # # Num = 5
    # np.random.seed(2382)
    # RA = np.random.normal(pi, pi, Num)
    # DE = np.random.normal(0, pi/2, Num)
    # RA = pi
    # DE = 0

    # print("RA:  ", RA)
    # print("Dec: ", DE)

    # param = np.random.normal(0, 10, 16)
    # param = np.array([1.76405235, 0.40015721, 0.97873798,
    #                   2.2408932, 1.86755799, -0.97727788,
    #                   0.95008842,
    #                   -0.15135721, -0.10321885,
    #                   0.4105985, 0.14404357,
    #                   1.45427351,
    #                   0.76103773, 0.12167502,
    #                   0.44386323, 0.33367433])
    # Only degree 01
    # param = np.array([1.76405235, 0.40015721, 0.97873798,
    #                   2.2408932, 1.86755799, -0.97727788,
    #                   0.0, 0.0, 0.0, 0.0, 0.0,
    #                   0.0, 0.0, 0.0, 0.0, 0.0])
    # Only degree 02
    # param = np.array([0.0, 0.0, 0.0,
    #                   0.0, 0.0, 0.0,
    #                   0.95008842, -0.15135721, -0.10321885,
    #                   0.4105985, 0.14404357,
    #                   1.45427351, 0.76103773, 0.12167502,
    #                   0.44386323, 0.33367433]) * 10
    # print("Parameters: ")
    # print("Dipole: ", param[3:6])
    # print("Rotation: ", param[:3])
    # print("Quadrupolar parameters: ", param[6:])
    # param[9:11], param[7:9], param[6],
    # param[14:16], param[12:14], param[11])

    # dra1, ddec1 = vsh_func01(RA, DE, *param[:6])
    # print("dRA1:  ", dra1)
    # print("dDec1: ", ddec1)

    # dra2, ddec2 = vsh_func02(RA, DE, *param[6:])
    # print("dRA2:  ", dra2)
    # print("dDec2: ", ddec2)

    # Mock data
    # dra, ddec = vsh_func_calc02(RA, DE, param)

    # R1, R2, R3 = 1.3, 3.2, 5.6
    # dRA = R1 * cos(RA) * sin(DE) + R2 * sin(RA) * sin(DE) - R3 * cos(DE) \
    #     + np.random.normal(mu, sigma, Num) * 0.5
    # dDE = -R1 * sin(RA) + R2 * cos(DE) + np.random.normal(mu, sigma, Num) * 0.8
    # dRA = dra  # + np.random.normal(mu, sigma, Num) * 0.005
    # print(np.max(dRA))
    # dDE = ddec  # + np.random.normal(mu, sigma, Num) * 0.008
    # err1, err2 = np.arange(1, Num+1) * 0.3, np.arange(100, Num+100) * 0.4
    # err1, err2 = np.ones(Num), np.ones(Num)
    # cor1 = np.random.normal(mu, sigma, Num) * 0.1
    # cor = np.zeros_like(err1)
    # cor = None

    # Using VSH degree 02.
    # print('VSH deg02:')
    # w, sig, corrcoef, outliers, _, _ = VSHdeg02_fitting(
    #     dRA, dDE, err1, err2, cor, RA, DE, flog, elim_flag="None")
    # w, sig, corrcoef = VSHdeg02_fitting(
    # dRA, dDE, err1, err2, cor, RA, DE, flog, elim_flag="None")
    # print("Estimations: ")
    # print("Dipole: ", w[:3])
    # print("Rotation: ", w[3:6])
    # print("Quadrupolar parameters: ", w[6:])
    # print("Quadrupolar parameters: ", w)
    # print('sigma: ', sig)
    # print("correlation: ", corrcoef)
    # print(" Outlier: ", outliers.size)

    # # Using VSH degree 02.
    # print('VSH deg02: full covariance')
    # w, sig, corrcoef, _, _, _ = VSHdeg02_fitting(
    #     dRA, dDE, err1, err2, cor1, RA, DE, flog)
    # print('w = ', w)
    # print('sigma: ', sig)

    # Check the result with Titov &Lambert (2013)
    # Log file.
    flog = open('../logs/Titov_Lambert2013_check_vsh02.log', 'w')

    # Read data
    RAdeg, DEdeg, pmRA, pmDE, e_pmRA, e_pmDE = np.genfromtxt(
        "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/list429.dat",
        usecols=range(2, 8), unpack=True)

    # degree -> rad
    RArad, DErad = np.deg2rad(RAdeg), np.deg2rad(DEdeg)

    # pmRA_nor = np.fabs(pmRA / e_pmRA)
    # pmDE_nor = np.fabs(pmDE / e_pmDE)

    # pmRA_new = pmRA[pmRA_nor <= 7]
    # print(pmRA_new.size)
    # # pmDE_new = pmDE[pmDE_nor <= 7]
    # # print(pmDE_new.size)
    # # pm_nor = np.sqrt(pmRA**2 + pmDE**2) / np.sqrt(e_pmRA**2 + e_pmDE**2)
    # # pm_new = pmRA[pm_nor <= 7]
    # # print(pm_new.size)
    # RArad = RArad[pmRA_nor <= 7]
    # DErad = DErad[pmRA_nor <= 7]
    # pmRA = pmRA[pmRA_nor <= 7]
    # pmDE = pmDE[pmRA_nor <= 7]
    # e_pmRA = e_pmRA[pmRA_nor <= 7]
    # e_pmDE = e_pmDE[pmRA_nor <= 7]

    # empty
    cor = np.zeros_like(RArad)

    # Using VSH degree 02.
    print('VSH deg02:')
    w, sig, corrcoef, outliers, _, _ = VSHdeg02_fitting(
        pmRA, pmDE, e_pmRA, e_pmDE, cor, RArad, DErad,
        flog, elim_flag="sigma", N=5.5)
    # flog, elim_flag="none")
    print("Estimations: ")
    print("Dipole: ", w[:3])
    print("Rotation: ", w[3:6])
    print("Quadrupolar parameters: ", w[6:])
    print('sigma: ', sig)
    # print("correlation: ", corrcoef)
    print(" Outlier: ", outliers.size)

    # flog.close()
    # print('Done!')

    '''
    Result in paper:
    glide:  (+0.7 +/- 0.8, -6.2 +/- 0.9, -3.3 +/- 1.0)
    rotation: -(-0.5 +/- 1.2, +1.8 +/- 1.3, +0.8 +/- 1.4)
    Quadrupolar parameters:
        ER_22: +0.5 +/- 0.5
        EI_22: -1.0 +/- 0.4
        ER_21: +1.7 +/- 1.0
        EI_21: +0.9 +/- 1.1
        E_20:  +1.3 +/- 1.0
        MR_22: +2.3 +/- 0.6
        MI_22: -1.3 +/- 0.6
        MR_21: -1.8 +/- 1.0
        MI_21: +2.3 +/- 0.9
        M_20:  +0.3 +/- 0.8

    My result
    VSH deg02:
    Estimations:
    Dipole:    [ 0.55424135 -6.47696074 -3.33950471]
    Rotation:  [ 0.38691097 -0.55422749 -0.72421312]
    Quadrupolar parameters:
    [ 0.88369348 -0.77000852  2.59649833 -0.42160213  1.00972686
      0.71851149 -0.16461046 -1.29718691  1.36164007  0.46597315]
    sigma:  [0.88503584 1.02999466 0.99604722
             1.07569391 1.06509899 0.71849847
             0.45962383 0.42690982 1.19400724 1.24006569 1.06557833
             0.5587668 0.5610063  1.19137649 1.04173452 0.84731648]
    '''


test_code()
# -------------------- END -----------------------------------
