#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: GA_effect_plot.py
"""
Created on Wed Jan  3 15:07:46 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
# cos = np.cos
# sin = np.sin


# -----------------------------  FUNCTIONS -----------------------------
def GA_apriori(A, alpha_G, delta_G):
    '''Calculate the glide vector.

    Parameters
    ----------
    A : float
        amplitude of GA effect, micro-arcsec/yr
    alpha_G, delta_G : float
        direction of Galactic center, degree

    Return
    ----------
    g : array of (3,), float
        glide vector (g1, g2, g3), micro-arcsec/yr
    '''

    alpha_G = np.deg2rad(alpha_G)
    delta_G = np.deg2rad(delta_G)

    g1 = A * cos(alpha_G) * cos(delta_G)
    g2 = A * sin(alpha_G) * cos(delta_G)
    g3 = A * sin(delta_G)

    return np.array([g1, g2, g3])


def GA_effect_calc(RA, DC, g, flag="arc"):
    '''Calculate the positional offset caused by GA effect


    Parameters
    ----------
    RA, DC : array, float
        Right ascension / Declination of radio sources, degree
    g : array of (3,0), float
        glide vector
    flag : string
        determine whether return the dRA or dRA * cos(DC)

    Returns
    ----------
    dRA, dDC : array, float
        positional offsets caused by GA effect
    '''

    RAr = np.deg2rad(RA)
    DCr = np.deg2rad(DC)

    g1, g2, g3 = g

    if flag is "arc":
        dRA = -g1 * sin(RAr) + g2 * cos(RAr)
    else:
        dRA = (-g1 * sin(RAr) + g2 * cos(RAr)) / cos(DCr)

    dDC = -g1 * cos(RAr) * sin(DCr) - g2 * sin(RAr) * sin(DCr) + g3 * cos(DCr)

    return dRA, dDC


def GA_SF_calc(RA, DC, g):
    '''Calculate the scale factor F caused by GA effect.


    F = g1*cos(RA)*cos(DC) + g2*sin(RA)*cos(DC) + g3*sin(DC)


    Parameters
    ----------
    RA, DC : array, float
        Right ascension / Declination of radio sources, degree
    g : array of (3,0), float
        glide vector

    Returns
    ----------
    F : array, float
        scale factor caused by GA effect
    '''

    RAr = np.deg2rad(RA)
    DCr = np.deg2rad(DC)

    g1, g2, g3 = g

    F = g1 * cos(RAr) * cos(DCr) + g2 * sin(RAr) * cos(DCr) + g3 * sin(DCr)

    return F


def GAEffectPlot(RA, DC, dRA, dDC, plot_dir):
    '''
    '''

    fig, (ax0, ax1) = plt.subplots(nrows=2)

    ax0.plot(RA, dRA, ".", ms=1)
    ax0.set_xlabel("R.A. ($\circ$)")
    ax0.set_xlim([0, 360])
    ax0.set_ylim([-5, 5])
    ax0.set_xticks(np.arange(0, 361, 30))

    ax1.plot(DC, dRA, ".", ms=1)
    ax1.set_xlabel("Dec. ($\circ$)")
    ax1.set_xlim([-90, 90])
    ax1.set_ylim([-5, 5])
    ax1.set_xticks(np.arange(-90, 91, 30))

    plt.suptitle("$\mu _{\\alpha^*}\,(\mu as/yr)$")
    # plt.tight_layout()
    plt.subplots_adjust(hspace=0.28, top=0.92)

    # plt.tight_layout()
    plt.savefig("%s/GA_dRA.eps" % plot_dir)
    plt.close()

    fig, (ax0, ax1) = plt.subplots(nrows=2)

    ax0.plot(RA, dDC, ".", ms=1)
    ax0.set_xlabel("R.A. ($\circ$)")
    ax0.set_xlim([0, 360])
    ax0.set_xticks(np.arange(0, 360, 30))
    ax0.set_xticks(np.arange(0, 361, 30))

    ax1.plot(DC, dDC, ".", ms=1)
    ax1.set_xlabel("Dec. ($\circ$)")
    ax1.set_xlim([-90, 90])
    ax1.set_xticks(np.arange(-90, 90, 30))
    ax1.set_xticks(np.arange(-90, 91, 30))

    plt.suptitle("$\mu _\delta\,(\mu as/yr)$")
    # plt.tight_layout()
    plt.subplots_adjust(hspace=0.28, top=0.92)

    # plt.tight_layout()
    plt.savefig("%s/GA_dDC.eps" % plot_dir)
    plt.close()


def GASFPlot(RA, DC, F, plot_dir):
    '''
    '''

    fig, (ax0, ax1) = plt.subplots(nrows=2)

    ax0.plot(RA, F, ".", ms=1)
    ax0.set_xlabel("R.A. ($\circ$)")
    ax0.set_xlim([0, 360])
    ax0.set_ylim([-5, 5])
    ax0.set_xticks(np.arange(0, 361, 30))

    ax1.plot(DC, F, ".", ms=1)
    ax1.set_xlabel("Dec. ($\circ$)")
    ax1.set_xlim([-90, 90])
    ax1.set_ylim([-5, 5])
    ax1.set_xticks(np.arange(-90, 91, 30))

    plt.suptitle("Scale factor")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.28, top=0.92)

    plt.savefig("%s/GA_SF.eps" % plot_dir)
    plt.close()


# A priori value for Galactic Aberration Effect
A = 5.0
alpha_G = 267.0
delta_G = -29.0
g = GA_apriori(A, alpha_G, delta_G)
# print(g)
# Initial RA/Dec
N = 10000
RA = np.random.random(N) * 360
DC = (np.random.random(N) - 0.5) * 180

# path to store the plots, without ending with a slash('/')
# plot_dir = "/home/nliu/solutions/GalacticAberration"
plot_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots"

#
dRA, dDC = GA_effect_calc(RA, DC, g, "non")
GAEffectPlot(RA, DC, dRA, dDC, plot_dir)

F = GA_SF_calc(RA, DC, g)
GASFPlot(RA, DC, F, plot_dir)
# --------------------------------- END --------------------------------
