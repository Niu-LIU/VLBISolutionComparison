#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: opa-sx-20180425.py
"""
Created on Wed Jun 13 17:27:38 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import matplotlib.pyplot as plt
# My modules
from glide_calc import glide_gen
from glide_plot import glide_plot
from read_sou import read_cat


# -----------------------------  FUNCTIONS -----------------------------
def glide_decomposed_plot(gv, tag):
    """Plot for GA and non-GA component for a glide vector.

    Parameters
    ----------
    gv : array of float
        glide vector
    tag : string
        lsbel for solution

    Returns
    ----------
    g_dra : array of float
        RA offset induced by glide
    g_ddec : array of float
        Dec. offset induced by glide
    """

    GA_hat = glide_gen(1.0, 266.4, -28.9)

    # GA component
    g_GA = GA_hat * np.dot(gv, GA_hat)
    # non-GA component
    g_nonGA = gv - g_GA

    # G_nonGA = np.sqrt(np.dot(g_nonGA, g_nonGA))

    glide_plot(g_GA, "../plots/%s_GA.png" % tag,
               "GA component of %s" % tag)
    glide_plot(g_nonGA, "../plots/%s_nonGA.png" % tag,
               "Non-GA component of %s wrt. GaiaDR2" % tag)


# -------------------------------- MAIN --------------------------------
# OPA Solution
# vsh 01 parameter estimations
[x1, x1_err, x2, x2_err, x3, x3_err, x, x_err,
 g1, g1_err, g2, g2_err, g3, g3_err, g, g_err] = np.genfromtxt(
    "../logs/opa-sx-180425_Gaiadr2_vsh01.log",
    usecols=range(1, 17), unpack=True)

# Analysis of the obtained glide
# GA15
gv = np.array([g1[0], g2[0], g3[0]])
glide_decomposed_plot(gv, "OPAA")

# GA00
gv = np.array([g1[2], g2[2], g3[2]])
glide_decomposed_plot(gv, "OPAB")

# noGA
gv = np.array([g1[4], g2[4], g3[4]])
glide_decomposed_plot(g[4], "OPAC")


# Read catalogs
# opa-sx-180425-GA15
[ivs_name_nga, iers_name_nga,
 RA_nga, Dec_nga, RAc_err_nga, Dec_err_nga, corr_nga,
 num_ses_nga, num_obs_nga] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/opa/gaia/opa-sx-180425-GA00/"
    "opa-sx-180425-GA00.cat")
# --------------------------------- END --------------------------------
