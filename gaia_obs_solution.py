#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: gaia_obs_solution.py
"""
Created on Wed Jun 13 11:49:42 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import matplotlib.pyplot as plt
# My modules
from cross_match import pos_Xmatch, postional_difference_calc
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
        label for solution

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

    G_nonGA = np.sqrt(np.dot(g_nonGA, g_nonGA))

    glide_plot(g_GA, "../plots/%s_GA.png" % tag,
               "GA component of %s" % tag)
    glide_plot(g_nonGA, "../plots/%s_nonGA.png" % tag,
               "Non-GA component of %s wrt. GaiaDR2" % tag)


# -------------------------------- MAIN --------------------------------
# Gaia time span solution
# GaiaDR1/2 GA00/GA15/noGA
# vsh 01 parameter estimations
[x1, x1_err, x2, x2_err, x3, x3_err, x, x_err,
 g1, g1_err, g2, g2_err, g3, g3_err, g, g_err] = np.genfromtxt(
    "../logs/gaia_obs_Gaiadr2_vsh01.log",
    usecols=range(1, 17), unpack=True)

# GaiaDR2-GA00
gv = np.array([g1[3], g2[3], g3[3]])
glide_decomposed_plot(gv, "GaiaDR2-GA00")

# GaiaDR2-GA15
gv = np.array([g1[4], g2[4], g3[4]])
glide_decomposed_plot(gv, "GaiaDR2-GA15")

# GaiaDR2-noGA
gv = np.array([g1[5], g2[5], g3[5]])
glide_decomposed_plot(gv, "GaiaDR2-noGA")

# GaiaDR2-GA00b
gv = np.array([g1[6], g2[6], g3[6]])
glide_decomposed_plot(gv, "GaiaDR2-GA00b")

# GaiaDR2-GA15b
gv = np.array([g1[7], g2[7], g3[7]])
glide_decomposed_plot(gv, "GaiaDR2-GA15b")

# GaiaDR2-noGAb
gv = np.array([g1[8], g2[8], g3[8]])
glide_decomposed_plot(gv, "GaiaDR2-noGAb")


# Read Catalog
# gaiadr2-timespan-noGA
[ivs_name_nga, iers_name_nga,
 RA_nga, Dec_nga, RAc_err_nga, Dec_err_nga, corr_nga,
 num_ses_nga, num_obs_nga] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GalacticAberration/"
    "gaiadr2-timespan-noGA/gaiadr2-timespan-noGA.cat")

# gaiadr2-timespan-noGAb
[ivs_name_ngab, iers_name_ngab,
 RA_ngab, Dec_ngab, RAc_err_ngab, Dec_err_ngab, corr_ngab,
 num_ses_ngab, num_obs_ngab] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GalacticAberration/"
    "gaiadr2-timespan-noGAb/gaiadr2-timespan-noGAb.cat")

# gaiadr2-timespan-GA00
[ivs_name_g00, iers_name_g00,
 RA_g00, Dec_g00, RAc_err_g00, Dec_err_g00, corr_g00,
 num_ses_g00, num_obs_g00] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GalacticAberration/"
    "gaiadr2-timespan-GA00/gaiadr2-timespan-GA00.cat")

# gaiadr2-timespan-GA15
[ivs_name_g15, iers_name_g15,
 RA_g15, Dec_g15, RAc_err_g15, Dec_err_g15, corr_g15,
 num_ses_g15, num_obs_g15] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GalacticAberration/"
    "gaiadr2-timespan-GA15/gaiadr2-timespan-GA15.cat")

# Use the solution gaiadr2-timespan-noGA as reference
# gaiadr2-timespan-noGAb - gaiadr2-timespan-noGA
[soucom,
 RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
 RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n] = pos_Xmatch(
    iers_name_nga, RA_nga, RAc_err_nga, Dec_nga, Dec_err_nga, corr_nga,
    iers_name_ngab, RA_ngab, RAc_err_ngab, Dec_ngab, Dec_err_ngab, corr_ngab)
[dRA_1, dDC_1, dRA_err_1, dDC_err_1, cov_1,
 ang_sep_1, X_a_1, X_d_1, X_1, X2_1] = postional_difference_calc(
    RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
    RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n)

# gaiadr2-timespan-noGAb - gaiadr2-timespan-GA00
[soucom,
 RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
 RA3n, RAc_err3n, Dec3n, Dec_err3n, cor3n] = pos_Xmatch(
    iers_name_nga, RA_nga, RAc_err_nga, Dec_nga, Dec_err_nga, corr_nga,
    iers_name_g00, RA_g00, RAc_err_g00, Dec_g00, Dec_err_g00, corr_g00)
[dRA_2, dDC_2, dRA_err_2, dDC_err_2, cov_2,
 ang_sep_2, X_a_2, X_d_2, X_2, X2_2] = postional_difference_calc(
    RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
    RA3n, RAc_err3n, Dec3n, Dec_err3n, cor3n)

# gaiadr2-timespan-noGAb - gaiadr2-timespan-GA15
[soucom,
 RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
 RA4n, RAc_err4n, Dec4n, Dec_err4n, cor4n] = pos_Xmatch(
    iers_name_nga, RA_nga, RAc_err_nga, Dec_nga, Dec_err_nga, corr_nga,
    iers_name_g15, RA_g15, RAc_err_g15, Dec_g15, Dec_err_g15, corr_g15)
[dRA_3, dDC_3, dRA_err_3, dDC_err_3, cov_3,
 ang_sep_3, X_a_3, X_d_3, X_3, X2_3] = postional_difference_calc(
    RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
    RA4n, RAc_err4n, Dec4n, Dec_err4n, cor4n)

# Plot for positional differences
# RA difference
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
ax1.errorbar(Dec2n, dRA_1, yerr=dRA_err_1, fmt='.',
             elinewidth=0.05, markersize=1)
ax2.errorbar(Dec3n, dRA_2, yerr=dRA_err_2, fmt='.',
             elinewidth=0.05, markersize=1)
ax3.errorbar(Dec4n, dRA_3, yerr=dRA_err_3, fmt='.',
             elinewidth=0.05, markersize=1)
# Ticks
ax1.set_xticks(np.arange(-90, 91, 30))
ax1.set_ylabel("noGAb (mas)")
ax1.set_ylim([-0.15, 0.15])

# ax2.set_ylim([-0.15, 0.15])
ax2.set_ylabel("GA00 (mas)")

# ax3.set_ylim([-0.25, 0.25])
ax3.set_ylabel("GA15 (mas)")
ax3.set_xlabel("Declination (degree)")

ax1.set_title("$\Delta\\alpha^*$")
fig.subplots_adjust(hspace=0)
# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.savefig("../plots/gaia_obs_RA_diff.png")
plt.close()

# Declination differences
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
ax1.errorbar(Dec2n, dDC_1, yerr=dDC_err_1, fmt='.',
             elinewidth=0.05, markersize=1)
ax2.errorbar(Dec3n, dDC_2, yerr=dDC_err_2, fmt='.',
             elinewidth=0.05, markersize=1)
ax3.errorbar(Dec4n, dDC_3, yerr=dDC_err_3, fmt='.',
             elinewidth=0.05, markersize=1)
# Ticks
ax1.set_xticks(np.arange(-90, 91, 30))
ax1.set_ylabel("noGAb (mas)")
ax1.set_ylim([-0.15, 0.15])

# ax2.set_ylim([-0.15, 0.15])
ax2.set_ylabel("GA00 (mas)")

# ax3.set_ylim([-0.25, 0.25])
ax3.set_ylabel("GA15 (mas)")
ax3.set_xlabel("Declination (degree)")

ax1.set_title("$\Delta\delta$")
fig.subplots_adjust(hspace=0)
# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.savefig("../plots/gaia_obs_Dec_diff.png")
plt.close()
# --------------------------------- END --------------------------------
