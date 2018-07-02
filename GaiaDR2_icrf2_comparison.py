#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: GaiaDR2_icrf2_comparison.py
"""
Created on Tue May  8 17:31:57 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from nor_sep import vlbi_gaia_sep, pos_max_calc


# -----------------------------  FUNCTIONS -----------------------------
def source_list(sou, flist):
    """
    Parameters
    ----------
    sou : array of list


    Outputs
    ----------
    """

    linefmt = "%8s " * 8 + "\\"
    i = 0

    while i + 8 <= sou.size:
        s0, s1, s2, s3, s4, s5, s6, s7 = sou[i:i+8]
        print(linefmt % (s0, s1, s2, s3, s4, s5, s6, s7), file=flist)
        i += 8

    for soui in sou[i:]:
        print("%8s " % soui, end="", file=flist)


# --------------------------------  MAIN -------------------------------
datafile = (
    "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/icrf2_gaiadr2.fits")


hdulist = fits.open(datafile, memmap=True)
tbdat = hdulist[1].data

# Gaia DR2 data
# source_id = tbdat.field("source_id")
# ref_epoch = tbdat.field("ref_epoch")
ra_gdr2 = tbdat.field("ra")
ra_error_gdr2 = tbdat.field("ra_error")
dec_gdr2 = tbdat.field("dec")
dec_error_gdr2 = tbdat.field("dec_error")
ra_dec_corr_gdr2 = tbdat.field("ra_dec_corr")

# ICRF2 data
ra_icrf2 = tbdat.field("ra_icrf2")
ra_error_icrf2 = tbdat.field("ra_err_icrf2")
dec_icrf2 = tbdat.field("dec_icrf2")
dec_error_icrf2 = tbdat.field("dec_err_icrf2")
ra_dec_corr_icrf2 = tbdat.field("ra_dec_corr_icrf2")
sou_type = tbdat.field("sou_type")


# Calculate the Gaia-VLBI positional difference
ang_sep, X_a, X_d, X = vlbi_gaia_sep(
    ra_gdr2, dec_gdr2, ra_error_gdr2, dec_error_gdr2, ra_dec_corr_gdr2,
    ra_icrf2, dec_icrf2, ra_error_icrf2, dec_error_icrf2, ra_dec_corr_icrf2)

# uas -> mas
ang_sep *= 1.e-3


# Simple statistics
# max_rho = 10
# D_outlier = ang_sep_D[ang_sep_D >= max_rho].size
# N_outlier = ang_sep_N[ang_sep_N >= max_rho].size


# print("For %d Defining sources, angular distance of %d sources "
#       "is larger than %.1f mas" % (ang_sep_D.size, D_outlier, max_rho))
# print("For %d Non-defining sources, angular distance of %d sources "
#       "is larger than %.1f mas" % (ang_sep_N.size, N_outlier, max_rho))
"""Output

For 268 Defining sources, angular distance of 18 sources is larger than 1.0 mas
For 2059 Non-defining sources, angular distance of 575 sources is larger than 1.0 mas

For 268 Defining sources, angular distance of 0 sources is larger than 5.0 mas
For 2059 Non-defining sources, angular distance of 85 sources is larger than 5.0 mas

For 268 Defining sources, angular distance of 0 sources is larger than 10.0 mas
For 2059 Non-defining sources, angular distance of 33 sources is larger than 10.0 mas
"""

# print sources with angular difference <= 5 mas
# souIERS_i = tbdat.field("iers_name_icrf2")
# souIERS_g = tbdat.field("iers_name_gaia")
# souIERS = tbdat.field("iers_name")

# max_rho = 5
# sou_selected = souIERS[ang_sep <= max_rho]

# # print(ang_sep.size)

# # print(sou_selected, sou_selected.size)

# flist = open("../data/icrf2_gaiadr2_rho5.iers_list", "w")
# source_list(sou_selected, flist)
# flist.close()

# # check
# if ang_sep_D.size + ang_sep_N.size != ang_sep.size:
#     print("%d(Defining) + %d(Non-defining) != %d(All)"
#           % (ang_sep_D.size, ang_sep_N, ang_sep.size))
#     print("Wrong number!")
#     exit()


# Plot
# angular distance
# fig, ax = plt.subplots(figsize=(8, 4))
# n_bins = 5000
# # plot the cumulative histogram
# ax.hist(ang_sep_D, n_bins, histtype='step',
#         cumulative=True, label='268Defining')
# ax.hist(ang_sep_N, n_bins, histtype='step',
#         cumulative=True, label='2059Non-Def')
# # ax.set_xlim([0, 10])
# ax.set_xscale('log')
# # ax.set_yscale('log')
# # tidy up the figure
# ax.grid(True)
# ax.legend(loc='right')
# ax.set_xlabel('Angular distance $\\rho$ (mas)')
# ax.set_ylabel('No. source')

# # plt.show()
# plt.savefig("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
#             "GaiaDR2_icrf2_rho.eps")
# plt.close()


# Orientation
from vsh_deg1_cor import VSHdeg01_fitting

# positional difference
# degree -> mirco-arc-sec
dra = (ra_gdr2 - ra_icrf2) * 3.6e9
ddec = (dec_gdr2 - dec_icrf2) * 3.6e9
dra_error = np.sqrt(ra_error_gdr2**2 + ra_error_icrf2**2) * 1.e3
ddec_error = np.sqrt(dec_error_gdr2**2 + dec_error_icrf2**2) * 1.e3
dra_ddec_cov = (ra_error_gdr2 * dec_error_gdr2 * ra_dec_corr_gdr2 +
                ra_error_icrf2 * dec_error_icrf2 * ra_dec_corr_icrf2) * 1.e6

# degree -> rad
ra_rad, dec_rad = np.deg2rad(ra_gdr2), np.deg2rad(dec_gdr2)

# Log file.
# flog = open('../logs/GaiaDR2_ICRF2_comparison.log', 'w')

# # Full 6-parameters
# w, sig, corrcoef, _, _, _ = VSHdeg01_fitting(
#     dra, ddec, dra_error, ddec_error, ra_rad, dec_rad,
#     flog, cov=dra_ddec_cov, elim_flag="None")

# print("Estimations: ")
# print("Dipole:   ", w[:3])
# print("          ", sig[:3])
# print("Rotation: ", w[3:])
# print("          ", sig[3:])
# # print('sigma: ', sig)
# print("correlation: ", corrcoef)

# flog.close()

# # Positional uncertainties
# sig_pos_max_gdr2 = pos_max_calc(
#     ra_error_gdr2, dec_error_gdr2, ra_dec_corr_gdr2)
# sig_pos_max_icrf2 = pos_max_calc(
#     ra_error_icrf2, dec_error_icrf2, ra_dec_corr_icrf2)


# # Divide sources into two subsets: Defining, Non-defining
# # defining
cond_D = sou_type == "D"
# ang_sep_D = ang_sep[cond_D]
# sig_pos_max_gdr2_D = sig_pos_max_gdr2[cond_D]
# sig_pos_max_icrf2_D = sig_pos_max_icrf2[cond_D]
dra_D = dra[cond_D]
ddec_D = ddec[cond_D]
dra_error_D = dra_error[cond_D]
ddec_error_D = ddec_error[cond_D]
dec_D = dec_gdr2[cond_D]

# # non-defining
cond_N = sou_type != "D"
# ang_sep_N = ang_sep[cond_N]
# sig_pos_max_gdr2_N = sig_pos_max_gdr2[cond_N]
# sig_pos_max_icrf2_N = sig_pos_max_icrf2[cond_N]
dra_N = dra[cond_N]
ddec_N = ddec[cond_N]
dra_error_N = dra_error[cond_N]
ddec_error_N = ddec_error[cond_N]
dec_N = dec_gdr2[cond_N]


# Positional difference vs declination
fig, (ax0, ax1) = plt.subplots(nrows=2)

ax0.errorbar(dec_D, dra_D, yerr=dra_error_D, ms=1, fmt='o', elinewidth=0.1)
ax0.set_xlabel("$Dec.\\,(^\\circ)$")
ax0.set_ylabel("$\\Delta\\alpha^*\\,(\\mu\\,as)$")
ax0.set_ylim([-1000, 1000])
# ax0.set_yscale('log')
# ax0.set_xscale('log')
# ax0.grid(True)

ax1.errorbar(dec_D, ddec_D, yerr=ddec_error_D, ms=1, fmt='o', elinewidth=0.1)
ax1.set_xlabel("$Dec.\\,(^\\circ)$")
ax1.set_ylabel("$\\Delta\\delta\\,(\\mu\\,as)$")
ax1.set_ylim([-1000, 1000])

plt.tight_layout()
plt.show()


# # angular distance vs positional uncertainty
# fig, (ax0, ax1) = plt.subplots(nrows=2)

# ax0.plot(sig_pos_max_gdr2_D, ang_sep_D, "rx", ms=1)
# ax0.plot(sig_pos_max_gdr2_N, ang_sep_N, "b*", ms=1)
# ax0.set_xlabel("$\\sigma_{pos,max}\\,GaiaDR2\\,(mas)$")
# ax0.set_ylabel("$\\rho\\,(mas)$")
# ax0.set_yscale('log')
# ax0.set_xscale('log')
# ax0.grid(True)

# ax1.plot(sig_pos_max_icrf2_D, ang_sep_D, "rx", ms=1)
# ax1.plot(sig_pos_max_icrf2_N, ang_sep_N, "b*", ms=1)
# ax1.set_xlabel("$\\sigma_{pos,max}\\,ICRF2\\,(mas)$")
# ax1.set_ylabel("$\\rho\\,(mas)$")
# ax1.set_yscale('log')
# ax1.set_xscale('log')
# ax1.grid(True)

# # plt.suptitle("Scale factor")
# plt.tight_layout()
# # plt.subplots_adjust(hspace=0.25, top=0.92)
# # plt.show()

# plt.savefig("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
#             "GaiaDR2_icrf2_rho_poserr.eps")
# plt.close()


# Normalised seperation
# num_bins = 50

# X /= 1.e6
# fig, ax = plt.subplots()
# ax.plot(X, '.')
# ax.set_yscale('log')
# plt.show()
# print(X[X <= 0])

# --------------------------------- END --------------------------------
