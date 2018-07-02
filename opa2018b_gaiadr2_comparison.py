#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: opa2018b_gaiadr2_comparison.py
"""
Created on Sun May 27 15:28:04 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import time

from read_sou import read_cat
from read_GaiaDR2 import read_gaiadr2_iers_position
from cross_match import list_crossmatch, pos_max_calc, \
    overall_err_calc, postional_difference_calc
from VSH_analysis import vsh_analysis


# -----------------------------  FUNCTIONS -----------------------------

# -------------------------------  MAINS -------------------------------
# Load VLBI solution
# vlbi_cat = ("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
#             "opa-sx-180425-GA15.cat")
vlbi_cat = "../data/opa-sx-180425-GA15.cat"
[ivs_name_v, iers_name_v, ra_v, dec_v, ra_error_v, dec_error_v,
 ra_dec_corr_v, num_ses, num_obs] = read_cat(vlbi_cat)

# ellipe semi-major axis
sig_pos_max_v = pos_max_calc(ra_error_v, dec_error_v, ra_dec_corr_v)

# overall formal uncertainty
overall_err_v = overall_err_calc(ra_error_v, dec_error_v, ra_dec_corr_v)


# Load Gaia DR2 data
gaia_cat = "/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/gaiadr2_iers.fits"
[iers_name_g, ra_g, ra_error_g,
 dec_g, dec_error_g, ra_dec_corr_g] = read_gaiadr2_iers_position(gaia_cat)

# ellipe semi-major axis
sig_pos_max_g = pos_max_calc(ra_error_g, dec_error_g, ra_dec_corr_g)

# overall formal uncertainty
overall_err_g = overall_err_calc(ra_error_g, dec_error_g, ra_dec_corr_g)


# Cross-match between between two catalogs
comsou_iers, ind_v, ind_g = list_crossmatch(iers_name_v, iers_name_g)
# print(comsou_iers.size)
"""2766
"""

# Verify the result of cross-match
soucom_v = np.take(iers_name_v, ind_v)
soucom_g = np.take(iers_name_g, ind_g)
for i, (comsoui, soucom_vi, soucom_gi) in enumerate(
        zip(comsou_iers, soucom_v, soucom_g)):

    if comsoui != soucom_vi:
        print("%dth source %s are not consistented in list1 %s." %
              (i, comsoui, soucom_vi))

    if comsoui != soucom_gi:
        print("%dth source %s are not consistented in list2 %s." %
              (i, comsoui, soucom_gi))


# Extract data for common source
# VLBI data
comsou_ivs = np.take(ivs_name_v, ind_v)
ra_com_v = np.take(ra_v, ind_v)
dec_com_v = np.take(dec_v, ind_v)
ra_error_com_v = np.take(ra_error_v, ind_v)
dec_error_com_v = np.take(dec_error_v, ind_v)
ra_dec_corr_com_v = np.take(ra_dec_corr_v, ind_v)
num_ses_com = np.take(num_ses, ind_v)
num_obs_com = np.take(num_obs, ind_v)
sig_pos_max_com_v = np.take(sig_pos_max_v, ind_v)
overall_err_com_v = np.take(overall_err_v, ind_v)

# Gaia data
ra_com_g = np.take(ra_g, ind_g)
dec_com_g = np.take(dec_g, ind_g)
ra_error_com_g = np.take(ra_error_g, ind_g)
dec_error_com_g = np.take(dec_error_g, ind_g)
ra_dec_corr_com_g = np.take(ra_dec_corr_g, ind_g)
sig_pos_max_com_g = np.take(sig_pos_max_g, ind_g)
overall_err_com_g = np.take(overall_err_g, ind_g)


# compute the gaia-vlbi seperation
# positional difference
[dra, ddec, dra_err, ddec_err, dra_ddec_cov,
 ang_sep, Xa, Xd, X, X2] = postional_difference_calc(
    ra_com_v, dec_com_v, ra_error_com_v, dec_error_com_v, ra_dec_corr_com_v,
    ra_com_g, dec_com_g, ra_error_com_g, dec_error_com_g, ra_dec_corr_com_g)

# Positional offset
foffset = open("../logs/opa180425GA15-gaiadr2-positional-offset.log", "w")

print("# Positional offset between OPA2018b and GaiaDR2\n"
      "# Original file:\n"
      "#               %s\n"
      "#               %s\n"
      "# Columns  Units   Meaning\n"
      "#    1     --      IVS designation\n"
      "#    2     --      IERS designation\n"
      "#    3     deg     Right ascension\n"
      "#    4     deg     Declination\n"
      "#    5     mas     Semi-major axis of error ellipse for GaiaDR2\n"
      "#    6     mas     Semi-major axis of error ellipse for VLBI\n"
      "#    7     --      Number of sessions\n"
      "#    8     mas     Positional offset in R.A. cos(Decl.)\n"
      "#    9     mas     Positional offset in Decl.\n"
      "#    10    mas     Formal uncertainty of positional offset in"
      " R.A. (*cos(Dec))\n"
      "#    11    mas     Formal uncertainty of positional offset in"
      " declination\n"
      "#    12    mas^2   Correlation between positional offset between "
      "right ascension and declination\n"
      "#    13    mas     Angular separation between VLBI and Gaia positions\n"
      "#    14    --      Normalized positional offset in R.A.\n"
      "#    15    __      Normalized positional offset in Decl.\n"
      "#    16    --      Normalized positional offset\n"
      "# Created date: %s\n#"
      % (vlbi_cat, gaia_cat,
         time.strftime("%d/%m/%Y", time.localtime())), file=foffset)

for (comsou_ivsi, comsou_iersi, ra_com_gi, dec_com_gi,
     sig_pos_max_com_gi, sig_pos_max_com_vi, num_ses_comi,
     drai, ddeci, dra_erri, ddec_erri,
     dra_ddec_covi, ang_sepi, Xai, Xdi, Xi, X2i) in zip(
        comsou_ivs, comsou_iers, ra_com_g, dec_com_g,
        sig_pos_max_com_g, sig_pos_max_com_v, num_ses_com,
        dra, ddec, dra_err, ddec_err,
        dra_ddec_cov, ang_sep, Xa, Xd, X, X2):

    print("%-8s  %8s  %14.9f  %+14.9f  %8.3f  %8.3f  %5d  %+8.3f  %+8.3f  "
          "%8.3f  %8.3f  %+5.2f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f" %
          (comsou_ivsi, comsou_iersi, ra_com_gi, dec_com_gi,
           sig_pos_max_com_gi, sig_pos_max_com_vi, num_ses_comi,
           drai, ddeci, dra_erri, ddec_erri,
           dra_ddec_covi, ang_sepi, Xai, Xdi, Xi, X2i), file=foffset)

foffset.close()

flg = np.empty_like(comsou_iers)

pos_offset = [comsou_iers, ra_com_g, dec_com_g,
              dra, dra_err, ddec, ddec_err, dra_ddec_cov,
              ang_sep, Xa, Xd, X, flg]

# VSH analysis
vsh_analysis(pos_offset, vlbi_cat, label="opa-sx-180425-GA15")
# --------------------------------- END --------------------------------
