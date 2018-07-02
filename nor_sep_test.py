#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: nor_sep_test.py
"""
Created on Fri May 25 12:03:43 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Verify my code of computation of Optical-Radio offset

"""

import numpy as np

from read_sou import read_cat
from read_GaiaDR2 import read_gaiadr2_iers_position
from list_crossmatch import list_crossmatch
from nor_sep import pos_max_calc, overall_err_calc, vlbi_gaia_sep


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
[iers_name_g, ra_g, ra_error_g,
 dec_g, dec_error_g, ra_dec_corr_g] = read_gaiadr2_iers_position(
    "/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/gaiadr2_iers.fits")

# ellipe semi-major axis
sig_pos_max_g = pos_max_calc(ra_error_g, dec_error_g, ra_dec_corr_g)

# overall formal uncertainty
overall_err_g = overall_err_calc(ra_error_g, dec_error_g, ra_dec_corr_g)


# Cross-match between between two catalogs
comsou, ind_v, ind_g = list_crossmatch(iers_name_v, iers_name_g)
# print(comsou.size)
"""2766
"""

# Verify the result of cross-match
soucom_v = np.take(iers_name_v, ind_v)
soucom_g = np.take(iers_name_g, ind_g)
for i, (comsoui, soucom_vi, soucom_gi) in enumerate(
        zip(comsou, soucom_v, soucom_g)):

    if comsoui != soucom_vi:
        print("%dth source %s are not consistented in list1 %s." %
              (i, comsoui, soucom_vi))

    if comsoui != soucom_gi:
        print("%dth source %s are not consistented in list2 %s." %
              (i, comsoui, soucom_gi))


# Extract data for common source
# VLBI data
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
[dra, ddec, dra_err, ddec_err, dra_ddec_corr,
 ang_sep, Xa, Xd, X, X2] = vlbi_gaia_sep(
    ra_com_v, dec_com_v, ra_error_com_v, dec_error_com_v, ra_dec_corr_com_v,
    ra_com_g, dec_com_g, ra_error_com_g, dec_error_com_g, ra_dec_corr_com_g)

# Log file
flog = open("../logs/vlbi-gaia-seperation-test.log", "w")

for (comsoui, ra_com_gi, dec_com_gi,
     sig_pos_max_com_gi, sig_pos_max_com_vi, num_ses_comi,
     drai, ddeci, dra_erri, ddec_erri,
     dra_ddec_corri, ang_sepi, Xai, Xdi, Xi, X2i) in zip(
        comsou, ra_com_g, dec_com_g,
        sig_pos_max_com_g, sig_pos_max_com_v, num_ses_com,
        dra, ddec, dra_err, ddec_err,
        dra_ddec_corr, ang_sep, Xa, Xd, X, X2):

    print("%8s  %14.9f  %+14.9f  %8.3f  %8.3f  %5d  %+8.3f  %+8.3f  "
          "%8.3f  %8.3f  %+5.2f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f" %
          (comsoui, ra_com_gi, dec_com_gi,
           sig_pos_max_com_gi, sig_pos_max_com_vi, num_ses_comi,
           drai, ddeci, dra_erri, ddec_erri,
           dra_ddec_corri, ang_sepi, Xai, Xdi, Xi, X2i), file=flog)

flog.close()
# --------------------------------- END --------------------------------
