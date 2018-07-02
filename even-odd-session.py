#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: even-odd-session.py
"""
Created on Wed May 30 15:13:12 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np

from rewrite_sou import rewrite_sou
from nor_sep import pos_max_calc, overall_err_calc
from list_crossmatch import list_crossmatch
from error_plot import error_vs_numses, error_vs_numses2, \
    error_vs_numobs, error_vs_numobs2, maxerror_vs_numses, \
    maxerror_vs_numobs, maxerror_vs_numses2, maxerror_vs_numobs2, \
    overallerror_vs_numses, overallerror_vs_numobs, \
    overallerror_vs_numses2, overallerror_vs_numobs2

# -----------------------------  FUNCTIONS -----------------------------

# --------------------------------- MAIN -------------------------------
# Rewrite catalogs
print("Rewrite .sou file:")
# the solution of odd sessions
rewrite_sou("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
            "opa2018b-ga15-odd/opa2018b-ga15-odd.sou")
# the solution of even sessions
rewrite_sou("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
            "opa2018b-ga15-even/opa2018b-ga15-even.sou")

# Load data
# odd sessions
dat_dir1 = ("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
            "opa2018b-ga15-odd")

cat1 = "%s/opa2018b-ga15-odd.cat" % dat_dir1
sou1 = np.genfromtxt(cat1, dtype=str, usecols=(0,))
RA1, Dec1, RAc_err1, Dec_err1, corr1 = np.genfromtxt(
    cat1, usecols=range(2, 7), unpack=True)
num_ses1, num_obs1 = np.genfromtxt(cat1, usecols=range(10, 12),
                                   dtype=int, unpack=True)

# print("range of N_session: [%d, %d]" % (min(num_ses1), max(num_ses1)))
# print("range of N_observation: [%d, %d]" % (min(num_obs1), max(num_obs1)))

# ellipe semi-major axis
sig_pos_max1 = pos_max_calc(RAc_err1, Dec_err1, corr1)

# overall formal uncertainty
overall_err1 = overall_err_calc(RAc_err1, Dec_err1, corr1)

# Error plot for odd solution
error_vs_numses(RAc_err1, Dec_err1, num_ses1,
                "%s/plots/ga15-odd-err-Nses.eps" % dat_dir1)
#                 "%s/plots/ga15-odd-err-Nses.png" % dat_dir1)
error_vs_numobs(RAc_err1, Dec_err1, num_obs1,
                "%s/plots/ga15-odd-err-Nobs.eps" % dat_dir1)
#                 "%s/plots/ga15-odd-err-Nobs.png" % dat_dir1)

maxerror_vs_numses(sig_pos_max1, num_ses1,
                   "%s/plots/ga15-odd-maxerr-Nses.eps" % dat_dir1)
maxerror_vs_numobs(sig_pos_max1, num_obs1,
                   "%s/plots/ga15-odd-maxerr-Nobs.eps" % dat_dir1)

overallerror_vs_numses(overall_err1, num_ses1,
                       "%s/plots/ga15-odd-overallerr-Nses.eps" % dat_dir1)
overallerror_vs_numobs(overall_err1, num_obs1,
                       "%s/plots/ga15-odd-overallerr-Nobs.eps" % dat_dir1)


# even sessions
dat_dir2 = ("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/GaiaCRF/"
            "opa2018b-ga15-even")

cat2 = "%s/opa2018b-ga15-even.cat" % dat_dir2
sou2 = np.genfromtxt(cat2, dtype=str, usecols=(0,))
RA2, Dec2, RAc_err2, Dec_err2, corr2 = np.genfromtxt(
    cat2, usecols=range(2, 7), unpack=True)
num_ses2, num_obs2 = np.genfromtxt(cat2, usecols=range(10, 12),
                                   dtype=int, unpack=True)

# print("range of N_session: [%d, %d]" % (min(num_ses2), max(num_ses2)))
# print("range of N_observation: [%d, %d]" % (min(num_obs2), max(num_obs2)))

# ellipe semi-major axis
sig_pos_max2 = pos_max_calc(RAc_err2, Dec_err2, corr2)

# overall formal uncertainty
overall_err2 = overall_err_calc(RAc_err2, Dec_err2, corr2)

# Error plot for even solution
error_vs_numses(RAc_err2, Dec_err2, num_ses2,
                "%s/plots/ga15-even-err-Nses.eps" % dat_dir2)
#                 # "%s/plots/ga15-even-err-Nses.png" % dat_dir2)
error_vs_numobs(RAc_err2, Dec_err2, num_obs2,
                "%s/plots/ga15-even-err-Nobs.eps" % dat_dir2)
#                 # "%s/plots/ga15-even-err-Nobs.png" % dat_dir2)

maxerror_vs_numses(sig_pos_max2, num_ses2,
                   "%s/plots/ga15-even-maxerr-Nobs.eps" % dat_dir2)
maxerror_vs_numobs(sig_pos_max2, num_obs2,
                   "%s/plots/ga15-even-maxerr-Nses.eps" % dat_dir2)

overallerror_vs_numses(overall_err2, num_ses2,
                       "%s/plots/ga15-even-overallerr-Nses.eps" % dat_dir2)
overallerror_vs_numobs(overall_err2, num_obs2,
                       "%s/plots/ga15-even-overallerr-Nobs.eps" % dat_dir2)


# Plot for both solution
error_vs_numses2(RAc_err1, Dec_err1, num_ses1,
                 RAc_err2, Dec_err2, num_ses2,
                 "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                 "ga15-odd-even-err-Nses.eps")
# "ga15-odd-even-err-Nses.png")
error_vs_numobs2(RAc_err1, Dec_err1, num_obs1,
                 RAc_err2, Dec_err2, num_obs2,
                 "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                 "ga15-odd-even-err-Nobs.eps")
# "ga15-odd-even-err-Nobs.png")

maxerror_vs_numses2(sig_pos_max1, num_ses1,
                    sig_pos_max2, num_ses2,
                    "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                    "ga15-odd-even-maxerr-Nses.eps")
# "ga15-odd-even-err-Nses.png")
maxerror_vs_numobs2(sig_pos_max1, num_obs1,
                    sig_pos_max2, num_obs2,
                    "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                    "ga15-odd-even-maxerr-Nobs.eps")
# "ga15-odd-even-err-Nobs.png")

overallerror_vs_numses2(overall_err1, num_ses1,
                        overall_err2, num_ses2,
                        "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                        "ga15-odd-even-overallerr-Nses.eps")
# "ga15-odd-even-err-Nses.png")
overallerror_vs_numobs2(overall_err1, num_obs1,
                        overall_err2, num_obs2,
                        "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots/"
                        "ga15-odd-even-overallerr-Nobs.eps")
# "ga15-odd-even-err-Nobs.png")

# --------------------------------- END --------------------------------
