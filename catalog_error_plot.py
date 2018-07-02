#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: catalog_error_plot.py
"""
Created on Sun Jun 10 11:35:39 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
# My module
from read_icrf1 import read_icrf1_pos
from read_icrf2 import read_icrf2
from read_GaiaDR1 import read_gaiadr1
from read_GaiaDR2 import read_gaiadr2_iers_position
from nor_sep import pos_max_calc, overall_err_calc
from read_sou import read_cat


# -----------------------------  FUNCTIONS -----------------------------

# --------------------------------- MAIN -------------------------------
# ICRF1 catalog
[icrf_name_i1, iers_name_i1,
 RA_i1, Dec_i1, e_RA_i1, e_DE_i1, corr_i1, _] = read_icrf1_pos()


# ICRF2 catalog
[icrf_name_i2, ivs_name_i2, iers_name_i2,
 RA_i2, Dec_i2, e_RA_i2, e_DE_i2, corr_i2, _, _] = read_icrf2()

# our solutions
# OPA-A
[ivs_name_a, iers_name_a,
 RA_a, Dec_a, e_RA_a, e_DE_a, corr_a, _, _] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/opa/gaia/opa-sx-180425-noGA/"
    "opa-sx-180425-noGA.cat")

# OPA-B
[ivs_name_b, iers_name_b,
 RA_b, Dec_b, e_RA_b, e_DE_b, corr_b, _, _] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/opa/gaia/opa-sx-180425-GA00/"
    "opa-sx-180425-GA00.cat")

# OPA-C
[ivs_name_c, iers_name_c,
 RA_c, Dec_c, e_RA_c, e_DE_c, corr_c, _, _] = read_cat(
    "/Users/Neo/Astronomy/Data/VLBISolutions/opa/gaia/opa-sx-180425-GA15/"
    "opa-sx-180425-GA15.cat")

# Gaia DR1 catalog
[icrf_name_g1,
 RA_g1, Dec_g1, e_RA_g1, e_DE_g1, corr_g1, _] = read_gaiadr1()

# Gaia DR2 catalog
[iers_name_g2,
 RA_g2, e_RA_g2, Dec_g2, e_DE_g2, corr_g2] = read_gaiadr2_iers_position(
    "/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/gaiadr2_iers.fits")


# --------------------------------- END --------------------------------
