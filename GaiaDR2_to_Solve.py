#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: Gaiadr22Solve.py
"""
Created on Thu Apr 26 17:43:09 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from astropy.io import fits
from write_solvesrc import write_solvesrc, write_NNRS


# -----------------------------  FUNCTIONS -----------------------------
def write_icrf2_solve(datafile, Gaiasouivs, fout):
    souivs = np.genfromtxt(
        datafile, usecols=(1,), dtype=str, unpack=True)

    RAh, RAm, DCd, DCm = np.genfromtxt(datafile,
                                       usecols=(4, 5, 7, 8),
                                       dtype=str,
                                       unpack=True)
    RAs, DCs, DC_erras = np.genfromtxt(datafile,
                                       usecols=(6, 9, 11),
                                       unpack=True)
    # as -> mas
    DC_err = DC_erras * 1.e3
    x = 0

    linefmt = "    %8s  %2s %2s %11.8f   %3s %2s  %10.7f   %7.3f   %s"

    for (souni, RAhi, RAmi, RAsi, DCdi, DCmi, DCsi, DC_erri) in zip(
            souivs, RAh, RAm, RAs, DCd, DCm, DCs, DC_err):

        if souni not in Gaiasouivs:
            print(linefmt % (souni, RAhi, RAmi, RAsi,
                             DCdi, DCmi, DCsi, DC_erri, 'icrf2'),
                  file=fout)
            x += 1

    print("number from ICRF2 :", x)


# -------------------------------  MAIN -------------------------------
'''Write the Gaia DR1 source position into the Solve format.
'''

print('======================= BEGIN ============================\n'
      '# Write Gaia DR2 qso into Solve format')


# read data
# datafile = (
#     "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/icrf2_gaiadr2.fits")
datafile = (
    "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/icrf2_gaiadr2_all.fits")


hdulist = fits.open(datafile, memmap=True)
tbdat = hdulist[1].data

# ICRF2 data
# souIVS = tbdat.field("ivs_souname")
# souIVS = tbdat.field("ivs_name")
souIERS_i = tbdat.field("iers_name_icrf2")
RAdeg_i = tbdat.field("ra_icrf2")
DEdeg_i = tbdat.field("dec_icrf2")
DEdeg_err_i = tbdat.field("dec_err_icrf2")


# Gaia DR2 data
souIERS_g = tbdat.field("iers_name_gaia")
RAdeg_g = tbdat.field("ra")
DEdeg_g = tbdat.field("dec")
DEdeg_err_g = tbdat.field("dec_error")

sou_type = tbdat.field("sou_type")

print("# Load data: OK")

# souIERS = np.where(souIERS_i == ' '*8, souIERS_g, souIERS_i)


# # Write data.
# fout = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
#             "gaiadr2_all_solve.src", "w")
# # Write header
# print("$$ This file consists of SOLVE-format positions for 2820 "
#       "sources from GaiaDR2_IERS.\n$$ original file:", datafile, file=fout)
# # Write positions
# write_solvesrc(souIVS, RAdeg, DEdeg, DEdeg_err, "GaiaDR2_IERS", fout)
# # Close file.
# fout.close()

# # 2820 sources in Gaia DR1 and other ICRF2 sources in ICRF2 catalog
# comfile = "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/qso_icrf2.src"
# print(RAdeg)
# for i, RAdegi in enumerate(RAdeg):
#     if np.isnan(RAdegi):
#         print(i)

print(sou_type.size)

souIERS = np.where(souIERS_g == ' '*8, souIERS_i, souIERS_g)
RAdeg = np.where(np.isnan(RAdeg_g), RAdeg_i, RAdeg_g)
DEdeg = np.where(np.isnan(RAdeg_g), DEdeg_i, DEdeg_g)
DEdeg_err = np.where(np.isnan(RAdeg_g), DEdeg_err_i, DEdeg_err_g)

fcom = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
            "gaiadr2_or_icrf2.src", "w")
# Write header
print("$$ This file consists of SOLVE-format positions "
      "for 2820 sources in Gaia DR2 "
      "and other ICRF2 sources in ICRF2 catalog.", file=fcom)
# Write position of sources in Gaia DR2
write_solvesrc(souIERS, RAdeg, DEdeg, DEdeg_err, "GaiaDR2_IERS", fcom)
# # Write position of sources in ICRF2
# write_icrf2_solve(
#     # "/home/nliu/Data/icrf2.dat") # vlbi2
#     "/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat",  # My MacOS
#     souIVS, fcom)
# Close file.
fcom.close()

# # GaiaDR2 or ICRF2 combination
# souIERS = np.where(souIERS_i == ' '*8, souIERS_g, souIERS_i)


# # Source List
# Gaia DR2 all 2820 source
# fall = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
#             "gaiadr2_iers_all.ivs_list", "w")
# print("## List (IVS name) of all 2820 sources in GaiaDR2_IERS", file=fall)
# write_NNRS(souIVS, fall)
# fall.close()


# # ICRF2 defining source
# fdef = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
#             "gaiadr2_iers_def.ivs_list", "w")
# print("## List (IVS name) of all 268 ICRF2 defining sources in GaiaDR2_IERS",
#       file=fdef)
# soudef = souIVS[sou_type == "D"]


# # simple check
# if soudef.size != 268:
#     print("Wrong source number in ICRF2 defining list")
# write_NNRS(soudef, fdef)
# fdef.close()

# # Check data
# data = np.genfromtxt("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/"
#                      "data/GaiaDR2_icrf/qso_solve.src",
#                      usecols=(0,), comments="$$")
# if data.size == souICRF.size:
#     print("# Write successfully. Exit!")
# else:
#     print("# Error! Please check file %s" % outfile)

# print('=======================  END  ============================')

# --------------------------------- END --------------------------------
