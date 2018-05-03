#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: Convert_GaiaDR1.py
"""
Created on Sun Apr 29 16:00:35 2018

@author: Neo(liuniu@smail.nju.edu.cn)

COnvert GaiaDR1 data from ascii to .fits

"""


# -------------------------------  BEGIN -------------------------------
import numpy as np
from astropy.io import fits

# Load data
data_gaiadr1 = ("/Users/Neo/Astronomy/Data/catalogs/GaiaDR1_cds/"
                "gaiadr1_icrf2_1.dat")

icrf_name, ivs_name, iers_name, Typ = np.genfromtxt(
    data_gaiadr1, usecols=(0, 1, 2, 13), delimiter="|",
    dtype=str, unpack=True)
[ra_gaiadr1, dec_gaiadr1,
 ra_err_gaiadr1, dec_err_gaiadr1, ra_dec_corr_gaiadr1,
 ra_icrf2, dec_icrf2,
 ra_err_icrf2, dec_err_icrf2, ra_dec_corr_icrf2] = np.genfromtxt(
    data_gaiadr1, usecols=range(3, 13), delimiter="|", unpack=True)

# uas -> mas
ra_err_gaiadr1 *= 1.e-3
dec_err_gaiadr1 *= 1.e-3
ra_err_icrf2 *= 1.e-3
dec_err_icrf2 *= 1.e-3

# Creat a new fits file to store these data.
outfile = "%s.fits" % data_gaiadr1[:-4]

print("-------------- BEGIN --------------\n"
      "Now write data into fits format.\n"
      "Original file: %s\n"
      "Output file:   %s" % (data_gaiadr1, outfile))

tbhdu = fits.BinTableHDU.from_columns([
    fits.Column(name="icrf_name", format='A16', array=icrf_name),
    fits.Column(name="ivs_souname", format='A8', array=ivs_name),
    fits.Column(name="iers_name", format='A8', array=iers_name),
    fits.Column(name="ra_gaiadr1", format='F16.12',
                unit='deg', array=ra_gaiadr1),
    fits.Column(name='dec_gaiadr1', format='F16.12',
                unit='deg', array=dec_gaiadr1),
    fits.Column(name='ra_err_gaiadr1', format='F10.4',
                unit='mas', array=ra_err_gaiadr1),
    fits.Column(name='dec_err_gaiadr1', format='F10.4',
                unit='mas', array=dec_err_gaiadr1),
    fits.Column(name='ra_dec_corr_gaiadr1', format='F6.3',
                array=ra_dec_corr_gaiadr1),
    fits.Column(name="ra_icrf2", format='F16.12',
                unit='deg', array=ra_icrf2),
    fits.Column(name='dec_icrf2', format='F16.12',
                unit='deg', array=dec_icrf2),
    fits.Column(name='ra_err_icrf2', format='F10.4',
                unit='mas', array=ra_err_icrf2),
    fits.Column(name='dec_err_icrf2', format='F10.4',
                unit='mas', array=dec_err_icrf2),
    fits.Column(name='ra_dec_corr_icrf2', format='F6.3',
                array=ra_dec_corr_icrf2)
])

# header information of FITS
prihdr = fits.Header()
prihdr['Creator'] = 'Niu Liu'
prihdr['COMMENT'] = ("Cross-match result between Gaia DR1 and ICRF2 "
                     "catalogs.")
prihdr['COMMENT'] = ("Original file: %s" % data_gaiadr1)
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto(outfile)
print("--------------- END ---------------\n")
# --------------------------------- END --------------------------------
