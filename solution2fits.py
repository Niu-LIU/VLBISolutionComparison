#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: solution2fits.py
"""
Created on Tue Apr  3 16:22:42 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Read .sou file and rewrite the source position in .fits format.

History
N. Liu, 27/04/2018: Add IVS source name, ICRF designation, and IERS
                    designation;

"""

import numpy as np
import time
import sys
import os
from astropy.io import fits
from read_sou import read_sou_pos
from souname_xmatch import find_sou_designation


# -----------------------------  FUNCTIONS -----------------------------
def sou2fits(ifile, ofile=None):
    '''Transform GETPAR_SOU format into SOLVEsrc format


    Parameters
    ----------
    ifile: string
        full name of .sou file as the input file
    ofile: string
        full name of output file. if None, .fits file will be selected as
        default.

    Returns
    ----------
    None
    '''

    # First eliminate the sources with 0 observations
    os.system("~/Astronomy/Works/201711_GDR2_ICRF3/progs/sou_elim %s"
              % ifile)

    # Read file
    sou, RA, RA_err, DC, DC_err, corr = read_sou_pos(ifile)

    IVS, ICRF, IERS = find_sou_designation(sou, "IVS")

    if ofile is None:
        ofile = "%s.fits" % ifile[:-4]
    # fout = open(ofile, "w")

    print("-------------- BEGIN --------------\n"
          "Now write data into fits format.\n"
          "Original file: %s\n"
          "Output file:   %s" % (ifile, ofile))

    tbhdu = fits.BinTableHDU.from_columns([
        fits.Column(name="ivs_sourcename", format='A8', array=IVS),
        # fits.Column(name="icrf_sourcename", format='A16', array=ICRF),
        fits.Column(name="iers_sourcename", format='A8', array=IERS),
        fits.Column(name='ra_vlbi', format='F16.12', unit='deg', array=RA),
        fits.Column(name='dec_vlbi', format='F16.12', unit='deg', array=DC),
        fits.Column(name='ra_err_vlbi', format='F10.4',
                    unit='mas', array=RA_err),
        fits.Column(name='de_cerr_vlbi', format='F10.4',
                    unit='mas', array=DC_err),
        fits.Column(name='ra_dec_corr_vlbi', format='F6.3', array=corr)])

# header information of FITS
    prihdr = fits.Header()
    prihdr['Creator'] = 'Niu Liu'
    prihdr['COMMENT'] = "Quasar position taken from %s" % ifile
    prihdr['COMMENT'] = "Created time : %s" % time.strftime(
        "%a %b %d %H:%M:%S %Y", time.localtime())
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(ofile)
    print("--------------- END ---------------\n")


def test_code():
    sou2fits("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/"
             "GalacticAberration/opa2018a_nga/opa2018a_nga.sou")


# -------------------------------- MAIN --------------------------------
if len(sys.argv) is 2:
    sou2fits(sys.argv[1])
elif len(sys.argv) is 3:
    sou2fits(sys.argv[1], sys.argv[2])
else:
    exit()

# Test code
# test_code()
# --------------------------------- END --------------------------------
