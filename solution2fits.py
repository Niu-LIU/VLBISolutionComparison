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

    sou, RAas, RA_err, DCas, DC_err, cor = read_sou_pos(ifile)

    IVS, ICRF, IERS = find_sou_designation(sou, "IVS")

    # as -> deg
    RA = RAas / 3.6e3
    DC = DCas / 3.6e3

    RA_err = RA_err * np.cos(np.deg2rad(DC))

    if ofile is None:
        ofile = "%s.fits" % ifile[:-4]
    # fout = open(ofile, "w")

    print("-------------- BEGIN --------------\n"
          "Now write data into fits format.\n"
          "Original file: %s\n"
          "Output file:   %s" % (ifile, ofile))

    tbhdu = fits.BinTableHDU.from_columns([
        fits.Column(name="souname_IVS", format='A8', array=IVS),
        fits.Column(name="souname_ICRF", format='A16', array=ICRF),
        fits.Column(name="souname_IERS", format='A8', array=IERS),
        fits.Column(name='ra', format='F16.12', unit='deg', array=RA),
        fits.Column(name='dec', format='F16.12', unit='deg', array=DC),
        fits.Column(name='ra_err', format='F10.4', unit='mas', array=RA_err),
        fits.Column(name='de_cerr', format='F10.4', unit='mas', array=DC_err),
        fits.Column(name='ra_dec_corr', format='F6.3', array=cor)])

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
    solufits(sys.argv[1])
elif len(sys.argv) is 3:
    solufits(sys.argv[1], sys.argv[2])
else:
    exit()
# Test code
# test_code()
# --------------------------------- END --------------------------------
