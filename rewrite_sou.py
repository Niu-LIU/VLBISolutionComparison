#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: rewrite_sou.py
"""
Created on Thu Apr 26 15:46:01 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import os
import sys
import time
from read_sou import read_sou
from souname_xmatch import find_sou_designation

__all__ = ["write_cat", "rewrite_sou"]


# -----------------------------  FUNCTIONS -----------------------------
def write_cat(ifile, ofile=None):
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

    # Read file
    [sou, RA, RA_err, DC, DC_err, cor,
     ObsU, _, SesU, _, ObsB, ObsE] = read_sou(ifile,
                                              unit_deg=True, arcerr=True)
    ObsM = (ObsB + ObsE) / 2.

    IVS, _, IERS = find_sou_designation(sou, "IVS")

    fout = open(ofile, "w")
    # print header
    print("# VLBI Celestial Reference Frame Solution OPA2018b\n"
          "#\n"
          "# Columns  Units   Meaning\n"
          "#    1     --      IVS designation\n"
          "#    2     --      IERS designation\n"
          "#    3     deg     Right ascension\n"
          "#    4     deg     Declination\n"
          "#    5     mas     Formal uncertainty of the right ascension "
          "(*cos(Dec))\n"
          "#    6     mas     Formal uncertainty of the declination\n"
          "#    7     --      Correlation between right ascension and "
          "declination\n"
          "#    8     days    Average epoch of observation\n"
          "#    9     days    First epoch of observation\n"
          "#   10     days    Last epoch of observation\n"
          "#   11     --      Number of sessions\n"
          "#   12     --      Number of delays\n"
          "#   13     --      Number of delay rates\n"
          "#   14     --      Estimation flag (global or local)\n"
          "# Created date: %s\n#"
          % time.strftime("%d/%m/%Y", time.localtime()), file=fout)

    for (IVSi, IERSi, RAi, DCi, RA_erri, DC_erri, cori,
         ObsMi, ObsBi, ObsEi, ObsUi, SesUi) in zip(
            IVS, IERS, RA, DC, RA_err, DC_err, cor,
            ObsM, ObsB, ObsE, ObsU, SesU):

        if IERSi == " " * 8:
            print("%-8s  ********  "
                  "%14.10f  %+14.10f  %10.4f  %10.4f  %+10.4f  "
                  "%7.2f  %7.2f  %7.2f  %5d  %10d  0  GLO"
                  % (IVSi, RAi, DCi, RA_erri, DC_erri, cori,
                     ObsMi, ObsBi, ObsEi, SesUi, ObsUi), file=fout)
        else:
            print("%-8s  %8s  %14.10f  %+14.10f  %10.4f  %10.4f  %+10.4f  "
                  "%7.2f  %7.2f  %7.2f  %5d  %10d  0  GLO"
                  % (IVSi, IERSi, RAi, DCi, RA_erri, DC_erri, cori,
                     ObsMi, ObsBi, ObsEi, SesUi, ObsUi), file=fout)

    fout.close()


def rewrite_sou(ifile, ofile=None):
    '''convert .sou format


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

    if ofile is None:
        ofile = "%s.cat" % ifile[:-4]

    if os.path.exists(ofile):
        os.system("rm %s" % ofile)

    print("  Remove these sources with 0 observation used in Solve.")
    os.system("~/Astronomy/Works/201711_GDR2_ICRF3/progs/sou_elim %s"
              % ifile)

    write_cat(ifile, ofile)
    print("   %s --> %s: done!\n" % (ifile, ofile))


# -----------------------------  MAIN  ---------------------------------
# print("-------------------------- BEGIN ---------------------------")
# # Check input parameters from std
# if len(sys.argv) is 2:
#     ifile = sys.argv[1]
#     rewrite_sou(ifile)
# elif len(sys.argv) is 3:
#     ifile, ofile = sys.argv[1:3]
#     rewrite_sou(ifile, ofile)
# else:
#     print(sys.argv)
#     print("Input error!")
#     exit()
# --------------------------------- END --------------------------------
