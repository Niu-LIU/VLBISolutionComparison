#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: write_sou_solve.py
"""
Created on Thu Apr 26 17:43:41 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from write_solvesrc import write_solvesrc, write_NNRS
# from find_icrf2def import find_icrf2def
from souname_xmatch import read_soun
from read_gaiadr2 import read_gaiadr2_ascii


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

    for (souni, RAhi, RAmi, RAsi, DCdi, DCmi, DCsi, DC_erri
         ) in zip(souivs, RAh, RAm, RAs, DCd, DCm, DCs, DC_err):

        if souni not in Gaiasouivs:
            print(linefmt % (souni, RAhi, RAmi, RAsi,
                             DCdi, DCmi, DCsi, DC_erri, 'icrf2'),
                  file=fout)
        else:
            x += 1


def soudesign_trans(dsg1_sou, dsg1_all, dsg2_all):
    '''source designation1 -> source designation2
    e.g. ICRF designation -> IVS designation.
    '''

    deg2_sou = np.empty_like(dsg1_sou)
    # deg2_sou = []

    for i, sou in enumerate(dsg1_sou):
        if sou in dsg1_all:
            deg2_sou[i] = dsg2_all[dsg1_all == sou][0]
            # deg2_sou.append(dsg2_all[dsg1_all == sou][0])
        else:
            print("# Couldn't find the counterpart of %s" % sou)
            deg2_sou[i] = ' '
            # exit()

    # print("# ICRF designation -> IVS designation: OK")
    # deg2_sou = np.asarray(deg2_sou)

    return deg2_sou


def GaiaDR22solve():
    '''Write the Gaia DR1 source position into the Solve format.
    '''

    print('======================= BEGIN ============================\n'
          '# Write Gaia DR1 qso into Solve format')

    # Load data.
    # ivs, icrf, _ = read_soun(
    #     '/Users/Neo/Astronomy/Data/SOLVE/IVS_SrcNamesTable.txt')
    # ivs, icrf = np.genfromtxt(
    #     '/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat',
    #     usecols=(1, 0), dtype=str, unpack=True)

    ivs, icrf, iers = read_soun()

    _, RAdeg, _, DEdeg, DEdeg_err = read_gaiadr2_ascii(
        '/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/GaiaDR2_icrf/'
        'GaiaDR2_qso.dat')
    # # mas -> uas
    # DEdeg_err = DEdeg_err * 1.e3
    # souIERS = np.genfromtxt("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/"
    #                         "data/GaiaDR2_icrf/singlefile/sourcelist",
    #                         dtype=str)
    souIVS = np.genfromtxt("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/"
                           "data/GaiaDR2_icrf/singlefile/sourcelist",
                           dtype=str)
    print("# Load data: OK")

    # icrf designation -> ivs designation
    # souIVS = soudesign_trans(souIERS, iers, ivs)

    # Write data.
    # outfile = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/GaiaDR2_icrf/qso_solve.src"
    fout = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
                "GaiaDR2_icrf/gaiadr2_all_solve.src", "w")
    # Write header
    print("$$ This file consists of SOLVE-format positions for 2820 "
          "sources from GaiaDR2_IERS .", file=fout)
    # Write positions
    write_solvesrc(souIVS, RAdeg, DEdeg, DEdeg_err, "GaiaDR2_IERS", fout)
    # Close file.
    fout.close()

    # 2820 sources in Gaia DR1 and other ICRF2 sources in ICRF2 catalog
    comfile = "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/qso_icrf2.src"
    fcom = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
                "GaiaDR2_icrf/gaiadr2_icrf2.src", "w")
    # Write header
    print("$$ This file consists of SOLVE-format positions "
          "for 2820 sources in Gaia DR1 "
          "and other ICRF2 sources in ICRF2 catalog.", file=fcom)
    # Write position of sources in Gaia DR1
    write_solvesrc(souIVS, RAdeg, DEdeg, DEdeg_err, "GaiaDR2_IERS", fcom)
    # Write position of sources in ICRF2
    write_icrf2_solve(
        # "/home/nliu/Data/icrf2.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat",  # My MacOS
        souIVS, fcom)
    # Close file.
    fcom.close()

    # Source List
    fall = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
                "GaiaDR2_icrf/qso_all.list", "w")
    print("## List of all 2820 sources in GaiaDR2_IERS", file=fall)
    write_NNRS(souIVS, fall)
    fall.close()

    # # ICRF2 defining source
    # fdef = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
    #             "GaiaDR2_icrf/qso_def.list", "w")
    # print("## List of all 262 ICRF2 defining sources in GaiaDR2_IERS",
    #       file=fdef)
    # soudef = find_icrf2def(souIVS, "IVS")
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


GaiaDR22solve()
# --------------------------------- END --------------------------------
