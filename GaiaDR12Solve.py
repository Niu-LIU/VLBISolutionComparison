#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: format_trans.py
"""
Created on Mon Dec 18 10:27:19 2017

@author: Neo(liuniu@smail.nju.edu.cn)

History:

Who         When            What
N., Liu     22 Jan, 2018    Add function 'GaiaDR12solve', and
                            change function 'write_solvesrc' for later
                            use, e.g., Gaia DR2.

"""

import numpy as np
from write_solvesrc import write_solvesrc, write_NNRS
from find_icrf2def import find_icrf2def
from souname_xmatch import read_soun


__all__ = {"read_GaiaDR1qso", "read_soun" "soudesign_trans",
           "GaiaDR12solve"}


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


def read_GaiaDR1qso(datafile):
    '''Read Gaia DR1 qso data.
    '''
    souicrf = np.genfromtxt(
        datafile, dtype=str, usecols=(10,), delimiter='|')
    RAdeg, DEdeg, DEdeg_err = np.genfromtxt(
        datafile, usecols=(3, 5, 6), unpack=True, delimiter='|')

    return souicrf, RAdeg, DEdeg, DEdeg_err


# def read_soun(datafile):
#     IVS, ICRF, IERS = np.genfromtxt(
#         datafile, usecols=(0, 1, 3), dtype=str, unpack=True)

#     return IVS, ICRF, IERS


def soudesign_trans(dsg1_sou, dsg1_all, dsg2_all):
    '''source designation1 -> source designation2
    e.g. ICRF designation -> IVS designation.
    '''

    deg2_sou = np.empty_like(dsg1_sou)

    for i, sou in enumerate(dsg1_sou):
        if sou in dsg1_all:
            deg2_sou[i] = dsg2_all[dsg1_all == sou][0]
        else:
            print("# Couldn't find the counterpart of %s" % sou)
            exit()

    # print("# ICRF designation -> IVS designation: OK")

    return deg2_sou


def GaiaDR12solve():
    '''Write the Gaia DR1 source position into the Solve format.
    '''

    print('======================= BEGIN ============================\n'
          '# Write Gaia DR1 qso into Solve format')

    # Load data.
    # ivs, icrf, _ = read_soun(
    #     '/Users/Neo/Astronomy/Data/SOLVE/IVS_SrcNamesTable.txt')
    ivs, icrf = np.genfromtxt(
        '/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat',
        usecols=(1, 0), dtype=str, unpack=True)
    souICRF, RAdeg, DEdeg, DEdeg_err = read_GaiaDR1qso(
        '/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/qso.dat')
    print("# Load data: OK")

    # icrf designation -> ivs designation
    souIVS = soudesign_trans(souICRF, icrf, ivs)

    # Write data.
    outfile = "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/qso_solve.src"
    fout = open(outfile, "w")
    # Write header
    print("$$ This file consists of SOLVE-format positions for 2191 "
          "sources from GaiaDR1 qso.", file=fout)
    # Write positions
    write_solvesrc(souIVS, RAdeg, DEdeg, DEdeg_err, "Gaia DR1 qso", fout)
    # Close file.
    fout.close()

    # 2191 sources in Gaia DR1 and other ICRF2 sources in ICRF2 catalog
    comfile = "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/qso_icrf2.src"
    fcom = open(comfile, "w")
    # Write header
    print("$$ This file consists of SOLVE-format positions "
          "for 2191 sources in Gaia DR1 "
          "and other ICRF2 sources in ICRF2 catalog.", file=fcom)
    # Write position of sources in Gaia DR1
    write_solvesrc(souIVS, RAdeg, DEdeg, DEdeg_err, "Gaia DR1 qso", fcom)
    # Write position of sources in ICRF2
    write_icrf2_solve(
        # "/home/nliu/Data/icrf2.dat") # vlbi2
        "/Users/Neo/Astronomy/Data/catalogs/icrf/icrf2.dat",  # My MacOS
        souIVS, fcom)
    # Close file.
    fcom.close()

    # Source List
    allsoufile = "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/qso_all.list"
    fall = open(allsoufile, "w")
    print("## List of all 2191 sources in Gaia DR1 qso", file=fall)
    write_NNRS(souIVS, fall)
    fall.close()

    # ICRF2 defining source
    defsoufile = "/Users/Neo/Astronomy/Data/catalogs/Gaia_cds/qso_def.list"
    fdef = open(defsoufile, "w")
    print("## List of all 262 ICRF2 defining sources in Gaia DR1 qso",
          file=fdef)
    soudef = find_icrf2def(souIVS, "IVS")
    write_NNRS(soudef, fdef)
    fdef.close()

    # Check data
    data = np.genfromtxt(outfile, usecols=(0,), comments="$$")
    if data.size == souICRF.size:
        print("# Write successfully. Exit!")
    else:
        print("# Error! Please check file %s" % outfile)

    print('=======================  END  ============================')


# ----------------------------------------------------------------------
GaiaDR12solve()
# --------------------------------- END --------------------------------
