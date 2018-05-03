#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: souname_xmatch.py
"""
Created on Thu Apr 26 15:55:00 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Read radio source name from data list
"""

from astropy.io import fits
import numpy as np
import os
import time
from cross_match import list_crossmatch


# -----------------------------  FUNCTIONS -----------------------------
def rewrite_soun():
    '''Read source names.

    Parameters
    ----------
    None

    Returns
    ----------
    None
    '''

    ifile = open("/Users/Neo/Astronomy/Data/SOLVE/IVS_SrcNamesTable.txt",
                 "r")
    ofile = open("/Users/Neo/Astronomy/Data/SOLVE/IVS_SrcNamesTable.csv",
                 "w")

    print("#==================================\n"
          "# IVS source Name Translation Table\n"
          "#==================================\n"
          "  # File format      : 26/04/2018\n"
          "  # Date last revised: %s\n"
          "# Columns      Format             Content\n"
          "#   1          A8                 IVS source name\n"
          "#                                 Source name used "
          "now or in the past in IVS schedules and databases\n"
          "#   2          A16                J2000 source name long\n"
          "#                                 ICRF designations following"
          " therecommendations of the IAU Task Group on Designations\n"
          "#   3          A8                 IERS designations source name"
          " constructed from B1950 coordinates\n"
          % time.strftime("%d/%m/%Y", time.localtime()), file=ofile)

    for line in ifile.readlines():
        if not len(line) or line.startswith('#'):
            continue

        if line[40] is "-" or line[40] is " ":
            print("%s,%s,%s" % (line[:8], line[10:26], line[:8]),
                  file=ofile)
        else:
            print("%s,%s,%s" % (line[:8], line[10:26], line[40:48]),
                  file=ofile)

    ofile.close()


def read_soun():
    '''Read source names.

    Parameters
    ----------
    None

    Returns
    -------
    IVS : array of string
        IVS source name
    ICRF : array of string
        J2000 source name long
    IERS : array of string
        IERS designations source name
    '''

    datafile = "/Users/Neo/Astronomy/Data/SOLVE/IVS_SrcNamesTable.csv"

    if not os.path.exists(datafile):
        rewrite_soun()

    data_file = open(datafile, "r")

    # empty list ot store data
    IVS, ICRF, IERS = np.genfromtxt(datafile, dtype=str,
                                    delimiter=",", unpack=True)

    return IVS, ICRF, IERS


def rewrite_sourcename(datafile=None):
    '''Read source names from a file named 'source.names' and rewrite it.

    Parameters
    ----------
    datafile : str
        full path of source.names. if None, default value will be used.

    Returns
    -------
    IVS : array of string
        IVS source name
    ICRF : array of string
        J2000 source name long
    IERS : array of string
        IERS designations source name
    '''

    if datafile is None:
        datafile = ("/Users/Neo/Astronomy/Tools/Calc_Solve/mk5-opa/"
                    "save_files/source.names")

    # print(datafile)
    ifile = open(datafile, "r")

    ivs_name = []
    iers_name = []
    Typ = []

    for line in ifile.readlines():
        if not len(line) or line.startswith('#') or line.startswith('\n'):
            continue

        # print(line[0], len(line))
        ivs_name.append(line[:8])
        iers_name.append(line[32:40])
        Typ.append(line[42])

    ifile.close()

    # convert list to array
    ivs_name = np.asarray(ivs_name)
    iers_name = np.asarray(iers_name)
    Typ = np.asarray(Typ)

    print(ivs_name, iers_name, Typ, Typ.size)

    # if not os.path.exists(datafile):
    #     rewrite_soun()

    # data_file = open(datafile, "r")

    # # empty list ot store data
    # IVS, ICRF, IERS = np.genfromtxt(datafile, dtype=str,
    #                                 delimiter=",", unpack=True)

    # return IVS, ICRF, IERS


def find_sou_designation(sou, nametype):
    """Find the other corresponding radio source designations.

    Parameters
    ----------
    sou : array_like
        radio source name
    nametype : string
        type of input source name, could be set as "IVS", "ICRF", or
        "IERS"

    Returns
    -------
    list1, list2, list3 : ndarrays
        IVS, ICRF, and IERS designation of input sources.
    """

    IVS, ICRF, IERS = read_soun()

    # if nametype is "IVS":
    #     list0 = np.where(IVS in sou, IVS)
    #     list1 = np.where(IVS in sou, ICRF)
    #     list2 = np.where(IVS in sou, IERS)
    # elif nametype is "ICRF":
    #     list0 = np.where(ICRF in sou, IVS)
    #     list1 = np.where(ICRF in sou, ICRF)
    #     list2 = np.where(ICRF in sou, IERS)
    # elif nametype is "IERS":
    #     list0 = np.where(IERS in sou, IVS)
    #     list1 = np.where(IERS in sou, ICRF)
    #     list2 = np.where(IERS in sou, IERS)
    # else:
    #     print("ERROR! nametype can only be IVS, ICRF, or IERS!")
    #     exit()

    if nametype is "IVS" or "ivs":
        list0 = IVS
    elif nametype is "ICRF" or "icrf":
        list0 = ICRF
    elif nametype is "IERS" or "iers":
        list0 = IERS
    else:
        print("ERROR! nametype can only be IVS, ICRF, or IERS!")
        exit()

    list1 = np.empty_like(sou, dtype='S8')
    list2 = np.empty_like(sou, dtype='S16')
    list3 = np.empty_like(sou, dtype='S8')

    for i, soui in enumerate(sou):
        indarr = np.where(list0 == soui)[0]

        if indarr:
            j = indarr[0]
            list1[i] = IVS[j]
            list2[i] = ICRF[j]
            list3[i] = IERS[j]
        else:
            list1[i] = " " * 8
            list2[i] = " " * 16
            list3[i] = " " * 8

    # print(list2)

    return list1, list2, list3


def test():
    """
    """

    x1 = np.array(['a', '4', '3', 'd', 't'])
    x2 = np.array(['4', 'e', 'a', 'd', '5'])

    a = np.arange(5)

    y1 = np.where(x1 in x2, x2, a)

    print(y1)


# test()
# ivs, icrf, iers = read_soun()

# print(ivs[:10])
# print(icrf[:10])
# print(iers[:100])

# read_soun('/Users/Neo/Astronomy/Data/SOLVE/IVS_SrcNamesTable.txt')

# rewrite_sourcename()
# --------------------------------- END --------------------------------
