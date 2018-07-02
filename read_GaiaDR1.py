#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_GaiaDR1.py
"""
Created on Mon Apr  9 12:02:56 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np


# -----------------------------  FUNCTIONS -----------------------------
def read_gaiadr1(datafile="/Users/Neo/Astronomy/Data/catalogs/GaiaDR1_cds/"
                 "qso.dat"):
    '''Read Gaia DR1 position.

    Note:   1) e_RA = err(RA*cos(DC))
            2) unit for RA/DC : deg
                    for e_RA/e_DC : mas
    '''

    sourcename, Flag = np.genfromtxt(
        datafile, usecols=(10, 11), dtype=str,
        delimiter="|", unpack=True)

    RAdeg, e_RAmas, DCdeg, e_DCmas, corr = np.genfromtxt(
        datafile, usecols=range(3, 8), delimiter="|", unpack=True)

    return sourcename, RAdeg, DCdeg, e_RAmas, e_DCmas, corr, Flag

# --------------------------------- END --------------------------------
