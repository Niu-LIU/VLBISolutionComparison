#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: solution2solve.py
"""
Created on Wed Jan 24 15:08:21 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import sys
from os import system

# -----------------------------  FUNCTIONS -----------------------------


def solu2solve(datafile, label):
    '''Transform GETPAR_SOU format into SOLVEsrc format
    '''

    sou, RA, DC, DCerr = np.genfromtxt(datafile, dtype=str,
                                       usecols=(1, 3, 7, 9),
                                       unpack=True)

    outfile = "%s.src" % datafile[:-4]
    fout = open(outfile, "w")

    fmt = "    %8s  %17s  %18s  %9s    ! %s"
    for soui, RAi, DCi, DCerri in zip(sou, RA, DC, DCerr):
        print(fmt % (soui, RAi, DCi, DCerri, label), file=fout)

    fout.close()

    system("sed -i 's/_/ /g' %s" % outfile)


def test_code():
    solu2solve("/obs/nliu/solutions/GalacticAberration/opa2018a_ga/opa2018a_ga.sou",
               "opa2018a")

if len(sys.argv) is 3:
    solu2solve(sys.argv[1], sys.argv[2])
else:
    exit()
# --------------------------------- END --------------------------------
