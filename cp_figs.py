#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: cp_figs.py
"""
Created on Wed Jan 31 17:10:57 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Copy figures from the directory where it was plotted to the directory
where it can be used in Latex.

"""

# import numpy as np
from os import system

# -----------------------------  FUNCTIONS -----------------------------
indir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/plots"
outdir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/notes/20180130"

figlist = ["opa2018a_00Mnga_dU_dif05.eps",
           "opa2018a_00Mnga_dXp_dif10.eps",
           "opa2018a_00Mnga_dYp_dif10.eps",
           "opa2018a_15Mnga_dU_dif05.eps",
           "opa2018a_15Mnga_dXp_dif10.eps",
           "opa2018a_15Mnga_dYp_dif10.eps",
           "opa2018a_ref15M00_dU_dif05.eps",
           "opa2018a_ref15M00_dXp_dif10.eps",
           "opa2018a_ref15M00_dYp_dif10.eps",
           "opa2018a_00Mnga_ddX_dif30.eps",
           "opa2018a_00Mnga_ddY_dif30.eps",
           "opa2018a_15Mnga_ddX_dif30.eps",
           "opa2018a_15Mnga_ddY_dif30.eps",
           "opa2018a_ref15M00_ddX_dif30.eps",
           "opa2018a_ref15M00_ddY_dif30.eps",
           "opa2018a_00Mnga_sou_dif.eps",
           "opa2018a_15Mnga_sou_dif.eps",
           "opa2018a_ref15M00_sou_dif.eps"]

for fig in figlist:
    system("cp %s/%s %s/" % (indir, fig, outdir))
    print("Copy %s : done!" % fig)
# --------------------------------- END --------------------------------
