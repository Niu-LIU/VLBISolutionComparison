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


# -----------------------------  FUNCTIONS -----------------------------
def write_cat():
    pass


# -----------------------------  MAIN  ---------------------------------
# Check input parameters from std
if len(sys.argv) is 2:
    ifile = sys.argv[1]
else:
    print(sys.argv)
    print("Input error!")
    exit()

print("-------------------------- BEGIN ---------------------------")
print("  Remove these sources with 0 observation used in Solve.")
os.system("~/Astronomy/Works/201711_GDR2_ICRF3/progs/sou_elim %s.sou"
          % ifile)

# --------------------------------- END --------------------------------
