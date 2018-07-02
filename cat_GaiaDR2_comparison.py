#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: cat_GaiaDR2_comparison.py
"""
Created on Tue May 29 16:32:28 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import sys
import os


# -----------------------------  FUNCTIONS -----------------------------

# -----------------------------  MAIN  ---------------------------------
# Check input parameters
if len(sys.argv) is 3:
    sol_dir, sol_lab = sys.argv[1:3]
else:
    print(sys.argv)
    print("Input error!")
    exit()
# --------------------------------- END --------------------------------
