#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: VLBI_solution_check.py
"""
Created on Thu Jan 11 16:26:47 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import sys
from check_eob_wrt_c04 import check_eob_wrt_c04
from check_nutation_offset import nutation_offset_stat
from solution_icrf2_comparison import solution_icrf2_diff_analysis
from solution_GaiaDR1_comparison import solution_Gaia_diff_analysis


# -----------------------------  FUNCTIONS -----------------------------

# -----------------------------  MAIN  ---------------------------------
# # For opa2018a_ga
# sol_dir = "/home/nliu/solutions/GalacticAberration/opa2018a_ga"
# sol_lab = "opa2018a_ga"

# # For opa2018a_nonga
# sol_dir = "/home/nliu/solutions/GalacticAberration/opa2018a_nonga"
# sol_lab = "opa2018a_nonga"

# Check input parameters
if len(sys.argv) is 3:
    sol_dir, sol_lab = sys.argv[1:3]
else:
    print(sys.argv)
    print("Input error!")
    exit()

print("-------------  BEGIN  ----------------------------")
print("Check the VLBI global solution: %s\n" % sol_lab)

# Check EOP and Nutation offsets.o
print("#1) check .eob file wrt. C04 series:\n"
      "  .eob file: %s.eob" % sol_lab)
check_eob_wrt_c04("%s/%s.eob" % (sol_dir, sol_lab))
print("Check eob: done!\n")

# Check Nutation offsets wrt IAU 2006/2000A Precession-Nutation
print("#2) Check Nutation offsets wrt. IAU 2006/2000A Precession-Nutation:\n"
      "  .eob file: %s.eob" % sol_lab)
nutation_offset_stat("%s/%s.eob" % (sol_dir, sol_lab))
print("Check Nutation: done!\n")

print("#3) Check source position:\n"
      "  .sou file: %s.sou" % sol_lab)
# Check source position wrt. ICRF2 catalog
print("#3.1) Check source position wrt. ICRF2 catalog:")
solution_icrf2_diff_analysis(
    # "%s/%s.sou" %
    sol_dir, "%s.sou" % sol_lab, sol_lab)
# , "%s_icrf2" % sol_lab)

# Check source position wrt. GaiaDR1 catalog
print("#3.2) Check source position wrt. GaiaDR1 catalog:")
# sol_dir, sol_lab, "%s_Gaiadr1" % sol_lab)
solution_Gaia_diff_analysis(sol_dir, "%s.sou" % sol_lab, sol_lab)
print("Check source positions: done!\n")
# --------------------------------- END --------------------------------
