#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: VLBI_solutions_comparison.py
"""
Created on Thu Jan 11 18:31:49 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import sys
from eob_comparison import eob_comparison
from source_position_comparison import sou_position_com


# -----------------------------  FUNCTIONS -----------------------------
# # For opa2018a_ga
# sol_dir1 = "/home/nliu/solutions/GalacticAberration/opa2018a_ga"
# sol_lab1 = "opa2018a_ga"
# # For opa2018a_nonga
# sol_dir2 = "/home/nliu/solutions/GalacticAberration/opa2018a_nonga"
# sol_lab2 = "opa2018a_nonga"
# com_label = "opa2018a_gaMnon"

# Check input parameters
if len(sys.argv) is 6:
    sol_dir1, sol_lab1, sol_dir2, sol_lab2, com_label = sys.argv[1:6]
else:
    print("Input error!")
    exit()

print("-------------  BEGIN  ----------------------------")
print("Comparison between the VLBI global solutions: %s and %s \n" %
      (sol_lab1, sol_lab2))

# Compare EOP and Nutation offsets.
print("#1) Compare .eob files:\n"
      "  .eob file: %s.eob/%s.eob" % (sol_lab1, sol_lab2))
eob_comparison("%s/%s.eob" % (sol_dir1, sol_lab1),
               "%s/%s.eob" % (sol_dir2, sol_lab2),
               com_label)
print("Compare eob: done!\n")

# Compare source positions
print("#2) Compare source positions:\n"
      "  .sou file: %s.sou/%s.sou" % (sol_lab1, sol_lab2))
sou_position_com("%s/%s.sou" % (sol_dir1, sol_lab1),
                 "%s/%s.sou" % (sol_dir2, sol_lab2),
                 com_label)

# # Check source position wrt. GaiaDR1 catalog
# print("#3) Check source position wrt. GaiaDR1 catalog:\n"
#       "  .sou file: %s.sou" % sol_lab)
# solution_Gaia_diff_analysis("%s/%s.sou" % (sol_dir, sol_lab))

# --------------------------------- END --------------------------------
