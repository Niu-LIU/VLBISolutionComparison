#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: getpar_soupm.py
"""
Created on Tue Jun  5 16:26:25 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np


# -----------------------------  FUNCTIONS -----------------------------
# flag
find_pm = False
line_num = 0
line_data_end = 0

# number of lines to be read every time
no_line_read = 100000

# output
fout = open("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/"
            "ts_JC/180510-pm/180510-pm.spm", "w")

# read file line by line
# with open("180510-pm.spl", "r") as spl_file:
#     while True:
#         lines = spl_file.readlines(no_line_read)
#         if not lines:
#             break
#         for i, line in enumerate(lines):
#             if "RT. ASC.velocity" in line:
#                 print(lines[i:i+3], file=fout)
#                 temp = i
#             elif "DEC.   .velocity" in line:
#                 if i >= temp + 1:
#                     print(lines[i:i+2], file=fout)
#             pass


with open("/Users/Neo/Astronomy/Data/VLBISolutions/vlbi2_server/"
          "ts_JC/180510-pm/180510-pm.spl", "r") as spl_file:
    while True:
        line = spl_file.readline()
        line_found = False
        line_num += 1

        if not line:
            break

        if line.strip():
            if line_num <= line_data_end:
                line_found = True
            elif "RT. ASC.velocity" in line:
                find_pm = True
                line_found = True
                line_data_end = line_num + 2
            elif "DEC.    velocity" in line:
                find_pm = True
                line_found = True
                line_data_end = line_num + 1

            if line_found:
                print(line, end="", file=fout)
            elif find_pm:
                break


# close file handle
fout.close()
# --------------------------------- END --------------------------------
