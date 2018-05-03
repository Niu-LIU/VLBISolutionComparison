#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: get_gaiadr2_icrf.py
"""
Created on Wed Apr 25 18:44:32 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import os

# -----------------------------  FUNCTIONS -----------------------------
# -------------------------------  MAIN --------------------------------
# Read the GaidDR2 ID for AGNs
iers_ids, gaiadr2_ids = np.genfromtxt(
    "/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/"
    "aux_iers_gdr2_cross_id.csv",
    dtype=str, unpack=True, delimiter=",")

# print(iers_ids, gaiadr2_ids)

# The first line is not an ID
iers_ids = iers_ids[1:]
gaiadr2_ids = gaiadr2_ids[1:]

# script for download data
downloadpy = "/Users/Neo/Astronomy/Tools/python-cdsclient/find_gaia_dr2.py"

# directory to put data
data_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/GaiaDR2_icrf"

# log file
flog = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/logs/"
            "get_gaiadr2_icrf.log", "w")

# print(iers_ids.size)
for i, (iers_id, gaiadr2_id) in enumerate(zip(iers_ids, gaiadr2_ids)):
    if os.path.exists("%s/%s.votable" % (data_dir, iers_id)):
        print("#%s has been downloaded!" % iers_id)
        print("#%s has been downloaded!" % iers_id, file=flog)
    else:
        os.system("python %s --format votable --source_id=%s >%s/%s.votable" %
                  (downloadpy, gaiadr2_id, data_dir, iers_id))
        print(" Download data for %s, count=%4d" % (iers_id, i))
        print(" Download data for %s" % iers_id, file=flog)

flog.close()

print("Done!")
# --------------------------------- END --------------------------------
