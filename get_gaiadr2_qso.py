#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: get_gaiadr2_qso.py
"""
Created on Wed Apr 25 18:25:19 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import os


# -----------------------------  FUNCTIONS -----------------------------
# -------------------------------  MAIN --------------------------------
# Read the GaidDR2 ID for AGNs
allwise_ids, gaiadr2_ids = np.genfromtxt(
    "/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/"
    "aux_allwise_agn_gdr2_cross_id.csv",
    dtype=str, unpack=True, delimiter=",")

# print(allwise_ids, gaiadr2_ids)

# The first line is not an ID
allwise_ids = allwise_ids[1:]
gaiadr2_ids = gaiadr2_ids[1:]

# script for download data
downloadpy = "/Users/Neo/Astronomy/Tools/python-cdsclient/find_gaia_dr2.py"

# directory to put data
data_dir = "/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/GaiaDR2_qso/"

# log file
flog = open("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/logs/"
            "get_gaiadr2_qso.log", "w")

# print(allwise_ids.size)
for i, (allwise_id, gaiadr2_id) in enumerate(zip(allwise_ids, gaiadr2_ids)):
    if os.path.exists("%s/%s.votable" % (data_dir, allwise_id)):
        print("#%s has been downloaded!" % allwise_id)
        print("#%s has been downloaded!" % allwise_id, file=flog)
    else:
        os.system("python %s --format votable --source_id=%s >%s/%s.votable"
                  % (downloadpy, gaiadr2_id, data_dir, allwise_id))
        print(" Download data for %s, count=%10d" % (allwise_id, i))
        print(" Download data for %s" % allwise_id, file=flog)

flog.close()

print("Done!")
# --------------------------------- END --------------------------------
