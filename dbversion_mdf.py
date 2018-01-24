#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: dbversion_mdf.py
"""
Created on Fri Dec 22 10:33:16 2017

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np


# -----------------------------  FUNCTIONS -----------------------------
supfile = "/home/nliu/solutions/test/super.cat"
arcfile = "/home/nliu/solutions/test/arc.2017a"
arnfile = "/home/nliu/solutions/test/arc.2017a_new"
logfile = "/home/nliu/solutions/test/log"

dbl, ver = np.genfromtxt(supfile, usecols=(0, 1),
                         dtype=str, unpack=True, comments="*")

farn = open(arnfile, "w")
flog = open(logfile, "w")

for line in open(arcfile, "r"):
    if line[0] is not "*":
        line = line[:-1]
        dbli, veri = line[3:12], line[15:17]
        if dbli[-1] is " ":
            dbli = dbli[:-1]
        if veri[0] is " ":
            veri = veri[-1]

        ind = (dbl == dbli)
        vern = ver[ind]

        if vern:
            vern = vern[0]
            if veri == vern:
                print(line, file=farn)
            else:
                print("Modify the version of %s: %s --> %s" %
                      (dbli, veri, vern))
                print("%s%2s%s" % (line[:15], vern, line[17:]), file=farn)
        else:
            print("Please check this data %s" % dbli,
                  file=flog)

farn.close()
flog.close()

# --------------------------------- END --------------------------------
