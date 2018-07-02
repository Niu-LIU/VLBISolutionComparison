#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: cnt_batch.py
"""
Created on Tue May 15 15:37:30 2018

@author: Neo(liuniu@smail.nju.edu.cn)


generate the control files for applying different GA values in VLBI solutions

"""

import numpy as np
import os


# -----------------------------  FUNCTIONS -----------------------------
def cntfile_gen(xmin, xmax, xstep):
    """Generate new control files

    """

    for ga in np.arange(xmin, xmax, xstep):
        new_file = "GA%.1f_15.5.cnt" % ga
        # copy_cmd = "cp GA3.0_15.5.cnt %s" % new_file
        # modi_cmd = "sed -i '2s/3.0/%.1f' %s" % (ga, new_file)

        os.system("cp GA3.0_15.5.cnt %s" % new_file)
        os.system("sed -i '2s/3.0/%.1f/' %s" % (ga, new_file))
        os.system("sed -i '189s/3.0/%.1f/' %s" % (ga, new_file))


# cntfile_gen(3.5, 10.0, 0.5)
cntfile_gen(5.1, 5.5, 0.1)
cntfile_gen(5.6, 6.0, 0.1)
cntfile_gen(6.1, 6.5, 0.1)
cntfile_gen(6.6, 7.0, 0.1)
# --------------------------------- END --------------------------------
