
# coding: utf-8

# In[2]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: find_icrf2def.py
"""
Created on Thu Jun 21 09:30:04 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
# My module
from read_icrf2 import read_icrf2


def find_icrf2def(soulist, design="IVS"):
    """Find the ICRF2 defining sources in a source list.

    Parameters
    ----------
    soulist : list or array
        source name
    design : str
        source name designation type, 'ivs', 'iers', or 'icrf'

    Returns
    -------
    defsou : list or array
        ICRF2 defining sources in the source list
    """

    # ICRF2 data
    icrfn, ivsn, iersn, _, _, _, _, _, _, Flag = read_icrf2()

    if design is 'ivs' or design is 'IVS':
        icrf2_list = ivsn
    elif design is 'iers' or design is 'IERS':
        icrf2_list = iersn
    elif design is 'icrf' or design is 'ICRF':
        icrf2_list = icrfn

    # ICRF2 295 defining sources
    con = (Flag == "D")
    icrf2_def = icrf2_list[con]

    # Find defining sources
    defsou = []
    for sou in soulist:
        if sou in icrf2_def:
            defsou.append(sou)

    defsou = np.asarray(defsou)

    return defsou
