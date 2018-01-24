#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: verify_codes.py
"""
Created on Tue Jan  9 15:26:08 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from scipy import stats
from linear_regression import linear_regression


# -----------------------------  FUNCTIONS -----------------------------
x = np.random.normal(0, 1, 100)
y = 1.5 * x + 0.4 + np.random.normal(0, 0.5, 100)
# err = np.ones_like(x)
err = np.random.normal(0, 0.5, 100)
pi = err**-2

slope, intercept, r_value, p_value, std_err = stats.linregress(x, y * pi)
print(slope, intercept)

par, err, outlier, cor = linear_regression(x, y, err)
print(par)
# --------------------------------- END --------------------------------
