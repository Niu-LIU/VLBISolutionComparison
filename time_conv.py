#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 15:36:25 2017

@author: Neo
"""

from datetime import datetime


# ------------------------------  FUNCTIONS  ---------------------------
def count_date_num(y1, m1, d1, y2, m2, d2):
    tdelta = datetime(y1, m1, d1) - datetime(y2, m2, d2)
    return tdelta.days


def date2year(s):
    '''Turn a string of date into a float in unit of year.
    '''
    if s == ' ' * 10:
        return 0
    else:
        yyyy, mm, dd = s.split('.')
        Y, M, D = int(yyyy), int(mm), int(dd)
        dateno = count_date_num(Y, M, D, Y, 1, 1)
        dateyr = count_date_num(Y + 1, 1, 1, Y, 1, 1)
        return Y + dateno * 1.0 / dateyr
# ------------------------------ END -----------------------------------
