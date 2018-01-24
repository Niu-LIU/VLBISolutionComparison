#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:19:27 2017

@author: Neo
"""

import numpy as np


# ------------------------------  FUNCTIONS  ---------------------------
def get_dat(souname, soulist, data):
    '''Given a source name, return the corresponding data.
    '''
    return data[np.where(soulist == souname)]


def get_sou_pv(souname, soulist, RA, DC, RA_err, DC_err, cor):
    '''Given a source name, return the 2-D positions / proper-motions.
    '''
    RAsou = get_dat(souname, soulist, RA)
    DCsou = get_dat(souname, soulist, DC)
    RAsou_err = get_dat(souname, soulist, RA_err)
    DCsou_err = get_dat(souname, soulist, DC_err)
    corsou = get_dat(souname, soulist, cor)

    return RAsou, DCsou, RAsou_err, DCsou_err, corsou


def get_sta_pv(staname, stalist, X, Y, Z, X_err, Y_err, Z_err):
    '''Given a station name, return the 3-D positions / velocities.
    '''
    Xsta = get_dat(staname, stalist, X)
    Ysta = get_dat(staname, stalist, Y)
    Zsta = get_dat(staname, stalist, Z)
    Xsta_err = get_dat(staname, stalist, X_err)
    Ysta_err = get_dat(staname, stalist, Y_err)
    Zsta_err = get_dat(staname, stalist, Z_err)

    return Xsta, Ysta, Zsta, Xsta_err, Ysta_err, Zsta_err
# ------------------------------ END -----------------------------------
