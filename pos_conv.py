#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:14:30 2017

@author: Neo

Convert the string of RA/DEC into float.

"""


# ------------------------------  FUNCTIONS  ---------------------------
def RA_conv(RAbyt):
    '''Convert right ascesion string of HH_MM_SS.ssssssss into float.
    '''
    RAstr = str(RAbyt, encoding="utf8")
    hours, mins, secs = RAstr.split('_')
    RAs = (float(hours) * 60 + float(mins)) * 60 + float(secs)
    # second -> arcsec
    return RAs * 15


def DC_conv(DCbyt):
    '''Convert declination string of +DD_AM_AS.ssssssss into float.
    '''
    DCstr = str(DCbyt, encoding="utf8")
    degs, ams, ass = DCstr.split('_')
    # Determine the sign
    if DCstr[0] is '+':
        return (float(degs) * 60 + float(ams)) * 60 + float(ass)
    else:
        return (float(degs) * 60 - float(ams)) * 60 - float(ass)
# ------------------------------ END -----------------------------------
