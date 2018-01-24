#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:00:53 2017

@author: Neo
"""

import os


def ensure_dir(f):
    '''If the directory does not exit, creat it.
    '''
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.mkdirs(d)
# ------------------------------ END -----------------------------------
