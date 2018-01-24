#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:52:21 2017

@author: Neo

Retrieve the estimates of positions of global stations and the formal
uncertainties of these estimates from .sta file which is generated by
the program getpar.


   .sta file contains estimates of positions of global stations and the formal
uncertainties of these estimates. The list of station positions is sorted in
alphabetic order of station names. Stations before and after episodic motions
are treated as different stations. Correlations between station positions and
velocities are also written.

   File contains lines of four types:

1) Comment. The first character is #. Header comment contain the full name of
   the spool file.

2) Cartesian components of the vector of station position. The first
   8 characters of this line are STA_GCX:

   Field   Format Units Meaning
   1-8     A8     --    record type identifier: STA_GCX:
   11-25   A15    --    station name. Station name consist of 8-letters station
                        acronym and 6-letter epoch in format yymmdd. Epoch
                        is attached to the name only if the station had episodic
                        motion. Fields between the last letter of the station
                        name and the first letter of epoch are filled by _.
                        If the station didn't have episodic name then the space
                        after the last letter of the station name is left blank.
   28-29   A2     --    component identifier. One of "X:", "Y:" or "Z:"
   31-45   F15.2  mm    value of X-component of station position.
   50-59   F10.3  mm    formal uncertainty of X-component of station position.
   65-79   F15.2  mm    value of Y-component of station position.
   84-93   F10.3  mm    formal uncertainty of Y-component of station position.
   99-113  F15.2  mm    value of Z-component of station position.
   118-127 F10.3  mm    formal uncertainty of Z-component of station position.
   139-145 I7     --    the number of observations of this station used in
                        solution.
   156-162 I7     --    total number of observations of this station.
   174-178 I5     --    the number of sessions with this station used in
                        solution.
   189-193 I5     --    total number of sessions of this station.
   205-214 A10    --    the date of the first session with this station used
                        in solution. format: yyyy.mm.dd (as integer numbers)
   226-235 A10    --    the date of the last session with this station used
                        in solution. format: yyyy.mm.dd (as integer numbers)

3) Local topocentric components of the vector of station position: Up, East,
   North. The first 8 characters of this line are STA_GCU:

   Field   Format Units Meaning
   1-8     A8     --    record type identifier: STA_GCU:
   11-25   A15    --    station name. Station name consist of 8-letters station
                        acronym and 6-letter epoch in format yymmdd. Epoch
                        is attached to the name only if the station had episodic
                        motion. Fields between the last letter of the station
                        name and the first letter of epoch are filled by _.
                        If the station didn't have episodic name then the space
                        after the last letter of the station name is left blank.
   28-29   A2     --    component identifier. One of "U:", "E:" or "N:"
   31-45   F15.2  mm    value of U-component of station position.
   50-59   F10.3  mm    formal uncertainty of U-component of station position.
   65-79   F15.2  mm    value of E-component of station position.
   84-93   F10.3  mm    formal uncertainty of E-component of station position.
   99-113  F15.2  mm    value of N-component of station position.
   118-127 F10.3  mm    formal uncertainty of N-component of station position.

4) Correlations between station positions and velocities. Correlation matrix
   is defined as the matrix of 6x6 in the upper triangle representation without
   the main diagonal which. Elements in the columns or rows of the matrix are
   in the order: X-position, Y-position, Z-position, X-velocity, Y-velocity,
   Z-velocity.

   1-8     A8     --    record type identifier: STA_CRL:
   11-25   A15    --    station name. Station name consist of 8-letters station
                        acronym and 6-letter epoch in format yymmdd. Epoch
                        is attached to the name only if the station had episodic
                        motion. Fields between the last letter of the station
                        name and the first letter of epoch are filled by _.
                        If the station didn't have episodic name then the space
                        after the last letter of the station name is left blank.
   31-36   F6.3   d/l   Correlation between X-position and Y-position
   38-43   F6.3   d/l   Correlation between X-position and Z-position
   45-50   F6.3   d/l   Correlation between Y-position and Z-position
   52-57   F6.3   d/l   Correlation between X-position and X-velocity
   59-64   F6.3   d/l   Correlation between Y-position and X-velocity
   66-71   F6.3   d/l   Correlation between Z-position and X-velocity
   73-78   F6.3   d/l   Correlation between X-position and Y-velocity
   80-85   F6.3   d/l   Correlation between Y-position and Y-velocity
   87-92   F6.3   d/l   Correlation between Z-position and Y-velocity
   94-99   F6.3   d/l   Correlation between X-velocity and Y-velocity
   101-106 F6.3   d/l   Correlation between X-position and Z-velocity
   108-113 F6.3   d/l   Correlation between Y-position and Z-velocity
   115-120 F6.3   d/l   Correlation between Z-position and Z-velocity
   122-127 F6.3   d/l   Correlation between X-velocity and Z-velocity
   129-134 F6.3   d/l   Correlation between Y-velocity and Z-velocity

"""

import numpy as np
import sys
from time_conv import date2year


# ------------------------------  FUNCTIONS  ---------------------------
def read_sta(dataflea):
    '''Retrieve the result from .lso file.

    Parameters
    ----------
    datafile : string
        name of data file

    Returns
    ----------
    staname : string
        name of station
    X : array, float
        X component
    Y : array, float
        Y component
    Z : array, float
        Z component
    X_err : array, float
        formal uncertainty of X component
    Y_err : array, float
        formal uncertainty of Y component
    Z_err : array, float
        formal uncertainty of Z component
    U : array, float
        U component
    E : array, float
        E component
    N : array, float
        N component
    U_err : array, float
        formal uncertainty of U component
    E_err : array, float
        formal uncertainty of E component
    N_err : array, float
        formal uncertainty of N component
    ObsUse : array,
        Number of used observations of this source
    ObsTot : array, int
        Total number of observations of this source
    SesUse : array, int
        Number of used sessions for this source
    SesTot : array, int
        Total number of sessions for this source
    DateBeg : array, float
        Epoch of the first observation
    DateEnd : array, float
        Epoch of the last observation
    '''

    # empty list for store data
    staname = []
    # XYZ
    X = []
    X_err = []
    Y = []
    Y_err = []
    Z = []
    Z_err = []
    # UEN
    U = []
    U_err = []
    E = []
    E_err = []
    N = []
    N_err = []
    #
    ObsUse = []
    ObsTot = []
    SesUse = []
    SesTot = []
    DateBeg = []
    DateEnd = []
    # Correlation
    XpYp = []
    XpZp = []
    YpZp = []
    XpXv = []
    YpXv = []
    ZpXv = []
    XpYv = []
    YpYv = []
    ZpYv = []
    XvYv = []
    XpZv = []
    YpZv = []
    ZpZv = []
    XvZv = []
    YvZv = []

    for line in open(datafile, 'r'):
        if line[0] != '#':
            if line[:7] == 'STA_GCX':
                staname.append(line[10:25].rstrip())
                X.append(float(line[30:45]))
                X_err.append(float(line[49:59]))
                Y.append(float(line[64:79]))
                Y_err.append(float(line[83:93]))
                Z.append(float(line[98:113]))
                Z_err.append(float(line[117:127]))
                ObsUse.append(int(line[138:145]))
                ObsTot.append(int(line[155:162]))
                SesUse.append(int(line[173:178]))
                SesTot.append(int(line[188:193]))
                DateBeg.append(date2year(line[204:214]))
                DateEnd.append(date2year(line[225:235]))
            elif line[:7] == 'STA_GCU':
                U.append(float(line[30:45]))
                U_err.append(float(line[49:59]))
                E.append(float(line[64:79]))
                E_err.append(float(line[83:93]))
                N.append(float(line[98:113]))
                N_err.append(float(line[117:127]))
            elif line[:7] == 'STA_CRL':
                XpYp.append(float(line[30:36]))
                XpZp.append(float(line[37:43]))
                YpZp.append(float(line[44:50]))
                XpXv.append(float(line[51:57]))
                YpXv.append(float(line[59:64]))
                ZpXv.append(float(line[65:71]))
                XpYv.append(float(line[72:78]))
                YpYv.append(float(line[79:85]))
                ZpYv.append(float(line[86:92]))
                XvYv.append(float(line[93:99]))
                XpZv.append(float(line[100:106]))
                YpZv.append(float(line[107:113]))
                ZpZv.append(float(line[114:120]))
                XvZv.append(float(line[121:127]))
                YvZv.append(float(line[128:134]))
            else:
                print("Something must be wrong! Please check your file")
                exit()

    # List -> array, a rather stupid way.
    staname = np.array(staname)
    # XYZ
    X = np.array(X)
    X_err = np.array(X_err)
    Y = np.array(Y)
    Y_err = np.array(Y_err)
    Z = np.array(Z)
    Z_err = np.array(Z_err)
    # UEN
    U = np.array(U)
    U_err = np.array(U_err)
    E = np.array(E)
    E_err = np.array(E_err)
    N = np.array(N)
    N_err = np.array(N_err)
    #
    ObsUse = np.array(ObsUse)
    ObsTot = np.array(ObsTot)
    SesUse = np.array(SesUse)
    SesTot = np.array(SesTot)
    DateBeg = np.array(DateBeg)
    DateEnd = np.array(DateEnd)
    # Correlation
    XpYp = np.array(XpYp)
    XpZp = np.array(XpZp)
    YpZp = np.array(YpZp)
    XpXv = np.array(XpXv)
    YpXv = np.array(YpXv)
    ZpXv = np.array(ZpXv)
    XpYv = np.array(XpYv)
    YpYv = np.array(YpYv)
    ZpYv = np.array(ZpYv)
    XvYv = np.array(XvYv)
    XpZv = np.array(XpZv)
    YpZv = np.array(YpZv)
    ZpZv = np.array(ZpZv)
    XvZv = np.array(XvZv)
    YvZv = np.array(YvZv)

    return [staname, X, X_err, Y, Y_err, Z, Z_err,
            U, U_err, E, E_err, N, N_err,
            ObsUse, ObsTot, SesUse, SesTot, DateBeg, DateEnd,
            XpYp, XpZp, YpZp, XpXv, YpXv, ZpXv,
            XpYv, YpYv, ZpYv, XvYv, XpZv, YpZv, ZpZv, XvZv, YvZv]


# Retrieve estimates.
if len(sys.argv) == 1:
    datafile = 'result/test.sta'
else:
    datafile = sys.argv[1]
[staname, X, X_err, Y, Y_err, Z, Z_err,
 U, U_err, E, E_err, N, N_err,
 ObsUse, ObsTot, SesUse, SesTot, DateBeg, DateEnd,
 XpYp, XpZp, YpZp, XpXv, YpXv, ZpXv,
 XpYv, YpYv, ZpYv, XvYv, XpZv, YpZv, ZpZv, XvZv, YvZv] = read_sta(datafile)
print(staname[0],
      X[0],
      X_err[0],
      Y[0],
      Y_err[0],
      Z[0],
      Z_err[0],
      U[0],
      U_err[0],
      E[0],
      E_err[0],
      N[0],
      N_err[0],
      ObsUse[0],
      ObsTot[0],
      SesUse[0],
      SesTot[0],
      DateBeg[0],
      DateEnd[0],
      XpYp[0],
      XpZp[0],
      YpZp[0],
      XpXv[0],
      YpXv[0],
      ZpXv[0],
      XpYv[0],
      YpYv[0],
      ZpYv[0],
      XvYv[0],
      XpZv[0],
      YpZv[0],
      ZpZv[0],
      XvZv[0],
      YvZv[0])
# ------------------------------ END -----------------------------------