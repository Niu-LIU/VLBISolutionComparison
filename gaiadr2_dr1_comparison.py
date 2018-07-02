#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: gaiadr2_dr1_comparison.py
"""
Created on Thu May  3 15:39:10 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""


from astropy.io import fits
import numpy as np
from numpy import cos, sin, deg2rad
import time
from vsh_deg1_cor import VSHdeg01_fitting
from vector_direction import vec6_calc

# -----------------------------  FUNCTIONS -----------------------------


# -----------------------------  MAIN  ---------------------------------
# data file
datafile = ("/Users/Neo/Astronomy/Works/201711_GDR2_ICRF3/data/"
            "icrf2_gaiadr1_gaiadr2.fits")

# Log file.
flog = open('../logs/GaiaDR2_DR1_vsh.log', 'w')
print('## LOG FILE\n'
      '## Data: %s \n%s' %
      (datafile,
       time.strftime('##%Y-%m-%d %H:%M:%S Begins!',
                     time.localtime(time.time()))),
      file=flog)

# Load data.
hdulist = fits.open(datafile, memmap=True)
tbdat = hdulist[1].data

# Get the Gaia DR1 data
ra1 = tbdat.field("ra_gaiadr1")
ra_error1 = tbdat.field("ra_err_gaiadr1")
dec1 = tbdat.field("dec_gaiadr1")
dec_error1 = tbdat.field("dec_err_gaiadr1")
ra_dec_corr1 = tbdat.field("ra_dec_corr_gaiadr1")

# Get the Gaia DR2 data
ra2 = tbdat.field("ra")
ra_error2 = tbdat.field("ra_error")
dec2 = tbdat.field("dec")
dec_error2 = tbdat.field("dec_error")
parallax2 = tbdat.field("parallax")
parallax_error2 = tbdat.field("parallax_error")
pmra2 = tbdat.field("pmra")
pmra_error2 = tbdat.field("pmra_error")
pmdec2 = tbdat.field("pmdec")
pmdec_error2 = tbdat.field("pmdec_error")
ra_dec_corr2 = tbdat.field("ra_dec_corr")
ra_parallax_corr2 = tbdat.field("ra_parallax_corr")
ra_pmra_corr2 = tbdat.field("ra_pmra_corr")
ra_pmdec_corr2 = tbdat.field("ra_pmdec_corr")
dec_parallax_corr2 = tbdat.field("dec_parallax_corr")
dec_pmra_corr2 = tbdat.field("dec_pmra_corr")
dec_pmdec_corr2 = tbdat.field("dec_pmdec_corr")
parallax_pmra_corr2 = tbdat.field("parallax_pmra_corr")
parallax_pmdec_corr2 = tbdat.field("parallax_pmdec_corr")
pmra_pmdec_corr2 = tbdat.field("pmra_pmdec_corr")
phot_g_mean_mag2 = tbdat.field("phot_g_mean_mag")

# deg2uas
deg2uas = 3.6e9

# degree -> rad
ra_rad, dec_rad = deg2rad(ra2), deg2rad(dec2)
pmra_pmdec_cov = pmra_pmdec_corr2 * pmra_error2 * pmdec_error2
# pmra_pmdec_cov = None

# Estimation the residual spin
print("######## Residual spin", file=flog)
# Full 6-parameters
print("# VSH deg01(Full 6-parameters):", file=flog)
x, sig, corrcoef, _, _, _ = VSHdeg01_fitting(
    pmra2, pmdec2, pmra_error2, pmdec_error2, ra_rad, dec_rad,
    flog, cov=pmra_pmdec_cov, elim_flag="None")
x, sig = x * 1.e3, sig * 1.e3
gx,  gy,  gz,  wx,  wy,  wz = x
egx, egy, egz, ewx, ewy, ewz = sig
(r1, alr1, der1, errr1, erralr1, errder1,
 g1, alg1, deg1, errg1, erralg1, errdeg1) = vec6_calc(x, sig)
print('## Rotation component:\n',
      ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
      '=> %4d +/- %3d' % (r1, errr1), file=flog)
print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
      (alr1, erralr1, der1, errder1), file=flog)
print('## Glide component:\n',
      ' %+4d +/- %3d |' * 3 % (gx, egx, gy, egy, gz, egz),
      '=> %4d +/- %3d' % (g1, errg1), file=flog)
print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
      (alg1, erralg1, deg1, errdeg1), file=flog)
print('##   correlation coefficients are:\n', corrcoef, file=flog)
print("## tex output:\n",
      "&"
      " $%+4d \pm$ %3d &" * 3 % (wx, ewx, wy, ewy, wz, ewz),
      " $%4d \pm$ %3d" % (r1, errr1),
      " $%+4d \pm$ %3d &" * 3 % (gx, egx, gy, egy, gz, egz),
      " $%4d \pm$ %3d &(%.0f $\pm$ %.0f, $%+.0f \pm$ %.0f) \\" %
      (g1, errg1, alg1, erralg1, deg1, errdeg1), file=flog)

# print("Estimations: ")
# print("Dipole:   ", x[:3])
# print("          ", sig[:3])
# print("Spin:     ", x[3:])
# print("          ", sig[3:])
# print("correlation: ", corrcoef)
'''Output
VSH deg01:
Estimations:
Dipole:    [-6.64758806 -0.77685935 23.9621666 ]
           [5.69271791 5.59678426 4.94109605]
Spin:      [-18.6882817  -11.47197688  -0.12280909]
           [5.5966733  5.23016419 5.96978367]
correlation:  [[ 1.00000000e+00  6.16927228e-02 -3.24150686e-02  3.69331487e-02
  -3.70445428e-01  8.76878728e-02]
 [ 6.16927228e-02  1.00000000e+00 -6.91805674e-02  3.60423553e-01
  -9.46049539e-04  1.10424337e-02]
 [-3.24150686e-02 -6.91805674e-02  1.00000000e+00 -1.20426238e-01
  -7.84785791e-02 -1.36535746e-02]
 [ 3.69331487e-02  3.60423553e-01 -1.20426238e-01  1.00000000e+00
   9.06400860e-02 -2.24654248e-01]
 [-3.70445428e-01 -9.46049539e-04 -7.84785791e-02  9.06400860e-02
   1.00000000e+00 -2.08846034e-01]
 [ 8.76878728e-02  1.10424337e-02 -1.36535746e-02 -2.24654248e-01
  -2.08846034e-01  1.00000000e+00]]
'''


# Only rotation(spin)
print("\n\n# VSH deg01(Only rotation):", file=flog)
x, sig, corrcoef, _, _, _ = VSHdeg01_fitting(
    pmra2, pmdec2, pmra_error2, pmdec_error2, ra_rad, dec_rad,
    flog, cov=pmra_pmdec_cov, elim_flag="None", fit_type="rotation")
x, sig = x * 1.e3, sig * 1.e3
wx,  wy,  wz = x
ewx, ewy, ewz = sig

r1 = np.sqrt(np.dot(x, x))
errr1 = np.sqrt(np.dot(sig, sig))

print('## Rotation component:\n',
      ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
      '=> %4d +/- %3d' % (r1, errr1), file=flog)
print('##   correlation coefficients are:\n', corrcoef, file=flog)
print("## tex output:\n",
      # "&"
      " $%+4d \pm$ %3d &" * 3 % (wx, ewx, wy, ewy, wz, ewz),
      " $%4d \pm$ %3d" % (r1, errr1),
      "  &-  &-  &-  &-  &- \\\\", file=flog)

# print("Estimations: ")
# print("Spin:     ", x)
# print("          ", sig)
# print("correlation: ", corrcoef)
'''Output
VSH deg01:
Estimations:
Spin:      [-15.74023763 -11.47629987   0.79176579]
           [5.19251769 4.83421476 5.94636473]
correlation:  [[ 1.          0.10254288 -0.24986854]
 [ 0.10254288  1.         -0.1927077 ]
 [-0.24986854 -0.1927077   1.        ]]
'''


# direct comparison
# Positional difference (DR2 - DR1)
# degree --> uas
dra = (ra2 - ra1) * cos(dec_rad) * deg2uas
ddec = (dec2 - dec1) * deg2uas
# mas --> uas
dra_error = np.sqrt(ra_error1**2 + ra_error2**2) * 1.e3
ddec_error = np.sqrt(dec_error1**2 + dec_error2**2) * 1.e3
dra_ddec_cov = (ra_dec_corr1 * ra_error1 * dec_error1 +
                ra_dec_corr2 * ra_error2 * dec_error2) * 1.e6

print("\n\n########## Gaia DR2 - DR1 position (direct comparison)", file=flog)
print("# VSH deg01:", file=flog)
x, sig, corrcoef, _, _, _ = VSHdeg01_fitting(
    dra, ddec, dra_error, ddec_error, ra_rad, dec_rad,
    flog, cov=dra_ddec_cov, elim_flag="angsep", ang_sep=10e3)

gx,  gy,  gz,  wx,  wy,  wz = x
egx, egy, egz, ewx, ewy, ewz = sig
(r1, alr1, der1, errr1, erralr1, errder1,
 g1, alg1, deg1, errg1, erralg1, errdeg1) = vec6_calc(x, sig)
print('## Rotation component:\n',
      ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
      '=> %4d +/- %3d' % (r1, errr1), file=flog)
print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
      (alr1, erralr1, der1, errder1), file=flog)
print('## Glide component:\n',
      ' %+4d +/- %3d |' * 3 % (gx, egx, gy, egy, gz, egz),
      '=> %4d +/- %3d' % (g1, errg1), file=flog)
print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
      (alg1, erralg1, deg1, errdeg1), file=flog)
print('##   correlation coefficients are:\n', corrcoef, file=flog)
print("## tex output:\n",
      # "&"
      " $%+4d \pm$ %3d &" * 3 % (wx, ewx, wy, ewy, wz, ewz),
      " $%4d \pm$ %3d" % (r1, errr1),
      " $%+4d \pm$ %3d &" * 3 % (gx, egx, gy, egy, gz, egz),
      " $%4d \pm$ %3d &(%.0f $\pm$ %.0f, $%+.0f \pm$ %.0f) \\" %
      (g1, errg1, alg1, erralg1, deg1, errdeg1), file=flog)

# Epoch difference
deop = -0.5  # Julian year

# Epoch propagation
# mas/yr --> uas/yr
dra_n = dra + deop * pmra2 * 1.e3
ddec_n = ddec + deop * pmdec2 * 1.e3
dra_error_n = np.sqrt(dra_error**2 + deop ** 2 * pmra_error2**2 * 1.e6)
ddec_error_n = np.sqrt(ddec_error**2 + deop ** 2 * pmdec_error2**2 * 1.e6)

print("\n\n########## Gaia DR2 - DR1 position (same epoch J2015.0)", file=flog)
print("# VSH deg01:", file=flog)
x, sig, corrcoef, _, _, _ = VSHdeg01_fitting(
    dra_n, ddec_n, dra_error_n, ddec_error_n, ra_rad, dec_rad,
    flog, cov=dra_ddec_cov, elim_flag="angsep", ang_sep=10e3)

gx,  gy,  gz,  wx,  wy,  wz = x
egx, egy, egz, ewx, ewy, ewz = sig
(r1, alr1, der1, errr1, erralr1, errder1,
 g1, alg1, deg1, errg1, erralg1, errdeg1) = vec6_calc(x, sig)
print('## Rotation component:\n',
      ' %+4d +/- %3d |' * 3 % (wx, ewx, wy, ewy, wz, ewz),
      '=> %4d +/- %3d' % (r1, errr1), file=flog)
print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
      (alr1, erralr1, der1, errder1), file=flog)
print('## Glide component:\n',
      ' %+4d +/- %3d |' * 3 % (gx, egx, gy, egy, gz, egz),
      '=> %4d +/- %3d' % (g1, errg1), file=flog)
print('##       apex: (%.1f +/- %.1f, %.1f +/- %.1f)' %
      (alg1, erralg1, deg1, errdeg1), file=flog)
print('##   correlation coefficients are:\n', corrcoef, file=flog)
print("## tex output:\n",
      # "&"
      " $%+4d \pm$ %3d &" * 3 % (wx, ewx, wy, ewy, wz, ewz),
      " $%4d \pm$ %3d" % (r1, errr1),
      " $%+4d \pm$ %3d &" * 3 % (gx, egx, gy, egy, gz, egz),
      " $%4d \pm$ %3d &(%.0f $\pm$ %.0f, $%+.0f \pm$ %.0f) \\" %
      (g1, errg1, alg1, erralg1, deg1, errdeg1), file=flog)

flog.close()
print("Done!")
# --------------------------------- END --------------------------------
