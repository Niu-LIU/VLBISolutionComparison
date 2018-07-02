#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_GaiaDR2.py
"""
Created on Thu Apr 26 17:32:24 2018

@author: Neo(liuniu@smail.nju.edu.cn)

"""

import numpy as np
from astropy.io import fits


# -----------------------------  FUNCTIONS -----------------------------
def read_gaiadr2_fits(datafile):
    '''Read Gaia DR2 data file.

    Parameters
    ----------
    datafile : string
        Gaia DR2 file with the full path. The format is FITS.

    Returns
    ----------

    '''

    hdulist = fits.open(datafile, memmap=True)
    tbdat = hdulist[1].data

    source_id = tbdat.field("source_id")
    ref_epoch = tbdat.field("ref_epoch")
    ra = tbdat.field("ra")
    ra_error = tbdat.field("ra_error")
    dec = tbdat.field("dec")
    dec_error = tbdat.field("dec_error")
    parallax = tbdat.field("parallax")
    parallax_error = tbdat.field("parallax_error")
    pmra = tbdat.field("pmra")
    pmra_error = tbdat.field("pmra_error")
    pmdec = tbdat.field("pmdec")
    pmdec_error = tbdat.field("pmdec_error")
    ra_dec_corr = tbdat.field("ra_dec_corr")
    ra_parallax_corr = tbdat.field("ra_parallax_corr")
    ra_pmra_corr = tbdat.field("ra_pmra_corr")
    ra_pmdec_corr = tbdat.field("ra_pmdec_corr")
    dec_parallax_corr = tbdat.field("dec_parallax_corr")
    dec_pmra_corr = tbdat.field("dec_pmra_corr")
    dec_pmdec_corr = tbdat.field("dec_pmdec_corr")
    parallax_pmra_corr = tbdat.field("parallax_pmra_corr")
    parallax_pmdec_corr = tbdat.field("parallax_pmdec_corr")
    pmra_pmdec_corr = tbdat.field("pmra_pmdec_corr")
    frame_rotator_object_type = tbdat.field("frame_rotator_object_type")
    matched_observations = tbdat.field("matched_observations")
    phot_g_mean_mag = tbdat.field("phot_g_mean_mag")
    phot_bp_mean_mag = tbdat.field("phot_bp_mean_mag")
    phot_rp_mean_mag = tbdat.field("phot_rp_mean_mag")

    # Test
    # print(source_id)
    # print(source_id.size)
    # return frame_rotator_object_type

    return [source_id, ref_epoch,
            ra, ra_error, dec, dec_error, pmra, pmra_error, pmdec, pmdec_error,
            parallax, parallax_error,
            ra_dec_corr, ra_pmdec_corr, ra_pmra_corr, ra_pmdec_corr,
            dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr,
            parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr,
            phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
            frame_rotator_object_type, matched_observations]


def read_gaiadr2_iers_position(datafile):
    '''Read Gaia DR2 data file.

    Parameters
    ----------
    datafile : string
        Gaia DR2 file with the full path. The format is FITS.

    Returns
    ----------

    '''

    hdulist = fits.open(datafile, memmap=True)
    tbdat = hdulist[1].data

    # source_id = tbdat.field("source_id")
    # ref_epoch = tbdat.field("ref_epoch")
    iers_name = tbdat.field("iers_name")
    ra = tbdat.field("ra")
    ra_error = tbdat.field("ra_error")
    dec = tbdat.field("dec")
    dec_error = tbdat.field("dec_error")
    ra_dec_corr = tbdat.field("ra_dec_corr")

    # Test
    # print(source_id)
    # print(source_id.size)
    # return frame_rotator_object_type

    return iers_name, ra, ra_error, dec, dec_error, ra_dec_corr


def read_gaiadr2_ascii(datafile):
    '''Read Gaia DR2 data file.

    Parameters
    ----------
    datafile : string
        Gaia DR2 file with the full path. The format is ascii.

    Returns
    ----------
    (Later)
    '''

    RA, e_RA, DE, e_DE = np.genfromtxt(datafile, usecols=range(4),
                                       unpack=True)
    source = np.genfromtxt(datafile, usecols=(4,), dtype=str)

    return source, RA, e_RA, DE, e_DE


def test_code():
    '''For test
    '''

    # read_gaiadr2_fits("/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/"
    #                   "gaiadr2_iers.fits")
    # read_gaiadr2_fits("/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/"
    #                   "gaiadr2_agn.fits")

    [source_id, ref_epoch,
     ra, ra_error, dec, dec_error, pmra, pmra_error, pmdec, pmdec_error,
     parallax, parallax_error,
     ra_dec_corr, ra_pmdec_corr, ra_pmra_corr, ra_pmdec_corr,
     dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr,
     parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr,
     phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
     frame_rotator_object_type, matched_observations] = read_gaiadr2_fits(
        "/Users/Neo/Astronomy/Data/catalogs/Gaia_DR2/gaiadr2_qso_all.fits")
    # "gaiadr2_qso_all.fits") # For test on server

    # # For frame_rotator_object_type
    # print("Number of dataset: %10d" % frame_rotator_object_type.size)
    # # Count
    # for num in range(4):
    #     print("frame_rotator_object_type = %d: %10d" %
    #           (num, np.where(frame_rotator_object_type == num)[0].size))
    '''Output
    Number of dataset:     556869
    frame_rotator_object_type = 0:      81485
    frame_rotator_object_type = 1:          0
    frame_rotator_object_type = 2:       2820
    frame_rotator_object_type = 3:     472564
    '''

    # Table 2 in Mignard et al. 2018
    # 1) source selection
    # pm = np.sqrt(pmra**2 + pmdec**2)
    # G = phot_g_mean_mag
    # # condition
    # condition = [(pm <= 2) & (G < 18),
    #              (pm <= 2) & (G < 19),
    #              (pm <= 3) & (G < 20),
    #              (pm <= 3),
    #              ((ra * 1.e5).astype(int) % 2 == 0),
    #              ((ra * 1.e5).astype(int) % 2 == 0),
    #              (G > 19)]
    # label = ["u < 2mas/yr, G < 18",
    #          "u < 2mas/yr, G < 19",
    #          "u < 3mas/yr, G < 20",
    #          "u < 3mas/yr",
    #          "[10^5 * ra] mod 2 = 0",
    #          "[10^5 * ra] mod 2 = 1",
    #          "G > 19"]

    # print("Condition              |       No")
    # for cond, lab in zip(condition, label):
    #     source_id_n = source_id[cond]
    #     print("%-22s | %10d" % (lab, source_id_n.size))
    """Output
    Condition              |       No
    u < 2mas/yr, G < 18    |      27189
    u < 2mas/yr, G < 19    |     149146
    u < 3mas/yr, G < 20    |     400472
    u < 3mas/yr            |     513270
    [10^5 * ra] mod 2 = 0  |     278170
    [10^5 * ra] mod 2 = 1  |     278170  # --> 278169
    G > 19                 |     406356
    """

    from vsh_deg1_cor import VSHdeg01_fitting
    import datetime

    # max_num = np.arange(1000, 50000, 1000)
    # run_time = np.zeros_like(max_num)

    # Log file.
    flog = open('../logs/GaiaCRF2_test_vsh01.log', 'w')
    # flog = open('GaiaCRF2_test_vsh01.log', 'w')  # For test on server

    # degree -> rad
    ra_rad, dec_rad = np.deg2rad(ra), np.deg2rad(dec)
    # Full 6-parameters
    # print('VSH deg01:')
    w, sig, corrcoef, _, _, _ = VSHdeg01_fitting(
        pmra, pmdec, pmra_error, pmdec_error, ra_rad, dec_rad,
        flog, cov=pmra_pmdec_corr*pmra_error*pmdec_error,
        elim_flag="None")

    print("Estimations: ")
    print("Dipole:   ", w[:3])
    print("          ", sig[:3])
    print("Rotation: ", w[3:])
    print("          ", sig[3:])
    # print('sigma: ', sig)
    print("correlation: ", corrcoef)
    '''output
    Estimations:
    Dipole:    [-0.00702013  0.00465856  0.0121467 ]
               [0.00072152 0.00067442 0.00067208]
    Rotation:  [-0.00357653 -0.00224101 -0.00092421]
               [0.0007353  0.00064518 0.00080465]
    correlation:
    [[ 1.00000000e+00  5.76827353e-02 -8.88503550e-02 -1.47700531e-02
      -3.02289216e-02  4.51743055e-02]
     [ 5.76827353e-02  1.00000000e+00 -2.77233167e-02  1.74718608e-02
      -9.89511794e-04  6.65720258e-03]
     [-8.88503550e-02 -2.77233167e-02  1.00000000e+00 -5.02517072e-02
      -4.16227579e-02  2.23532319e-02]
     [-1.47700531e-02  1.74718608e-02 -5.02517072e-02  1.00000000e+00
       1.28544582e-01 -3.34799839e-01]
     [-3.02289216e-02 -9.89511794e-04 -4.16227579e-02  1.28544582e-01
       1.00000000e+00 -1.92617582e-01]
     [ 4.51743055e-02  6.65720258e-03  2.23532319e-02 -3.34799839e-01
      -1.92617582e-01  1.00000000e+00]]
    '''

    # Only rotation(spin)
    print('VSH deg01 (only rotation):')
    w, sig, corrcoef, _, _, _ = VSHdeg01_fitting(
        pmra, pmdec, pmra_error, pmdec_error, ra_rad, dec_rad,
        flog, cov=pmra_pmdec_corr*pmra_error*pmdec_error,
        elim_flag="None", fit_type="rotation")
    print("Estimations: ")
    print("Rotation: ", w)
    print("          ", sig)
    # print('sigma: ', sig)
    print("correlation: ", corrcoef)
    '''output
    Rotation:    [-0.00312404 -0.00193361 -0.00096809]
                 [0.00073412 0.00064424 0.00080354]
    correlation:
    [[ 1.          0.12618357 -0.33383656]
     [ 0.12618357  1.         -0.19061228]
     [-0.33383656 -0.19061228  1.        ]]
    '''

    flog.close()
    print("Done!")


# test_code()
# --------------------------------- END --------------------------------
