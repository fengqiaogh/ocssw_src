#!/usr/bin/env python3
import datetime
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import sys
import time
from typing import Tuple, Union
import multiprocessing as mp
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.optimize import check_grad
import netCDF4 as nc
PROGRAM = "GPIG"
VERSION = "1.0.0"

G1 = 0.0949  # g1 and g2 are values from Gordon et al., 1988
G2 = 0.0794
LNOT = 400  # reference lambda wavelength (nm)

def rrs_inversion_pigments(Rrs : np.ndarray, Rrs_unc : np.ndarray, wl : np.ndarray, temp : np.ndarray, sal : np.ndarray, correction_file : str = "Sullivan_pure_water_temp_sal.csv" , tol_err = 1e-7, start_chunk: int = 0, check_gradient: bool = False) -> Tuple[np.ndarray, np.ndarray, str, np.ndarray]:
    '''Inversion algorithm to estimate phytoplankton pigments from Rrs spectra

    Ali Chase, University of Maine, 2017 - 2019

    See the following publication for details on the method:

        Chase, A., E. Boss, I. Cetinic, and W. Slade. 2017. "Estimation of Phytoplankton
        Accessory Pigments from Hyperspectral Reflectance Spectra: Toward a
        Global Algorithm."
        Journal of Geophysical Research: Oceans, doi: 10.1002/2017JC012859.

    NOTE: This code was developed to estimate phytoplankton pigments from
    hyperspectral remote-sensing reflectance (Rrs) data measured in situ at
    the ocean surface. The Rrs data used in this algorithm development were
    quality controlled, corrected for Raman scattering, normalized to
    eliminate the angular effect of the sun position in the sky relative to
    nadir. Please see above reference for details of these steps.

    Parameters:
    -----------
    Rrs : Numpy array
        Remote-sensing reflectance measurements, defined as Lw/Ed
    Rrs_unc : Numpy array
        Uncertainy values in Rrs measurements (e.g. could be the standard deviation
        in Rrs for a given five-minute sample collection), must be on the same 
        wavlength grid as the Rrs data
    wl : Numoy array
        wavelengths associated with Rrs (and Rrs_unc) data
    tem : int or float
        water temperature at the time of radiometery data collection
    sal : int or float
        water salinity at the time of radiometery data collection
    correction_file : str
        Path to the correction file for pure water absorption based on temperature and salinity
    tol_err : float
        Tolerance for the least squares inversion convergence
    start_chunk : int
        Starting pixel index for processing (used for parallel processing)
    check_gradient : bool
        If True, perform gradient checking for the Jacobian function.
        Default is False.

    Returns:
    --------
    numpy.ndarray
        Estimated pigment concentrations
    numpy.ndarray
        Uncertainties in estimated pigments,
        calculated using a Monte Carlo method that in turn uses the 
        reported uncertainties in the A and B coefficients reported
        in Chase et al. (2017). 
    str
        The names and units of the estimated pigments:
        chlorophyll a (Chla), chlorophyll b (Chlb), chlorophyll c1
        +c2 (Chlc12), and photoprotective carotenoids (PPC) defined
        as abcarotene+zeaxanthin+alloxanthin+diadinoxanthin. All
        pigments and uncertainties are in mg m^-3.
    numpy.ndarray
        Amplitudes of Gaussian absorption functions
        representing Chla, Chlb, Chlc12, and PPC. Can be used to derive
        updated relationships between Gaussians and HPLC pigments.
    '''

    # Get the absorption and backscattering by water for the temperature and
    # salinity measured coincidently with Rrs
    a_sw, bb_sw = get_water_iops(wl, correction_file, temp, sal)

    # Calculate rrs from Rrs (method from Lee et al., 2002)
    rrs = Rrs / (0.52 + 1.7 * Rrs)
    rrs_unc = Rrs_unc / (0.52 + 1.7 * Rrs_unc)

    # Calculate the positive solution for U using rrs, where U = bb/(bb+a).
    # This model and g coefficients are from Gordon et al., 1988
    Upos = (-G1 + np.sqrt(G1 ** 2 + 4 * G2 * rrs)) / (2 * G2)
    Uunc = (-G1 + np.sqrt(G1 ** 2 + 4 * G2 * rrs_unc)) / (2 * G2)

    # Define the center peak locations (nm) and widths (nm) of the Gaussian functions
    # sig = sigma, where FWHM = sigma*2.355 and FWHM = full width at half max
    peaks = np.array([384, 413, 435, 461, 464, 490, 532, 583, 623, 644, 655, 676], dtype=float)
    sig = np.array([23, 9, 14, 11, 19, 19, 20, 20, 15, 12, 12, 9], dtype=float)

    #   Define the [lower bound, first guess, upper bound] for each parameter. These will be allowed to vary.
    s_nap = [.005, .011, .016]
    m_nap = [.0, .005, .05]
    s_cdom = [.005, .0185, .02]
    m_cdom = [.01, .1, .8]
    bbpbp_ratio = [.005, .01, .015]
    m_gaus = [.0, .01, 0.5]
    cpgam = [.0, 1.0, 1.3]
    m_cp = [.01, .1, 1.0]

    C_fluor = [.0, .001, .01]

    # first guess array
    Amp0 = [s_nap[1], m_nap[1], s_cdom[1], m_cdom[1], bbpbp_ratio[1], cpgam[1], m_cp[1], C_fluor[1]] + np.tile(
        m_gaus[1], len(peaks)).tolist() + peaks.tolist() + sig.tolist()

    # lower bound array
    LB = [s_nap[0], m_nap[0], s_cdom[0], m_cdom[0], bbpbp_ratio[0], cpgam[0], m_cp[0], C_fluor[0]] + np.tile(m_gaus[0],
                                                                                                             len(peaks)).tolist() + (
                     peaks - 1).tolist() + (sig - 1).tolist()

    # upper bound array
    UB = [s_nap[2], m_nap[2], s_cdom[2], m_cdom[2], bbpbp_ratio[2], cpgam[2], m_cp[2], C_fluor[2]] + np.tile(m_gaus[2],
                                                                                                             len(peaks)).tolist() + (
                     peaks + 1).tolist() + (sig + 1).tolist()
    pixels, _ = Upos.shape
    Amp2_all = np.zeros((pixels,44))
    if check_gradient:
        print("Checking gradient accuracy...")
        Upos_single = Upos[0, :]
        Uunc_single = Uunc[0, :]
        wl_single = wl.squeeze()
        bb_sw_single = bb_sw[0, :]
        a_sw_single = a_sw[0, :]
        def cost_fn(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT):
            """Scalar objective = 0.5 * sum(residuals^2)."""
            r = lsqnonlin_Amp_gen(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT)
            return 0.5 * np.sum(r * r)

        def cost_grad_fn(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT):
            """Gradient of cost: J^T r."""
            r = lsqnonlin_Amp_gen(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT)
            J = jacobian_Amp_gen(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT)
            return J.T @ r   # shape (44,)

        # --- Wrapper functions so check_grad sees only f(x) and grad(x) ---
        def f_wrapped(x):
            return cost_fn(x, Upos_single, Uunc_single, wl_single,
                        bb_sw_single, a_sw_single, LNOT)

        def g_wrapped(x):
            return cost_grad_fn(x, Upos_single, Uunc_single, wl_single,
                                bb_sw_single, a_sw_single, LNOT)
        error = check_grad(f_wrapped, g_wrapped,  np.array(Amp0),epsilon=1e-6)
        print(f"Gradient error: {error}")  # Should be < 1e-6

        J_analytical = jacobian_Amp_gen( np.array(Amp0), Upos_single, Uunc_single, wl_single, bb_sw_single, a_sw_single, LNOT)
        J_numerical = np.zeros_like(J_analytical)
        eps=1e-7
        for i in range(len(Amp0)):
            Amp_plus = Amp0.copy()
            Amp_minus = Amp0.copy()
            Amp_plus[i] += eps
            Amp_minus[i] -= eps
            
            r_plus = lsqnonlin_Amp_gen(np.array(Amp_plus), Upos_single, Uunc_single, wl_single, bb_sw_single, a_sw_single, LNOT)
            r_minus = lsqnonlin_Amp_gen(np.array(Amp_minus), Upos_single, Uunc_single, wl_single, bb_sw_single, a_sw_single, LNOT)
            
            J_numerical[:, i] = (r_plus - r_minus) / (2 * eps)
        
        print(np.allclose(J_analytical, J_numerical, rtol=1e-5, atol=1e-5))
 

        def verify_gradient_components(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT, eps=1e-7):
            """Check both residual Jacobian AND cost gradient separately."""
            
            # 1. Check residual Jacobian (you said this passes)
            J_analytical = jacobian_Amp_gen(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT)
            J_numerical = np.zeros_like(J_analytical)
            
            for i in range(len(x)):
                x_plus = x.copy()
                x_minus = x.copy()
                x_plus[i] += eps
                x_minus[i] -= eps
                
                r_plus = lsqnonlin_Amp_gen(x_plus, Upos, Uunc, wl, bb_sw, a_sw, LNOT)
                r_minus = lsqnonlin_Amp_gen(x_minus, Upos, Uunc, wl, bb_sw, a_sw, LNOT)
                
                J_numerical[:, i] = (r_plus - r_minus) / (2 * eps)
            
            jac_match = np.allclose(J_analytical, J_numerical, rtol=1e-5, atol=1e-5)
            print(f"Residual Jacobian matches: {jac_match}")
            print(f"Max Jacobian error: {np.max(np.abs(J_analytical - J_numerical))}")
            
            # 2. Check cost gradient directly
            r = lsqnonlin_Amp_gen(x, Upos, Uunc, wl, bb_sw, a_sw, LNOT)
            grad_analytical = J_analytical.T @ r
            
            grad_numerical = np.zeros(len(x))
            for i in range(len(x)):
                x_plus = x.copy()
                x_minus = x.copy()
                x_plus[i] += eps
                x_minus[i] -= eps
                
                cost_plus = 0.5 * np.sum(lsqnonlin_Amp_gen(x_plus, Upos, Uunc, wl, bb_sw, a_sw, LNOT)**2)
                cost_minus = 0.5 * np.sum(lsqnonlin_Amp_gen(x_minus, Upos, Uunc, wl, bb_sw, a_sw, LNOT)**2)
                
                grad_numerical[i] = (cost_plus - cost_minus) / (2 * eps)
            
            grad_match = np.allclose(grad_analytical, grad_numerical, rtol=1e-5, atol=1e-5)
            print(f"Cost gradient matches: {grad_match}")
            print(f"Max gradient error: {np.max(np.abs(grad_analytical - grad_numerical))}")
            
            # Show which parameters have largest errors
            grad_errors = np.abs(grad_analytical - grad_numerical)
            worst_params = np.argsort(grad_errors)[-5:]
            print(f"\nWorst 5 parameters (indices): {worst_params}")
            print(f"Their analytical gradients: {grad_analytical[worst_params]}")
            print(f"Their numerical gradients: {grad_numerical[worst_params]}")
            print(f"Their errors: {grad_errors[worst_params]}")
            
            return jac_match, grad_match

        verify_gradient_components(np.array(Amp0), Upos_single, Uunc_single, wl_single,
                            bb_sw_single, a_sw_single, LNOT)
    # Run the inversion using a non-linear least squares inversion function
    # time to run inversion for all pixels

    A_first_guess = np.array(Amp0)
    result = None
    for ip in range(pixels):
        try:
            result = least_squares(lsqnonlin_Amp_gen,A_first_guess, jac=jacobian_Amp_gen,  bounds=(LB, UB), ftol=tol_err, xtol=tol_err, gtol=tol_err,
                               args=(Upos[ip,:], Uunc[ip,:], wl.squeeze(), bb_sw[ip,:], a_sw[ip,:], LNOT))
            Amp2 = result.x
        except:
            if result:
                Amp2 = result.x
            else:
                Amp2 = A_first_guess
        # A_first_guess = result.x
        Amp2_all[ip,:] = Amp2
        if ip == 0:
            start_time = time.time()
        if (ip+1) % 5000 == 0:
            end_time = time.time()
            print(f"#Proc ID {start_chunk//pixels}:Inversion completed in {end_time - start_time:.2f} seconds from {start_chunk} to {start_chunk + ip + 1} pixels")

    # Estimate pigment concentrations and their uncertainties with coefficients reported in
    # Chase et al., 2017 (JGR-Oceans) and a Monte Carlo method. 
    # NOTE: these coefficients have been recomputed using the same method in Chase et al., 2017 (JGR-Oceans) using new data
    # matrix:   A   A_unc   B   B_unc

    coeffs = np.array([
        [0.022, 0.008, 0.563, 0.068],
        [0.014, 0.009, 0.14, 0.059],
        [0.026, 0.013, 0.286, 0.074],
        [0.026, 0.024, 0.194, 0.105]])

    np.seterr(divide='ignore')

    pigmedian = np.zeros((Amp2_all.shape[0],4))
    pigunc = np.zeros((Amp2_all.shape[0],4))
    chunk_size = 5000
    chunks = [range(i, min(i + chunk_size, pixels)) for i in range(0, pixels, chunk_size)]
    for chunk in chunks:
        start_time = time.time()
        for ii in range(4):
            mc = np.random.randn(len(chunk), 10000, 1) * coeffs[ii, [1, 3]]
            As = coeffs[ii, 0] + mc[:, :, 0]
            As[As < 0] = 0  # prevent imaginary pigment values
            Bs = coeffs[ii, 2] + mc[:, :, 1]
            pigest = (Amp2_all[chunk, 10 + ii][:, np.newaxis]  / As) ** (1.0 / Bs)
            pigmedian[chunk, ii] = np.median(pigest, axis=1)
            prc = np.percentile(pigest, [16, 84], axis=1)
            pigunc[chunk, ii] = (prc[1,:] - prc[0,:]) / 2
        end_time = time.time()
        print(f"#Proc ID {start_chunk//pixels}: Monte Carlo Integration completed for chunk {chunk.start + start_chunk} to {chunk.stop + start_chunk} processed in {end_time - start_time:.2f} seconds")

    vars_units = 'Chla, Chlc1+c2, Chlb, PPC; all in mg m^-3'
    amps = Amp2_all[:, 10:14]

    return pigmedian, pigunc, vars_units, amps

def jacobian_Amp_gen(Amp0, Upos, Uunc, wvns, bb_sw_r, a_sw_r, lnot):
    """
    Corrected Jacobian for Amp_gen.
    Assumes:
      - wvns is shape (N,) or (N,1)
      - Upos, Uunc, bb_sw_r, a_sw_r are shape (N,)
      - Amp0 length = 44
      - amps, peak_locs, sig come in as (12,) and will be reshaped to (1,12)
    """

    N = wvns.size
    Jacobian = np.zeros((N, 44), dtype=float)

    # Unpack parameters
    slope_nap, mag_nap = Amp0[0], Amp0[1]
    slope_cdom, mag_cdom = Amp0[2], Amp0[3]
    bbp_ratio = Amp0[4]
    slope_cp = Amp0[5]
    mag_cp = Amp0[6]
    fluor_amp = Amp0[7]

    amps = Amp0[8:20]    # (12,)
    peak_locs = Amp0[20:32] # (12,)
    sig = Amp0[32:44]       # (12,)

    # broadcasting helpers
    wvns_col = wvns[:, np.newaxis]  # (N,1)
    peak_locs_b = peak_locs[np.newaxis,:]         # (1,12)
    sig_b = sig[np.newaxis, :]                     # (1,12)
    amps_b = amps[np.newaxis, :]                     # (1,12)

    # CDOM and NAP absorption
    CDOM = mag_cdom * np.exp(-slope_cdom * (wvns - lnot))   # (N,)
    NAP = mag_nap * np.exp(-slope_nap * (wvns - lnot))      # (N,)

    # Gaussian absorption (APHI)
    g = np.exp(-0.5 * ((wvns_col - peak_locs_b) / sig_b) ** 2)  # (N,12)
    APHI = np.sum(g * amps_b, axis=1)                           # (N,)

    # Total particulate absorption
    AP = NAP + APHI                                             # (N,)

    # cp (total particle concentration)
    CP = mag_cp * (wvns / lnot) ** (-slope_cp)                  # (N,)

    # BBP (backscatter from particles)
    BBP = bbp_ratio * (CP - AP)                                 # (N,)

    # Fluorescence Gaussian
    fluor = np.exp(-0.5 * ((wvns - 685.0) / 10.6) ** 2)          # (N,)
    F = fluor_amp * fluor                                       # (N,)

    # Modeled U = numer / denom
    denom = APHI + NAP + CDOM + a_sw_r + BBP + bb_sw_r         # (N,)
    numer = BBP + bb_sw_r + F                                  # (N,)
    Unew = numer / denom                                       # (N,)

    # Weighted residual (vector)
    # spec_min = (Upos - Unew) / Uunc   # not needed here, but used by gradient assembly

    # Pre-broadcast denom, numer, Uunc to (N,1) for multiplication with (N,12)
    denom_col = denom[:, np.newaxis]   # (N,1)
    numer_col = numer[:, np.newaxis]   # (N,1)
    Uunc_col = Uunc[:, np.newaxis]     # (N,1)

    # --- Derivatives for first 8 params (0..7) ---
    # slope_nap derivative
    # NAP depends on slope_nap: dNAP = -(wvns - lnot) * NAP
    dNAP_dslope = - (wvns - lnot) * NAP    # (N,)
    # dBBP/dNAP = -bbp_ratio
    # ddenom/ds = dAPHI/ds + dNAP_dslope + dCDOM/ds + dBBP/ds
    # but here derive using chain rule: do direct dU/dp formula

    # 0: slope_nap
    dnum_0 = -bbp_ratio * dNAP_dslope          # dnumer/dslope_nap via BBP
    ddenom_0 = dNAP_dslope - bbp_ratio * dNAP_dslope  # ddenom = dAPHI(0) + dNAP + dBBP (APHI 0 for this param)
    # simplify: ddenom_0 = (1 - bbp_ratio) * dNAP_dslope
    ddenom_0 = (1 - bbp_ratio) * dNAP_dslope
    Jacobian[:, 0] =-(dnum_0 * denom - numer * ddenom_0) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_0, ddenom_0)

    # 1: mag_nap
    dNAP_dmag = NAP / mag_nap
    dnum_1 = -bbp_ratio * dNAP_dmag
    ddenom_1 = (1 - bbp_ratio) * dNAP_dmag
    Jacobian[:, 1] =-(dnum_1 * denom - numer * ddenom_1) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_1, ddenom_1)

    # 2: slope_cdom
    dCDOM_dslope = - (wvns - lnot) * CDOM
    dnum_2 = 0.0 * dCDOM_dslope   # numer does not include CDOM
    ddenom_2 = dCDOM_dslope
    Jacobian[:, 2] =-(dnum_2 * denom - numer * ddenom_2) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_2, ddenom_2)

    # 3: mag_cdom
    dCDOM_dmag = np.exp(-slope_cdom * (wvns - lnot))
    dnum_3 = 0.0 * dCDOM_dmag
    ddenom_3 = dCDOM_dmag
    Jacobian[:, 3] =-(dnum_3 * denom - numer * ddenom_3) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_3, ddenom_3)

    # 4: bbp_ratio
    # BBP = bbp_ratio * (CP - AP)
    dBBP_dbp = (CP - AP)
    dnum_4 = dBBP_dbp
    # denom has +BBP, so ddenom same
    ddenom_4 = dBBP_dbp
    Jacobian[:, 4] =-(dnum_4 * denom - numer * ddenom_4) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_4, ddenom_4)

    # 5: slope_cp
    # CP = mag_cp * (wvns / lnot) ** (-slope_cp)
    # dCP/dslope_cp = -CP * ln(wvns/lnot)
    # dBBP/dslope_cp = bbp_ratio * dCP/dslope_cp
    dCP_ds = - mag_cp * (wvns / lnot) ** (-slope_cp) * np.log(wvns / lnot)
    dBBP_ds = bbp_ratio * dCP_ds
    dnum_5 = dBBP_ds
    ddenom_5 = dBBP_ds
    Jacobian[:, 5] =-(dnum_5 * denom - numer * ddenom_5) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_5, ddenom_5)

    # 6: mag_cp
    dCP_dmag = (wvns / lnot) ** (-slope_cp)
    dBBP_dmag = bbp_ratio * dCP_dmag
    dnum_6 = dBBP_dmag
    ddenom_6 = dBBP_dmag
    Jacobian[:, 6] =-(dnum_6 * denom - numer * ddenom_6) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_6, ddenom_6)

    # 7: fluor_amp
    dnum_7 = fluor  # dnumer/dfluor_amp = fluor
    ddenom_7 = 0.0  # denom does not depend on fluor_amp
    Jacobian[:, 7] =-(dnum_7 * denom - numer * ddenom_7) / (denom * denom)/ Uunc # dr_from_dnum_ddenom(dnum_7, ddenom_7)

    # ---------------------------------------------------------------------
    # Now the Gaussian parameters (amps, peaks, sigs) — 12 columns each
    # For these parameters APHI changes; numer changes via BBP (through AP)
    # We compute the common multiplier:
    #
    # M = (bbp_ratio*denom + (1 - bbp_ratio)*numer) / (denom**2 * Uunc)
    #
    # and each jacobian column = dAPHI_param * M
    # ---------------------------------------------------------------------
    # Recompute denom_col, numer_col, Uunc_col in case they were mutated earlier:
    denom_col = denom[:, np.newaxis]   # (N,1)
    numer_col = numer[:, np.newaxis]   # (N,1)
    Uunc_col = Uunc[:, np.newaxis]     # (N,1)

    # M as derived (shape (N,1))
    M = (bbp_ratio * denom_col + (1.0 - bbp_ratio) * numer_col) / (denom_col * denom_col * Uunc_col)
    # M is positive in typical cases

    # g is (N,12), amps_b is (1,12) -> broadcast to (N,12)
    # dAPHI / d(amps_j) = g[:, j]
    dAPHI_amps = g * 1.0                        # (N,12)

    # dAPHI / d(peak_j) = amps_j * g * (wvns - peak_j) / sig_j^2
    dAPHI_peaks = g * amps_b * (wvns_col - peak_locs_b) / (sig_b * sig_b)   # (N,12)

    # dAPHI / d(sig_j) = amps_j * g * (wvns - peak_j)^2 / sig_j^3
    dAPHI_sigs = g * amps_b * ((wvns_col - peak_locs_b) ** 2) / (sig_b ** 3)  # (N,12)

    # Multiply by M (N,1) -> broadcasts to (N,12)
    Jacobian[:, 8:20]  = dAPHI_amps  * M
    Jacobian[:, 20:32] = dAPHI_peaks * M
    Jacobian[:, 32:44] = dAPHI_sigs  * M
    return Jacobian

def lsqnonlin_Amp_gen(Amp0, Upos, Uunc, wvns, bb_sw_r, a_sw_r, lnot):
    '''
    The following function uses a non-linear least squares solver to minimize
    the difference between the measured Rrs and the modeled Rrs
    '''

    # Unpack parameters
    slope_nap, mag_nap = Amp0[0], Amp0[1]
    slope_cdom, mag_cdom = Amp0[2], Amp0[3]
    bbp_ratio = Amp0[4]
    slope_cp = Amp0[5]
    mag_cp = Amp0[6]
    fluor_amp = Amp0[7]

    amps = Amp0[8:20]  # Amplitudes of Gaussians
    peak_locs = Amp0[20:32]  # Centers
    sig = Amp0[32:44]  # Widths

    # Ensure broadcasting works
    wvns_col = wvns[:, np.newaxis]  # Shape (N, 1)
    peak_locs = peak_locs[np.newaxis, :]  # Shape (1, 12)
    sig = sig[np.newaxis, :]
    amps = amps[np.newaxis, :]

    # CDOM and NAP absorption
    CDOM = mag_cdom * np.exp(-slope_cdom * (wvns - lnot))
    NAP = mag_nap * np.exp(-slope_nap * (wvns - lnot))

    # Gaussian absorption (a_phi)
    gaus = np.exp(-0.5 * ((wvns_col - peak_locs) / sig) ** 2) * amps
    APHI = np.sum(gaus, axis=1)

    # Total particulate absorption
    AP = NAP + APHI

    # cp (total particle concentration)
    CP = mag_cp * (wvns / lnot) ** (-slope_cp)

    # BBP (backscatter from particles)
    BBP = bbp_ratio * (CP - AP)

    # Fluorescence Gaussian
    fluor = np.exp(-0.5 * ((wvns - 685) / 10.6) ** 2)
    F = fluor_amp * fluor

    # Modeled U = bb / (a + bb)
    denom = APHI + NAP + CDOM + a_sw_r + BBP + bb_sw_r
    numer = BBP + bb_sw_r + F
    Unew = numer / denom

    # Weighted residual
    spec_min = (Upos - Unew) / Uunc

    return spec_min


def RInw(
    lambda_: Union[int, float, np.ndarray],
    Tc: Union[int, float],
    S: Union[int, float],
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """Refractive index of air is from Ciddor (1996,Applied Optics)
    Refractive index of seawater is from Quan and Fry (1994, Applied Optics)
    """
    n_air = (
        1.0 + (5792105.0 / (238.0185 - 1 / (lambda_ / 1e3) ** 2)
        + 167917.0 / (57.362 - 1 / (lambda_ / 1e3) ** 2)) / 1e8
    )
    n0 = 1.31405
    n1 = 1.779e-4
    n2 = -1.05e-6
    n3 = 1.6e-8
    n4 = -2.02e-6
    n5 = 15.868
    n6 = 0.01155
    n7 = -0.00423
    n8 = -4382
    n9 = 1.1455e6
    nsw = (
        n0 + (n1 + n2 * Tc + n3 * Tc ** 2) * S + n4 * Tc ** 2
        + (n5 + n6 * S + n7 * Tc) / lambda_ + n8 / lambda_ ** 2 + n9
        / lambda_ ** 3
    )
    nsw = nsw * n_air
    dnswds = (n1 + n2 * Tc + n3 * Tc ** 2 + n6 / lambda_) * n_air

    return nsw, dnswds


def BetaT(Tc: Union[int, float], S: Union[int, float]) -> float:
    """Pure water secant bulk Millero (1980, Deep-sea Research).
    Isothermal compressibility from Kell sound measurement in pure water.
    Calculate seawater isothermal compressibility from the secant bulk.
    """
    kw = (
            19652.21 + 148.4206 * Tc - 2.327105 * Tc ** 2 + 1.360477e-2
            * Tc ** 3 - 5.155288e-5 * Tc ** 4
    )
    Btw_cal = 1 / kw
    a0 = 54.6746 - 0.603459 * Tc + 1.09987e-2 * Tc ** 2 - 6.167e-5 * Tc ** 3
    b0 = 7.944e-2 + 1.6483e-2 * Tc - 5.3009e-4 * Tc ** 2

    Ks = kw + a0 * S + b0 * S ** 1.5
    IsoComp = 1 / Ks * 1e-5

    return IsoComp


def rho_sw(Tc: Union[int, float], S: Union[int, float]) -> float:
    """Density of water and seawater, unit is Kg/m^3, from UNESCO,38,1981.
    TODO: Compare to GSW Oceanographic Toolbox code (or other updated) density equations.
    """
    a0 = 8.24493e-1
    a1 = -4.0899e-3
    a2 = 7.6438e-5
    a3 = -8.2467e-7
    a4 = 5.3875e-9
    a5 = -5.72466e-3
    a6 = 1.0227e-4
    a7 = -1.6546e-6
    a8 = 4.8314e-4
    b0 = 999.842594
    b1 = 6.793952e-2
    b2 = -9.09529e-3
    b3 = 1.001685e-4
    b4 = -1.120083e-6
    b5 = 6.536332e-9

    # density for pure water
    density_w = b0 + b1 * Tc + b2 * Tc ** 2 + b3 * Tc ** 3 + b4 * Tc ** 4 + b5 * Tc ** 5
    # density for pure seawater
    density_sw = (
            density_w + ((a0 + a1 * Tc + a2 * Tc ** 2 + a3 * Tc ** 3 + a4 * Tc ** 4) * S
                         + (a5 + a6 * Tc + a7 * Tc ** 2) * S ** 1.5 + a8 * S ** 2)
    )
    return density_sw


def dlnasw_ds(Tc: Union[int, float], S: Union[int, float]) -> float:
    """Water activity data of seawater is from Millero and Leung (1976,American
    Journal of Science,276,1035-1077). Table 19 was reproduced using
    Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
    dlnawds is partial derivative of natural logarithm of water activity
    w.r.t.salinity.

    lnaw =  (-1.64555e-6-1.34779e-7*Tc+1.85392e-9*Tc.^2-1.40702e-11*Tc.^3)+......
            (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3).*S+......
            (1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^1.5+......
            (-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S.^2;

    density derivative of refractive index from PMH model
    """

    dlnawds = (
            (-5.58651e-4 + 2.40452e-7 * Tc - 3.12165e-9 * Tc ** 2 + 2.40808e-11 * Tc ** 3)
            + 1.5 * (1.79613e-5 - 9.9422e-8 * Tc + 2.08919e-9 * Tc ** 2 - 1.39872e-11 * Tc ** 3) * S ** 0.5
            + 2 * (-2.31065e-6 - 1.37674e-9 * Tc - 1.93316e-11 * Tc ** 2) * S
    )
    return dlnawds


def PMH(n_wat: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    n_wat2 = n_wat ** 2
    n_density_derivative = (
            (n_wat2 - 1) * (1 + 2 / 3 * (n_wat2 + 2)
                            * (n_wat / 3 - 1 / 3 / n_wat) ** 2)
    )
    return n_density_derivative


def betasw124_ZHH2009(lambda_, S, Tc, delta=0.039):
    """Scattering by pure seawater: Effect of salinity
    Xiaodong Zhang, Lianbo Hu, and Ming-Xia He, Optics Express, 2009, accepted
    lambda (nm): wavelength
    Tc: temperauter in degree Celsius, must be a scalar
    S: salinity, must be scalar
    delta: depolarization ratio, if not provided, default = 0.039 will be used (from Farinato and Roswell (1976))
    betasw: volume scattering at angles defined by theta. Its size is [x y],
    where x is the number of angles (x = length(theta)) and y is the number
    of wavelengths in lambda (y = length(lambda))
    beta90sw: volume scattering at 90 degree. Its size is [1 y]
    bw: total scattering coefficient. Its size is [1 y]
    for backscattering coefficients, divide total scattering by 2
    Xiaodong Zhang, March 10, 2009

    MODIFIED on 17/05/2011 to be able to process bbp profiles with coincident T and sal profiles
    MODIFIED on 05 Apr 2013 to use 124 degs instead of 117 degs
    """

    # values of the constants
    Na = 6.0221417930e23  # Avogadro's constant
    Kbz = 1.3806503e-23  # Boltzmann constant
    Tk = Tc + 273.15  # Absolute temperature
    M0 = 18e-3  # Molecular weight of water in kg/mol

    theta = np.linspace(0.0, 180.0, 18_001)

    rad = theta * np.pi / 180  # angle in radians as a 1-d array

    # nsw: absolute refractive index of seawater
    # dnds: partial derivative of seawater refractive index w.r.t. salinity
    nsw, dnds = RInw(lambda_, Tc, S)

    # isothermal compressibility is from Lepple & Millero (1971,Deep Sea-Research), pages 10-11
    # The error ~ +/-0.004e-6 bar^-1
    IsoComp = BetaT(Tc, S)

    # density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
    density_sw = rho_sw(Tc, S)

    # water activity data of seawater is from Millero and Leung (1976,American
    # Journal of Science,276,1035-1077). Table 19 was reproduced using
    # Eq.(14,22,23,88,107) then were fitted to polynominal equation.
    # dlnawds is partial derivative of natural logarithm of water activity
    # w.r.t.salinity
    dlnawds = dlnasw_ds(Tc, S)

    # density derivative of refractive index from PMH model
    DFRI = PMH(nsw)  # PMH model

    # volume scattering at 90 degree due to the density fluctuation
    beta_df = (
            np.pi * np.pi / 2 * ((lambda_ * 1e-9) ** (-4)) * Kbz * Tk * IsoComp * DFRI ** 2
            * (6 + 6 * delta) / (6 - 7 * delta)
    )

    # volume scattering at 90 degree due to the concentration fluctuation
    flu_con = S * M0 * dnds ** 2 / density_sw / (-dlnawds) / Na
    beta_cf = (
            2 * np.pi * np.pi * ((lambda_ * 1e-9) ** (-4)) * nsw ** 2
            * (flu_con) * (6 + 6 * delta) / (6 - 7 * delta)
    )

    # total volume scattering at 90 degree
    beta90sw = beta_df + beta_cf
    bsw = 8 * np.pi / 3 * beta90sw * (2 + delta) / (1 + delta)

    for i, value in enumerate(rad):  # TODO: Is there a better way to do this?
        if np.rad2deg(value) >= 124:
            rad124 = i
            break

    betasw124 = (
            beta90sw * (1 + ((np.cos(rad[rad124])) ** 2) * (1 - delta) / (1 + delta))
    )

    return betasw124, bsw, beta90sw, theta


def tempsal_corr(lambda_, p):
    """
    Function to obtain temperature and salinity correction factors for pure water absorption
    from Sullivan et al 2006
    """
    pwts = pd.read_csv(p) # p - the correction file path 

    if lambda_.min() < 400 or lambda_.min() > 750:
        raise NotImplementedError(
            "Wavelengths < 400 or > 750 not currently supported because these "
            "values are outside the range for temperature salinitity correction "
            "provided in Sullivan et al 2006."
        )

    # next two lines changed from 'spline' on 14 Dec 2015
    psiT = np.interp(lambda_, pwts["lambda"], pwts["PsiT"])
    psiS = np.interp(lambda_, pwts["lambda"], pwts["PsiS"])

    return psiT, psiS


def get_water_iops(lambda_,correction_file, T=22, S=35):
    """Function to obtain seawater absorption and backscattering spectra.
    Pure water absorption from Mason et al 2016 for 250-550, Pope and Frye for
    550-730 nm, and Smith and Baker for 730-800 nm salt water backscattering from Zhang et al 2009
    corrected for in situ temperature and salinity conditions Sullivan et al. 2006

    Parameters
    ----------
    wavel: (array)
        Wavelenghts (nm) for iops to be returned at.
    T: (float or int)
        Temperature (degC) of the water.
    S: (float or int)
        Salinity (psu) of the water.

    Returns
    -------
    numpy.ndarray
        Absorption of seawater at each wavelength in wavel for the specified temperature and salinity values.
    numpy.ndarray
        Backscattering of seawater at each wavelenght in wavel for the specified temperature and salinity values.
    """
    lambda_water1 = [250, 260, 270, 280, 290] + list(range(300, 551, 2))
    lambda_water2 = list(np.linspace(552.5, 800.0, 100))
    lambda_water = lambda_water1 + lambda_water2

    aw1 = [
        58.7100, 51.5000, 43.5700, 22.3000, 9.3900, 4.6700, 4.3300, 3.5600, 3.1300, 2.7500, 2.3600, 2.0500,
        1.8500, 1.7600, 1.6300, 1.4700, 1.3300, 1.2400, 1.1800, 1.1200, 1.0700, 1.0100, 0.9900, 0.9500,
        0.9100, 0.8500, 0.8200, 0.8100, 0.8200, 0.8400, 0.8900, 0.9400, 0.9700, 0.9800, 0.9900, 1.0600,
        1.1500, 1.2000, 1.2100, 1.2200, 1.2400, 1.2700, 1.2900, 1.3300, 1.3700, 1.4300, 1.4700, 1.5100,
        1.5500, 1.6200, 1.7000, 1.7500, 1.8500, 1.9600, 2.0800, 2.2200, 2.3700, 2.4800, 2.5700, 2.5900,
        2.6600, 2.7100, 2.8000, 2.8800, 3.0000, 3.1200, 3.2200, 3.3100, 3.4400, 3.5800, 3.7600, 3.9500,
        4.1700, 4.4200, 4.8000, 5.2200, 5.7400, 6.2600, 6.9100, 7.5100, 8.0800, 8.4200, 8.6300, 8.7700,
        8.9300, 9.0900, 9.3300, 9.5500, 9.7900, 9.9900, 10.3000, 10.6500, 11.0000, 11.3800, 11.7700, 12.1400,
        12.5400, 12.9400, 13.3600, 13.9100, 14.6000, 15.4500, 16.4800, 17.7400, 19.2600, 20.7300, 22.4200, 24.2400,
        26.6800, 29.7100, 33.0000, 35.6900, 37.3800, 38.2100, 38.7800, 39.1700, 39.6200, 40.1700, 40.8800, 41.6200,
        42.4200, 43.3000, 44.3600, 45.4100, 46.4500, 47.5400, 48.8200, 50.4000, 52.2400, 54.2500, 56.2900,
    ]
    aw2 = [
        0.0593, 0.0596, 0.0606, 0.0619, 0.064, 0.0642, 0.0672, 0.0695, 0.0733, 0.0772, 0.0836, 0.0896, 0.0989, 0.11,
        0.122, 0.1351, 0.1516, 0.1672, 0.1925, 0.2224, 0.247, 0.2577, 0.2629, 0.2644, 0.2665, 0.2678, 0.2707, 0.2755,
        0.281, 0.2834, 0.2904, 0.2916, 0.2995, 0.3012, 0.3077, 0.3108, 0.322, 0.325, 0.335, 0.34, 0.358, 0.371, 0.393,
        0.41, 0.424, 0.429, 0.436, 0.439, 0.448, 0.448, 0.461, 0.465, 0.478, 0.486, 0.502, 0.516, 0.538, 0.559, 0.592,
        0.624, 0.663, 0.704, 0.756, 0.827, 0.914, 1.007, 1.119, 1.231, 1.356, 1.489, 1.678, 1.7845, 1.9333, 2.0822,
        2.2311,
        2.3800, 2.4025, 2.4250, 2.4475, 2.4700, 2.4900, 2.5100, 2.5300, 2.5500, 2.5400, 2.5300, 2.5200, 2.5100, 2.4725,
        2.4350, 2.3975, 2.3600, 2.3100, 2.2600, 2.2100, 2.1600, 2.1375, 2.1150, 2.0925, 2.0700,
    ]
    a_water = list(np.asarray(aw1) * 1e-3) + aw2  # 1e-3 comes from Mason et al. 2016, Table 2

    a_pw = np.interp(lambda_, lambda_water, a_water)

    # use salt water scattering from Zhang et al 2009
    _, b_sw, _, _ = betasw124_ZHH2009(lambda_, S, T)

    bb_sw = b_sw / 2  # divide total scattering by 2 to get the backscattering component

    # use the temp and salinity corrections from Sullivan et al. 2006
    psiT, psiS = tempsal_corr(lambda_,correction_file)

    # temperature and salinity corrections:
    T_norm = 22.0
    a_sw = (a_pw + psiT * (T - T_norm) + psiS * S)

    return a_sw, bb_sw

def cli_gpig():
        """Parse command-line arguments for the GPIG Code.

        Required positional arguments:
            - --ifile: path to an L2 file containing `Rrs` and `Rrs_unc` (netCDF)
            - --sst: path to a sea-surface temperature file (netCDF)
            - --sal: path to a salinity file (netCDF)
            - --ofile: path to an output file to write results to

        Optional:
            - --ncores: number of cores to use (default: 1)
            - --doi: DOI for the dataset (default: unspecified)
            - --processing_vers: Processing version for the dataset (default: unspecified)
            - --correction: Correction file for pure water absorption
                (default: Sullivan_pure_water_temp_sal.csv)
            - --tol: Tolerance level of the least square solution (default: 1e-
        The function only parses and returns the parsed args (argparse.Namespace).
        It does not change program behaviour or call the inversion; that can be
        performed by the caller using the returned values.
        """
        import argparse
        parser = argparse.ArgumentParser(description="GPIG Rrs inversion code")
        parser.add_argument('--ifile',"-i", help='Path to L2 netCDF file containing Rrs and Rrs_unc', required=True)
        parser.add_argument('--sst',"-t",default="$OCDATAROOT/common/sst_climatology.nc", help='Path to SST file (netCDF). Default is climatology')
        parser.add_argument('--sal',"-s",default="$OCDATAROOT/common/sss_climatology_woa2009.nc", help='Path to salinity file (netCDF). Default is climatology')
        parser.add_argument('--ofile',"-o", help='Path to output file (netCDF)',required=True)
        parser.add_argument('--ncores',"-n", type=int, default=1, help='Number of cores to use (default: 1)')
        parser.add_argument('--doi',"-d", type=str, default="unspecified", help='DOI for the dataset (default: unspecified)')
        parser.add_argument('--processing_version',"-p", type=str, default="unspecified", help='Processing version for the dataset (default: unspecified)')
        parser.add_argument('--correction',"-c", type=str, default="$OCDATAROOT/common/Sullivan_pure_water_temp_sal.csv", help='Correction file for pure water absorption (default: Sullivan_pure_water_temp_sal.csv)')
        parser.add_argument("--tol","-a",type=float,default=1e-7,help="Tolerance level of the least square solution. Default is 1e-7")
        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit(0)
        return parser.parse_args()


def read_salinity_file(path_sal, month) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    # read from SMAS salinity file or climatology HDF file
    with nc.Dataset(path_sal) as f:
        if "smap_sss"  in f.variables:
            sal_data = f['smap_sss'][:]
            sal_lon = f['longitude'][:]
            sal_lat = f['latitude'][:]
        elif "sss01" in f.variables:
            sal_data = f[f'sss{month:02d}'][:]
            sal_data[sal_data < 0] = np.nan  # set invalid values to NaN
            sal_lon = f['Longitude'][:]
            sal_lat = f['Latitude'][:]
        else:
            raise KeyError(f"Salinity variable not found in the provided file {path_sal}. Expected 'smap_sss' or 'sss01'.")


    return sal_data, sal_lon, sal_lat

def read_sst_file(path_sst,month) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # read from SST NetCDF file or climatology HDF file
    with nc.Dataset(path_sst) as f:
        if "analysed_sst"  in f.variables:
            sst_data = f['analysed_sst'][:].squeeze()
            sst_lon = f['lon'][:]
            sst_lat = f['lat'][:]
        elif "data01" in f.variables:
            sst_data = f[f'data{month:02d}'][:]
            slope = f[f'data{month:02d}'].slope
            intercept = f[f'data{month:02d}'].intercept
            sst_data = sst_data * slope + intercept + 273.15  # convert to Kelvin
            sst_data[sst_data < 0] = np.nan  # set invalid values to NaN
            sst_lon = f['Longitude'][:]
            sst_lat = f['Latitude'][:]
        else:
            raise KeyError(f"SST variable not found in the provided file {path_sst}. Expected 'analysed_sst' or 'data01'.")
    return sst_data, sst_lon, sst_lat


def interpolate_sst_salinity(path_sst, path_sal, path_rrs):
    """
    Interpolate SST and salinity from regular grids to scatter points.
    NaN values are replaced with nearest non-NaN neighbor.
    """
    # Load target points
    with nc.Dataset(path_rrs) as f:
        target_lon = f['navigation_data/longitude'][:]
        target_lat = f['navigation_data/latitude'][:]
        time_coverage_start = f.getncattr('time_coverage_start')
        month = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%S.%fZ").month
    # Load salinity data
    sal_data, sal_lon, sal_lat = read_salinity_file(path_sal, month)    
    # Load SST data
    sst_data, sst_lon, sst_lat = read_sst_file(path_sst, month)  
    
    # Interpolate salinity
    sal_interp = interpolate_with_nan_handling(
        sal_data, sal_lon, sal_lat, target_lon, target_lat
    )
    
    # Interpolate SST
    sst_interp = interpolate_with_nan_handling(
        sst_data, sst_lon, sst_lat, target_lon, target_lat
    )
    
    return sal_interp, sst_interp
def interpolate_with_nan_handling(grid_data, grid_lon, grid_lat, target_lon, target_lat):
    """
    Interpolate grid data to target points, replacing NaNs with nearest non-NaN value.
    
    Parameters:
    -----------
    grid_data : array-like (nlat, nlon) or (time, nlat, nlon)
        Grid data to interpolate
    grid_lon : array-like (nlon,)
        Longitude coordinates of grid
    grid_lat : array-like (nlat,)
        Latitude coordinates of grid
    target_lon : array-like
        Target longitude points
    target_lat : array-like
        Target latitude points
    
    Returns:
    --------
    result : numpy array
        Interpolated values at target points
    """
    # Handle multidimensional data (squeeze time dimension if present)
    original_shape = grid_data.shape
    if grid_data.ndim > 2:
        grid_data = np.squeeze(grid_data)
    
    # Convert masked arrays to regular arrays with NaN
    if hasattr(grid_data, 'filled'):
        grid_data = grid_data.filled(np.nan)
    
    # Flatten target points
    target_lon_flat = np.asarray(target_lon).flatten()
    target_lat_flat = np.asarray(target_lat).flatten()
    target_shape = np.asarray(target_lon).shape
    
    # Create interpolator for non-NaN regions
    result = np.full(len(target_lon_flat), np.nan)
    
    # First pass: regular interpolation
    if not np.all(np.isnan(grid_data)):
        # Create a version with NaNs filled temporarily for interpolation
        grid_filled = grid_data.copy()
        grid_filled[np.isnan(grid_filled)] = np.nanmean(grid_filled)
        
        interp = RegularGridInterpolator(
            (grid_lat, grid_lon), 
            grid_filled,
            method='nearest',
            bounds_error=False,
            fill_value=np.nan
        )
        
        points = np.column_stack([target_lat_flat, target_lon_flat])
        result = interp(points)
    
    # Second pass: for NaN results, find nearest non-NaN grid point
    nan_mask = np.isnan(result)
    if np.any(nan_mask) and not np.all(np.isnan(grid_data)):
        # Get coordinates of non-NaN grid points
        non_nan_mask = ~np.isnan(grid_data)
        lat_grid, lon_grid = np.meshgrid(grid_lat, grid_lon, indexing='ij')
        
        non_nan_coords = np.column_stack([
            lat_grid[non_nan_mask],
            lon_grid[non_nan_mask]
        ])
        non_nan_values = grid_data[non_nan_mask]
        
        if len(non_nan_coords) > 0:
            # Build KD-tree for fast nearest neighbor search
            tree = cKDTree(non_nan_coords)
            
            # Find nearest non-NaN point for each NaN result
            nan_points = np.column_stack([
                target_lat_flat[nan_mask],
                target_lon_flat[nan_mask]
            ])
            
            _, indices = tree.query(nan_points)
            result[nan_mask] = non_nan_values[indices]
    
    # Reshape to original target shape
    result = result.reshape(target_shape)
    
    return result
if __name__ == "__main__":
    print(PROGRAM,VERSION)
    args = cli_gpig()
    path_nc = args.ifile
    path_out = args.ofile
    path_salinity = os.path.expandvars(args.sal)
    path_sst = os.path.expandvars(args.sst)
    Sullivan_pure_water_temp_sal = os.path.expandvars(args.correction)  
    ncores = args.ncores
    tol_err = args.tol
    # check if all input files exist
    if not os.path.isfile(path_nc) or not os.path.isfile(path_sst) or not os.path.isfile(path_salinity) or not os.path.isfile(Sullivan_pure_water_temp_sal):
        raise FileNotFoundError("One or more input files do not exist.")
    sal_interp, sst_interp = interpolate_sst_salinity(path_sst, path_salinity, path_nc)
    # fill nan values in sst and salinity:
    sal_interp[np.isnan(sal_interp)] = 35.0
    sst_interp[np.isnan(sst_interp)] = 293.15 # 20 degC in Kelvin
     # open input file and output file
    if os.path.isfile(path_out):
        os.remove(path_out)
    try:
        fout = nc.Dataset(path_out, "w")
    except Exception as e:
        print(f"-E-: Could not create output file {path_out}. Error: {e}")
        sys.exit(1)
    with nc.Dataset(path_nc) as f:
        geophysical_data  : nc.Group = f["geophysical_data"]
        if "wavelength" not in geophysical_data.variables:
            wavelength = f["sensor_band_parameters/wavelength_3d"][:]
        else:
            wavelength = geophysical_data["wavelength"][:]
        # select wavelengths between 400 and 600 nm        
        mask = (wavelength > 400) & (wavelength < 600)
        wavelength = wavelength.data[mask]
        Rrs  : np.ndarray = geophysical_data["Rrs"][:,:,mask]
        Rrs_unc : np.ndarray  = geophysical_data["Rrs_unc"][:,:,mask]
        # select only positve Rrs values
        common_valid_mask = (~Rrs.mask).all(axis=2) & (Rrs.data > 0).all(axis=2)
        common_valid_mask = common_valid_mask & (~np.isnan(sst_interp)) & (~np.isnan(sal_interp))
        sstref = sst_interp - 273.15 # Convert from Kelvin to Celsius
        salref = sal_interp
        Rrs = Rrs.data[common_valid_mask, :]
        Rrs_unc =  Rrs_unc.data[common_valid_mask, :]
        if Rrs.shape[0] == 0:
            print("-E-: No valid pixels found for processing. Exiting.")
            fout.close()
            os.remove(path_out)
            sys.exit(110)
        sstref = sstref[common_valid_mask]
        salref = salref[common_valid_mask]
        start_time_all = time.time()
        multiprocessing_enabled = ncores > 1
        print("Total number of pixels to process:", Rrs.shape[0])
        if multiprocessing_enabled:
            # implement load balancing for multiprocessing
            perm = np.random.permutation(Rrs.shape[0])
            # randomly permute the input data to even out load across processes
            Rrs = Rrs[perm,:]
            Rrs_unc = Rrs_unc[perm,:]
            sstref = sstref[perm]
            salref = salref[perm]
            chunk_size = Rrs.shape[0] // ncores + 1
            chunks = [range(i, min(i + chunk_size, Rrs.shape[0])) for i in range(0, Rrs.shape[0], chunk_size)]
            args_list = [(Rrs[chunk, :], Rrs_unc[chunk, :], wavelength[None, :], sstref[chunk, None], salref[chunk, None],Sullivan_pure_water_temp_sal,tol_err, chunk[0]) for chunk in chunks]
            with mp.Pool(processes=ncores) as pool:
                results = pool.starmap(rrs_inversion_pigments, args_list)
            pigmedian = np.vstack([result[0] for result in results])
            pigunc = np.vstack([result[1] for result in results])
            vars_units = results[0][2]
            amps = np.vstack([result[3] for result in results])
            # restore original order
            inv_perm = np.argsort(perm)
            pigmedian = pigmedian[inv_perm,:]
            pigunc = pigunc[inv_perm,:]
            amps = amps[inv_perm,:]
            Rrs = Rrs[inv_perm,:]
            Rrs_unc = Rrs_unc[inv_perm,:]
            sstref = sstref[inv_perm]
            salref = salref[inv_perm]
        else:
        # Alternatively, run without multiprocessing for testing
            pigmedian, pigunc, vars_units, amps = rrs_inversion_pigments(Rrs, Rrs_unc, wavelength[None,:],sstref[:,None],salref[:,None],Sullivan_pure_water_temp_sal, tol_err)
        end_time_all = time.time()
        print(f"Total time for inversion of {Rrs.shape[0]} pixels: {end_time_all - start_time_all:.2f} seconds")
        print(f"Average time per pixel: {(end_time_all - start_time_all)/Rrs.shape[0]:.4f} seconds")
        output_dict = dict()
        keys = ["chla_gpig","chlc_gpig","chlb_gpig","ppc_gpig","chla_unc_gpig","chlc_unc_gpig","chlb_unc_gpig","ppc_unc_gpig"]
        fillvalue =  -32767.0
        chunk_sizes = (256, 1272)
        for key in keys:
            output_dict[key] = np.full(common_valid_mask.shape, fillvalue)
        output_dict["chla_gpig"][common_valid_mask] = pigmedian[:,0]
        output_dict["chlc_gpig"][common_valid_mask] = pigmedian[:,1]
        output_dict["chlb_gpig"][common_valid_mask] = pigmedian[:,2]
        output_dict["ppc_gpig"][common_valid_mask] = pigmedian[:,3]
        output_dict["chla_unc_gpig"][common_valid_mask] = pigunc[:,0]
        output_dict["chlc_unc_gpig"][common_valid_mask] = pigunc[:,1]
        output_dict["chlb_unc_gpig"][common_valid_mask] = pigunc[:,2]
        output_dict["ppc_unc_gpig"][common_valid_mask] = pigunc[:,3]                 
        number_of_lines = f.dimensions['number_of_lines'].size
        pixels_per_line = f.dimensions['pixels_per_line'].size
        fout.createDimension('number_of_lines', number_of_lines)
        fout.createDimension('pixels_per_line', pixels_per_line)
        fout.createDimension('wavelength', wavelength.shape[0])
        fout.createGroup('geophysical_data')
        geo_out : nc.Group = fout['geophysical_data']
        nav_out : nc.Group = fout.createGroup('navigation_data')
        long_names_dict = {
            "chla_gpig": "Chlorophyll-a concentration",
            "chlb_gpig": "Chlorophyll-b concentration",
            "chlc_gpig": "Chlorophyll-c concentration",
            "ppc_gpig": "Photoprotective carotenoids",
            "chla_unc_gpig": "Uncertainty in chlorophyll-a concentration",
            "chlb_unc_gpig": "Uncertainty in chlorophyll-b concentration",
            "chlc_unc_gpig": "Uncertainty in chlorophyll-c concentration",
            "ppc_unc_gpig":  "Uncertainty in photoprotective carotenoids",}
        standard_name_dict = {
                              "chla_gpig" : "chlorophyll_concentration_in_sea_water",
                              "chlb_gpig" : "chlorophyll_concentration_in_sea_water",
                              "chlc_gpig" : "chlorophyll_concentration_in_sea_water",
                              "ppc_gpig"  : "photoprotective_carotenoids_concentration_in_sea_water",
                              "chla_unc_gpig" : "chlorophyll_concentration_in_sea_water standard_error",
                              "chlb_unc_gpig" : "chlorophyll_concentration_in_sea_water standard_error",
                              "chlc_unc_gpig" : "chlorophyll_concentration_in_sea_water standard_error",
                              "ppc_unc_gpig"  : "photoprotective_carotenoids_concentration_in_sea_water standard_error",                             
                              }
        reference = '''Chase, A., E. Boss, I. Cetinic, and W. Slade. 2017. "Estimation of Phytoplankton
        Accessory Pigments from Hyperspectral Reflectance Spectra: Toward a
        Global Algorithm."
        Journal of Geophysical Research: Oceans, doi: 10.1002/2017JC012859.'''
        for key in output_dict:
            var = geo_out.createVariable(key, 'f4', ('number_of_lines', 'pixels_per_line'), zlib=True, fill_value=fillvalue,chunksizes=chunk_sizes)
            var[:] = output_dict[key]
            var.units = "mg m^-3"
            var.references = reference
            var.coordinate_reference = "longitude latitude"
            var.long_name = long_names_dict[key]
            var.standard_name = standard_name_dict[key]
            var.valid_min = 0
            var.valid_max = 100
        wavelength_out : nc.Variable = geo_out.createVariable("wavelength", 'f4', ('wavelength',), zlib=True, fill_value=fillvalue)
        wavelength_out[:] = wavelength
        wavelength_out.units = "nm"
        navigation_data : nc.Group = f["navigation_data"]
        latitude : nc.Variable = navigation_data["latitude"]
        longitude : nc.Variable = navigation_data["longitude"]
        lat_out : nc.Variable = nav_out.createVariable("latitude", 'f4', ('number_of_lines', 'pixels_per_line'), zlib=True, fill_value=fillvalue,chunksizes=chunk_sizes)
        lon_out : nc.Variable = nav_out.createVariable("longitude", 'f4', ('number_of_lines', 'pixels_per_line'), zlib=True, fill_value=fillvalue,chunksizes=chunk_sizes)
        lat_out[:] = latitude[:]
        lon_out[:] = longitude[:]
        lat_out.units = "degrees_north"
        lat_out.standard_name = "latitude"
        lat_out.valid_min = -90
        lat_out.valid_max = 90
        lon_out.units = "degrees_east"
        lon_out.standard_name = "longitude"
        lon_out.valid_min = -180
        lon_out.valid_max = 180
        fout.history = " ".join(sys.argv)
        time_coverage_start  = f.getncattr("time_coverage_start")
        time_coverage_end    = f.getncattr("time_coverage_end")
        fout.time_coverage_start = time_coverage_start
        fout.time_coverage_end   = time_coverage_end
        platform : str = f.getncattr("platform")
        fout.platform = platform
        instrument : str = f.getncattr("instrument")
        fout.instrument = instrument
        fout.title = f"{instrument} Level-2 GPIG"
        fout.suite = "GPIG"
        Conventions : str = f.getncattr("Conventions") 
        fout.Conventions = Conventions
        naming_authority : str = f.getncattr("naming_authority")
        fout.naming_authority = naming_authority
        geospatial_bounds : str = f.getncattr("geospatial_bounds")
        fout.geospatial_bounds = geospatial_bounds
        creator_name : str = f.getncattr("creator_name")
        fout.creator_name = creator_name
        creator_email : str = f.getncattr("creator_email")
        fout.creator_email = creator_email
        creator_url : str = f.getncattr("creator_url")
        fout.creator_url = creator_url
        project : str = f.getncattr("project")
        fout.project = project
        publisher_name : str = f.getncattr("publisher_name")
        fout.publisher_name = publisher_name
        publisher_email : str = f.getncattr("publisher_email")
        fout.publisher_email = publisher_email
        publisher_url : str = f.getncattr("publisher_url")
        fout.publisher_url = publisher_url
        geospatial_lat_max : float = f.getncattr("geospatial_lat_max")
        fout.geospatial_lat_max = geospatial_lat_max
        geospatial_lat_min : float = f.getncattr("geospatial_lat_min")
        fout.geospatial_lat_min = geospatial_lat_min
        geospatial_lon_max : float = f.getncattr("geospatial_lon_max")
        fout.geospatial_lon_max = geospatial_lon_max
        geospatial_lon_min : float = f.getncattr("geospatial_lon_min")
        fout.geospatial_lon_min = geospatial_lon_min
        processing_level : str = f.getncattr("processing_level")
        fout.processing_level = processing_level
        institution : str = f.getncattr("institution")
        fout.institution = institution
        earth_sun_distance_correction : float = f.getncattr("earth_sun_distance_correction")
        fout.earth_sun_distance_correction = earth_sun_distance_correction
        fout.date_created = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        license : str = f.getncattr("license")
        fout.license = license
        startDirection : str = f.getncattr("startDirection")
        fout.startDirection = startDirection
        endDirection : str = f.getncattr("endDirection")
        fout.endDirection = endDirection
        day_night_flag : str = f.getncattr("day_night_flag")
        fout.day_night_flag = day_night_flag
        l2_flags : nc.Variable = geophysical_data["l2_flags"]
        l2_flags_out : nc.Variable = geo_out.createVariable("l2_flags", 'i4', ('number_of_lines', 'pixels_per_line'), zlib=True,chunksizes=chunk_sizes)
        l2_flags_out[:] = l2_flags[:]
        l2_flags_out.long_name = "Level-2 Processing Flags"
        l2_flags_out.coordinate_reference = "longitude latitude"
        l2_flags_out.flag_meanings = l2_flags.flag_meanings
        l2_flags_out.flag_masks = l2_flags.flag_masks
        product_name = os.path.basename(path_out)
        fout.product_name = product_name
        geospatial_bounds_crs = f.getncattr("geospatial_bounds_crs")
        fout.geospatial_bounds_crs = geospatial_bounds_crs
        if keywords := f.ncattrs().count("keywords") > 0:
            keywords : str = f.getncattr("keywords")
            fout.keywords = keywords
            keywords_vocabulary : str = f.getncattr("keywords_vocabulary")
            fout.keywords_vocabulary = keywords_vocabulary
        fout.DOI = args.doi
        fout.processing_version = args.processing_version
        processing_control_group : nc.Group = fout.createGroup("processing_control")
        processing_control_group.software = PROGRAM
        processing_control_group.version = VERSION
        input_parameters : nc.Group = processing_control_group.createGroup("input_parameters")
        input_parameters.ifile = args.ifile
        input_parameters.sst = args.sst
        input_parameters.sal = args.sal
        input_parameters.ncores = ncores
        input_parameters.ofile = args.ofile
        input_parameters.correction = args.correction
        input_parameters.tol = tol_err
        
