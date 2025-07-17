#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 12:12:34 2025

@author: spiniello
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
from scipy.optimize import curve_fit
import pandas as pd
from mgefit.find_galaxy import find_galaxy
import matplotlib
matplotlib.use('Agg')

def saturation_fit(radius, f_max, k):
    return f_max * (1 - np.exp(-k * radius))

def process_galaxy(file_path, output_dir, max_radius=8.0, step_arcsec=0.1):
    hdu = fits.open(file_path)
    img = hdu[0].data
    header = hdu[0].header

    scale = abs(header['CDELT2']) * 3600
    skylev = np.median(img)
    img -= skylev
    
    
    y, x = np.ogrid[:img.shape[0], :img.shape[1]]

    # Assume galaxy is at center of image
    y_center, x_center = np.array(img.shape) / 2
    # Target: mask ~10,000 pixels
    n_target_pixels = 10000
    default_radius = int(np.sqrt(n_target_pixels / np.pi))

    # Create circular mask centered on image
    galaxy_mask = (x - x_center)**2 + (y - y_center)**2 < default_radius**2
    n_mask_pixels = np.sum(galaxy_mask)

    # Estimate sky sigma
    sky_pixels = img[~galaxy_mask]
    skysigma = np.std(sky_pixels)

    print(f"Requested mask area: ~{n_target_pixels} pixels")
    print(f"Calculated radius: {default_radius} pixels")
    print(f"Actual mask size: {n_mask_pixels} pixels")
    print(f"Sky sigma: {skysigma:.3f}")
    
    
    #If galaxy is not in the centre, FINDS galaxy automatically
    #f = find_galaxy(img, fraction=0.04, plot=False)
    #x_center, y_center = f.xpeak, f.ypeak   
    #galaxy_mask = (x - f.xpeak) ** 2 + (y - f.ypeak) ** 2 < (2 * f.majoraxis) ** 2
    #sky_pixels = img[~galaxy_mask]
    #skysigma = np.std(sky_pixels)
    #print(f" Sky: {skylev:.3f}, Sigma: {skysigma:.3f}")
    #  Number of pixels inside the mask
    #n_mask_pixels = np.sum(galaxy_mask)
    #print(f" Number of pixels in galaxy mask: {n_mask_pixels}")

    radii_arcsec = np.arange(step_arcsec, max_radius + step_arcsec, step_arcsec)
    radii_pixel = radii_arcsec / scale

    fluxes = []
    radii_arcsec_corrected = []

    for r in radii_pixel:
        mask = (x - x_center) ** 2 + (y - y_center) ** 2 <= r ** 2
        npix_j = np.sum(mask)
        r_j = np.sqrt(npix_j / np.pi) * scale
        radii_arcsec_corrected.append(r_j)
        flux = np.sum(img[mask])
        fluxes.append(flux)

    fluxes = np.array(fluxes)
    radii_arcsec_corrected = np.array(radii_arcsec_corrected)

    try:
        saturation_mask = (radii_arcsec_corrected >= 4.0) & (radii_arcsec_corrected <= 6.0)
        radii_fit = radii_arcsec_corrected[saturation_mask]
        fluxes_fit = fluxes[saturation_mask]
        popt, _ = curve_fit(saturation_fit, radii_fit, fluxes_fit, p0=[fluxes[-1], 1.0])
        radii_fine = np.linspace(4, 8, 500)
        fit_curve = saturation_fit(radii_fine, *popt) / popt[0]
        max_flux_fit, k = popt
    except RuntimeError as e:
        print(f"Fit failed for {file_path.stem}: {e}")
        return file_path.stem, None, None

    fit_fluxes = saturation_fit(radii_arcsec_corrected, *popt)
    residuals = np.abs(fluxes - fit_fluxes)
    uncertainty = np.mean(residuals[saturation_mask])
    norm_fluxes_fit = fluxes / max_flux_fit
    effective_radius = np.interp(0.5, norm_fluxes_fit, radii_arcsec_corrected)

    plt.figure(figsize=(8, 6))
    plt.plot(radii_arcsec_corrected, norm_fluxes_fit, label='Normalized Flux Data', color='blue', marker='o', linestyle='-')
    plt.plot(radii_fine, fit_curve, label='Saturation Fit (Extrapolated)', color='red', linestyle='--')
    plt.axhline(0.5, color='purple', linestyle='--', label='50% Flux Level')
    plt.axvline(effective_radius, color='green', linestyle=':', label=f'Reff = {effective_radius:.2f} arcsec')
    plt.xlabel('Radius (arcsec)', fontsize=14)
    plt.ylabel('Normalized Flux', fontsize=14)
    plt.title(f'Cumulative Flux vs Radius: {file_path.stem}', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plot_path = output_dir / f"{file_path.stem}_flux_profile.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

    return file_path.stem, effective_radius, uncertainty

def batch_process(psf_csv_path, output_dir, output_csv):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df_psf = pd.read_csv(psf_csv_path)
    fits_base_dir = Path("/Users/spiniello/Desktop/DESI/PS1_visible_fits")
    results = []

    for _, row in df_psf.iterrows():
        target_id = str(row["TARGETID"])
        file_path = fits_base_dir / f"{target_id}_r.fits"

        if not file_path.exists():
            print(f" File not found: {file_path}")
            continue

        try:
            galaxy_name, effective_radius, uncertainty = process_galaxy(file_path, output_dir)
            if effective_radius is not None:
                results.append({
                    "Galaxy": galaxy_name,
                    "Effective Radius (arcsec)": effective_radius,
                    "Uncertainty": uncertainty,
                    "PSF_FWHM (arcsec)": row["psf_fwhm"]
                })
        except Exception as e:
            print(f" Error processing {file_path.name}: {e}")

    df_out = pd.DataFrame(results)
    df_out.to_csv(output_csv, index=False)
    print(f" Results saved to {output_csv}")

# --- Main ---
if __name__ == "__main__":
    psf_csv = Path("UCMGs_VISIBLE_PSF_with_paths.csv")
    output_dir = Path("/Users/spiniello/Desktop/DESI/PS1_Flux_Results")
    output_csv = Path("/Users/spiniello/Desktop/DESI/ps1_effective_radii.csv")
    batch_process(psf_csv, output_dir, output_csv)
