import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
import pandas as pd

from mgefit.find_galaxy import find_galaxy
from mgefit.sectors_photometry import sectors_photometry
from mgefit.mge_fit_sectors import mge_fit_sectors
from mgefit.mge_print_contours import mge_print_contours
import jampy as jam
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless environments

# ----------------------------------------------------------------------------
def process_galaxy(file_path, psf_fwhm, output_dir, ngauss=3):
    hdu = fits.open(file_path)
    img = hdu[0].data * 1e10
    header = hdu[0].header
    scale = abs(header['CDELT2']) * 3600  # pixel to arcsec

    skylev = np.median(img)
    img -= skylev

    sigmapsf = psf_fwhm / 2.355
    normpsf = 1


    f = find_galaxy(img, fraction=0.04, plot=0)
    y, x = np.ogrid[:img.shape[0], :img.shape[1]]
    galaxy_mask = (x - f.xpeak) ** 2 + (y - f.ypeak) ** 2 < (2 * f.majoraxis) ** 2
    sky_pixels = img[~galaxy_mask]
    skysigma = np.std(sky_pixels)
    minlevel = skysigma

    print(f" Galaxy: {file_path.stem} | Sky Level: {skylev:.3f}, Sigma: {skysigma:.3f}, PSF_FWHM: {psf_fwhm:.3f}")

    # Sector photometry
    s = sectors_photometry(img, f.eps, f.theta, f.xpeak, f.ypeak, minlevel=minlevel, plot=0)

    # MGE fit
    m = mge_fit_sectors(s.radius, s.angle, s.counts, f.eps,
                        ngauss=ngauss, scale=scale, plot=0,
                        sigmapsf=sigmapsf, normpsf=normpsf)

    total_counts, sigma, q_obs = m.sol
    surf = total_counts / (2 * np.pi * q_obs * sigma**2)
    reff, reff_maj, eps_e, lum_tot = jam.mge.half_light_isophote(surf, sigma * scale, q_obs)

    galaxy_name = file_path.stem
    output_data_path = output_dir / f"{galaxy_name}_radii.txt"

    with open(output_data_path, 'w') as file:
        file.write(f"Galaxy: {galaxy_name}\n")
        file.write(f"Effective Radius (Reff): {reff:.3f} arcsec\n")
        file.write(f"Effective Major Axis Radius (Reff_maj): {reff_maj:.3f} arcsec\n")
        file.write(f"Eccentricity: {eps_e:.3f}\n")
        file.write(f"Total Luminosity: {lum_tot:.3f}\n")

    # Save contour plot
    output_image_path = output_dir / f"{galaxy_name}_contours.png"
    plt.figure(figsize=(10, 5))
    mge_print_contours(img, f.theta, f.xpeak, f.ypeak, m.sol, scale=scale,
                       sigmapsf=sigmapsf, normpsf=normpsf, binning=1, minlevel=minlevel)
    plt.title(f"Contours of {galaxy_name}")
    plt.savefig(output_image_path)
    plt.close()

    print(f" {galaxy_name}: Reff = {reff:.3f} arcsec")

    return {
        "Galaxy": galaxy_name,
        "Reff (arcsec)": reff,
        "Reff_maj (arcsec)": reff_maj,
        "Eccentricity": eps_e,
        "Luminosity": lum_tot,
        "fwhm": psf_fwhm,
    }

# ----------------------------------------------------------------------------
def batch_process(psf_csv_path, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df_psf = pd.read_csv(psf_csv_path)
    fits_dir = Path("/Users/spiniello/Desktop/DESI/PS1_visible_fits")
    summary = []

    for _, row in df_psf.iterrows():
        target_id = str(row["TARGETID"])
        file_path = fits_dir / f"{target_id}_r.fits"

        if not file_path.exists():
            print(f" File not found: {file_path}")
            continue

        psf_fwhm = row["psf_fwhm"]

        try:
            result = process_galaxy(file_path, psf_fwhm, output_dir)
            if result:
                summary.append(result)
        except Exception as e:
            print(f" Error processing {file_path.name}: {e}")

    # Save summary CSV
    if summary:
        df_summary = pd.DataFrame(summary)
        csv_path = output_dir / "ps1_mge_radii_summary.csv"
        df_summary.to_csv(csv_path, index=False)
        print(f" Summary CSV saved to: {csv_path}")

# ----------------------------------------------------------------------------
if __name__ == "__main__":
    psf_csv = Path("UCMGs_VISIBLE_PSF_with_paths.csv")
    output_dir = Path("/Users/spiniello/Desktop/DESI/PS1_MGE_Results")
    batch_process(psf_csv, output_dir)