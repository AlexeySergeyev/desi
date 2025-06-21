# %%
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.config import set_temp_cache
import os

# %% [markdown]
# # Get stellar data from the file

# %%
if not os.path.exists('data/dr1_galaxy_lowZ_stellarmass_UCMGs_10.3.fits'):
    raise FileNotFoundError("The stellar data file does not exist. Please download it first.")
    exit(1)

hdu_stellar = fits.open('data/dr1_galaxy_lowZ_stellarmass_UCMGs_10.3.fits')
hdu_stellar.info()

# %%
hdu_stellar[1].header

# %%
# Convert FITS BINTABLE to astropy Table
stellar_table = Table(hdu_stellar[1].data)

# Display basic info about the table
print(f"Number of rows: {len(stellar_table)}")
print(f"Number of columns: {len(stellar_table.colnames)}")
print(f"Column names: {stellar_table.colnames}")
print("\nFirst few rows:")
stellar_table[:5]

# %%
df_stellar = stellar_table.to_pandas()
# df_stellar

# %%
# df_stellar['HEALPIX'].unique().shape

# %% [markdown]
# # Get spectrum from DESI

# %%
# Define the spectroscopic product directory
desi_root = "https://data.desi.lbl.gov/public/dr1//spectro/redux/iron"


# %%
def save_spectrum_to_fits(targetid, wave, flux, ivar, mask):
    """
    Save the spectrum data to a FITS file.
    
    Parameters:
    - targetid: Unique identifier for the target
    - wave: Dictionary with wavelength data for each camera
    - flux: Dictionary with flux data for each camera
    - ivar: Dictionary with inverse variance data for each camera
    - mask: Dictionary with mask data for each camera
    """
    # Ensure the output directory exists
    out_folder = './data/spectra'
    if not os.path.exists(out_folder):
        print(f"Creating output directory: {out_folder}")
        os.makedirs(out_folder, exist_ok=True)
    
    # Define the output path
    out_path = os.path.join(out_folder, f'{targetid}.fits')
    
    # Create HDUList with PrimaryHDU and BinTableHDUs for each camera
    hdus = []
    hdus.append(fits.PrimaryHDU())
    
    for cam in ['B', 'R', 'Z']:
        cols = [
            fits.Column(name='WAVELENGTH', format='D', array=wave[cam]),
            fits.Column(name='FLUX'      , format='D', array=flux[cam]),
            fits.Column(name='IVAR'      , format='D', array=ivar[cam]),
            fits.Column(name='MASK'      , format='J', array=mask[cam]),
        ]
        table_hdu = fits.BinTableHDU.from_columns(cols, name=f'SPECTRUM_{cam}')
        hdus.append(table_hdu)
    
    # Write the HDUList to a FITS file
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(out_path, overwrite=True)
    print(f"Saved spectrum for TARGETID {targetid} to {out_path}")

# %%
def load_fits_spectrum(targetid):
    """
    Load the spectrum data from a FITS file.
    
    Parameters:
    - targetid: Unique identifier for the target
    
    Returns:
    - wave: Dictionary with wavelength data for each camera
    - flux: Dictionary with flux data for each camera
    - ivar: Dictionary with inverse variance data for each camera
    - mask: Dictionary with mask data for each camera
    """
    out_folder = './data/spectra'
    out_path = os.path.join(out_folder, f'{targetid}.fits')
    
    if not os.path.exists(out_path):
        raise FileNotFoundError(f"Spectrum file for TARGETID {targetid} not found.")
    
    hdu_spectrum = fits.open(out_path)
    
    wave = dict()
    flux = dict()
    ivar = dict()
    mask = dict()
    
    for cam in ['B', 'R', 'Z']:
        wave[cam] = hdu_spectrum[f'SPECTRUM_{cam}'].data['WAVELENGTH']
        flux[cam] = hdu_spectrum[f'SPECTRUM_{cam}'].data['FLUX']
        ivar[cam] = hdu_spectrum[f'SPECTRUM_{cam}'].data['IVAR']
        mask[cam] = hdu_spectrum[f'SPECTRUM_{cam}'].data['MASK']
    
    return wave, flux, ivar, mask


# %%
def get_values(hdu_coadd, target_index):
    """
    Extracts wavelength, flux, inverse variance, and mask data for a given target index
    from the coadd HDU.
    """
    # Extract data for the specified target index
    wave = dict()
    flux = dict()
    ivar = dict()
    mask = dict()
    for camera in ['B', 'R', 'Z']:
        wave[camera] = hdu_coadd[f'{camera}_WAVELENGTH'].data
        flux[camera] = hdu_coadd[f'{camera}_FLUX'].data[target_index]
        ivar[camera] = hdu_coadd[f'{camera}_IVAR'].data[target_index]
        mask[camera] = hdu_coadd[f'{camera}_MASK'].data[target_index]
    return wave, flux, ivar, mask


# %%
def plot_spectrum(wave, flux, ivar, mask, targetid, redshift, issave=False):
    """
    Plots the spectrum for a given target.
    """
    plt.figure(figsize=(8, 3))
    for camera in ['B', 'R', 'Z']:
        w = wave[camera]
        f = flux[camera]
        i = ivar[camera]
        m = np.bool(mask[camera])
        plt.plot(w[~m], f[~m], label=f'{camera} band', lw=0.5)
        # plt.fill_between(w[~m], f[~m] - 1/np.sqrt(i[~m]), f[~m] + 1/np.sqrt(i[~m]), alpha=0.2)
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Flux (arbitrary units)')
    plt.title(f'Spectrum for TARGETID {targetid} at redshift {redshift:.3f}')
    plt.legend()
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    if issave:
        if not os.path.exists('./figs/spectra/'):
            os.makedirs('./figs/spectra/')
        plt.savefig(f'./figs/spectra/{targetid}.png', dpi=300)
        plt.close()
    else:
        plt.show()


# %%
start = 40
num = 20
for i, row in df_stellar[start:start + num].iterrows():

    redshift = row['Z']
    targetid = row['TARGETID']
    survey = row['SURVEY']
    program = row['PROGRAM']
    healpix = row['HEALPIX']
    hpixgroup = healpix // 100
    print(f"Row {i}: SURVEY={survey}, PROGRAM={program}, HEALPIX={healpix}, HPixGroup={hpixgroup}")
    # Filename
    coadd_filepath = f'{desi_root}/healpix/{survey}/{program}/{hpixgroup}/{healpix}/coadd-{survey}-{program}-{healpix}.fits'
    
    if not os.path.exists('./data/stellars/'):
        os.makedirs('./data/stellars/')

    with set_temp_cache(path='./data/stellars/', delete=False):
        hdu_coadd = fits.open(coadd_filepath, cache=True)
    # hdu_coadd.info()

    df_fiber = pd.DataFrame(hdu_coadd['FIBERMAP'].data)
    cond = df_fiber['TARGETID'] == targetid
    target_index = df_fiber.index[cond][0]  # Get the first index where condition is True

    wave, flux, ivar, mask = get_values(hdu_coadd, target_index)

    plot_spectrum(wave, flux, ivar, mask, targetid, redshift, issave=True)
    # Save the spectrum to a FITS file
    save_spectrum_to_fits(targetid, wave, flux, ivar, mask)
    # Break after processing the first target for demonstration purposes
  