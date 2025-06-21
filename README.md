# DESI Stellar Spectra Analysis

This repository contains tools for analyzing stellar spectra from the Dark Energy Spectroscopic Instrument (DESI) Data Release 1 (DR1).

## Overview

The main notebook `stellar.ipynb` provides functionality to:
- Load and process DESI stellar catalog data
- Download and extract individual stellar spectra
- Visualize spectra across B, R, and Z camera bands
- Save spectra to FITS files for further analysis

## Data Requirements

The analysis requires the following DESI DR1 data files:
- `data/dr1_galaxy_lowZ_stellarmass_UCMGs_10.3.fits` - Stellar catalog with target information
- Access to DESI DR1 spectroscopic data via `https://data.desi.lbl.gov/public/dr1/spectro/redux/iron`

## Dependencies

Required Python packages:
```python
matplotlib
numpy
astropy
pandas
```

Install dependencies:
```bash
pip install matplotlib numpy astropy pandas
```

## Usage

### 1. Load Stellar Catalog

```python
from astropy.io import fits
from astropy.table import Table

# Load FITS catalog
hdu_stellar = fits.open('data/dr1_galaxy_lowZ_stellarmass_UCMGs_10.3.fits')

# Convert to astropy Table
stellar_table = Table(hdu_stellar[1].data)

# Convert to pandas DataFrame for easier manipulation
df_stellar = stellar_table.to_pandas()
```

### 2. Extract Individual Spectra

The notebook provides functions to download and process individual stellar spectra:

```python
# Process spectra for a range of targets
start = 20
num = 20
for i, row in df_stellar[start:start + num].iterrows():
    targetid = row['TARGETID']
    redshift = row['Z']
    survey = row['SURVEY']
    program = row['PROGRAM']
    healpix = row['HEALPIX']
    
    # Download and process spectrum
    # (see notebook for full implementation)
```

### 3. Visualization

Plot stellar spectra across all three DESI camera bands:

```python
def plot_spectrum(wave, flux, ivar, mask, targetid, redshift, issave=False):
    """
    Plots the spectrum for a given target across B, R, Z bands.
    """
    plt.figure(figsize=(8, 3))
    for camera in ['B', 'R', 'Z']:
        w = wave[camera]
        f = flux[camera]
        m = np.bool(mask[camera])
        plt.plot(w[~m], f[~m], label=f'{camera} band', lw=0.5)
    
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Flux (arbitrary units)')
    plt.title(f'Spectrum for TARGETID {targetid} at redshift {redshift:.3f}')
    plt.legend()
    plt.grid(linestyle='--', alpha=0.5)
    plt.show()
```

### 4. Save Spectra to FITS

Save processed spectra to individual FITS files:

```python
def save_spectrum_to_fits(targetid, wave, flux, ivar, mask):
    """
    Save spectrum data to a FITS file with separate extensions for each camera.
    """
    # Creates FITS file with PRIMARY HDU and SPECTRUM_B, SPECTRUM_R, SPECTRUM_Z extensions
    # (see notebook for full implementation)
```

### 5. Load Saved Spectra

Load previously saved spectrum FITS files:

```python
def load_fits_spectrum(targetid):
    """
    Load spectrum data from a saved FITS file.
    
    Returns:
    - wave, flux, ivar, mask dictionaries for each camera
    """
    # (see notebook for full implementation)
```

## Directory Structure

```
desi/
├── stellar.ipynb           # Main analysis notebook
├── README.md              # This file
├── data/
│   ├── dr1_galaxy_lowZ_stellarmass_UCMGs_10.3.fits  # Input catalog
│   ├── stellars/          # Cache directory for downloaded spectra
│   └── spectra/           # Saved individual spectrum FITS files
└── figs/
    └── spectra/           # Saved spectrum plots
```

## Example Analysis

Load and visualize a specific stellar spectrum:

```python
# Example target
test_targetid = 2842484900102147

# Get redshift from catalog
redshift = df_stellar.loc[df_stellar['TARGETID'] == test_targetid, 'Z'].values[0]

# Load saved spectrum
wave, flux, ivar, mask = load_fits_spectrum(test_targetid)

# Plot spectrum
plot_spectrum(wave, flux, ivar, mask, targetid=test_targetid, redshift=redshift, issave=False)
```

## Notes

- The notebook processes DESI DR1 data with iron pipeline reduction
- Spectra are cached locally to avoid repeated downloads
- All three camera bands (B: ~3600-5800Å, R: ~5600-7600Å, Z: ~7300-9800Å) are processed
- Masked pixels are excluded from visualization
- Output FITS files follow standard astronomical data format conventions

## License
This project is licensed under the MIT License. See the LICENSE file for details.