#!/usr/bin/env python3
"""
Compute and plot dependences of observed galaxy radii vs. intrinsic R_e
for four surveys (SDSS, Pan-STARRS1, KiDS, Euclid) at their typical
PSF FWHM and pixel scales.

Metrics plotted:
  - Observed half-light radius R50
  - Petrosian50 and Petrosian90 (within 2 r_P, classical Petrosian eta=0.2)

Method (no noise):
  1) Render a 2D Sérsic profile on a fine grid.
  2) Convolve with a circular Gaussian PSF.
  3) Measure circular profiles in annuli with dr = min(fine_pixel, pixel_scale/2).
  4) Extract radii from cumulative light or Petrosian definition.

Notes:
  * Uses matplotlib only (no seaborn). One chart per figure. No explicit colors set.
  * SciPy is optional (for gaussian_filter); otherwise FFT Gaussian is used.
  * You may want to increase "extent_factor" or samples for very extended profiles.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt, log
from pathlib import Path

# --- Optional SciPy convolution (faster), with fallback ---
try:
    from scipy.ndimage import gaussian_filter
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False

# ----------------------------- Sérsic helpers ------------------------------
def b_n(n: float) -> float:
    """Ciotti & Bertin (1999) asymptotic expansion, good for 0.5<n<10."""
    return 2*n - 1/3 + 4/(405*n) + 46/(25515*n**2) + 131/(1148175*n**3)

def sersic_I(r, Re, n, Ie=1.0):
    bn = b_n(n)
    return Ie * np.exp(-bn*((r/Re)**(1.0/n) - 1.0))

# --------------------------- Image + PSF helpers ---------------------------
def make_image(Re_arcsec, n, fwhm_arcsec, hires_pix, extent_factor=5.0):
    """Render a centered 2D Sérsic image on a fine grid."""
    L = extent_factor * max(Re_arcsec, fwhm_arcsec)
    N = int(np.ceil(2*L/hires_pix)) | 1  # enforce odd size
    ax = (np.arange(N) - (N//2)) * hires_pix
    xx, yy = np.meshgrid(ax, ax, indexing='xy')
    rr = np.hypot(xx, yy)
    img = sersic_I(rr, Re_arcsec, n, 1.0)
    return img, hires_pix

def convolve_psf(img, fwhm_pix):
    """Convolve with a circular Gaussian PSF of given FWHM in pixels."""
    if fwhm_pix <= 0:
        return img.copy()
    sigma = fwhm_pix / (2*sqrt(2*log(2)))
    if SCIPY_OK:
        return gaussian_filter(img, sigma=sigma, mode='nearest')
    # FFT fallback
    N = img.shape[0]
    ky = np.fft.fftfreq(N)[:, None]
    kx = np.fft.fftfreq(N)[None, :]
    H = np.exp(-2*(np.pi**2)*(sigma**2)*(kx**2 + ky**2))
    return np.fft.ifft2(np.fft.fft2(img) * H).real

# ---------------------- Radial profile + radii metrics ---------------------
def radial_profile(img, pixscale, dr):
    """
    Compute circularly averaged annular surface brightness and cumulative flux.
    Returns dict with arrays: r, I_ann, eta (Petrosian), cum (cumulative flux), tot.
    """
    N = img.shape[0]
    ax = (np.arange(N) - (N//2)) * pixscale
    xx, yy = np.meshgrid(ax, ax, indexing='xy')
    rr = np.hypot(xx, yy).ravel()
    vals = img.ravel()

    rmax = rr.max()
    nbins = max(10, int(np.ceil(rmax/dr)))
    edges = np.linspace(0, nbins*dr, nbins+1)
    inds  = np.digitize(rr, edges) - 1
    sum_I = np.bincount(inds, weights=vals, minlength=nbins)
    count = np.bincount(inds, minlength=nbins).astype(float)
    area  = count * (pixscale**2)

    r_cent = (edges[:-1] + edges[1:])/2

    # annular mean SB and cumulative light
    with np.errstate(invalid='ignore', divide='ignore'):
        I_ann = sum_I / area
    I_ann = np.nan_to_num(I_ann)

    cum_flux = np.cumsum(sum_I)
    total_flux = cum_flux[-1]

    # Mean SB within R for Petrosian eta
    area_cum = np.cumsum(area)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean_inside = np.divide(cum_flux, area_cum, out=np.zeros_like(cum_flux), where=area_cum>0)
        eta = np.divide(I_ann, mean_inside, out=np.zeros_like(I_ann), where=mean_inside>0)

    return dict(r=r_cent, I_ann=I_ann, eta=eta, cum=cum_flux, tot=total_flux)

def interp_radius_at_y(r, y, y0):
    """Linear interpolation to find r where y crosses y0 (first crossing)."""
    dif = y - y0
    s = np.sign(dif[:-1]) * np.sign(dif[1:])
    idx = np.where((s <= 0) & np.isfinite(y[:-1]) & np.isfinite(y[1:]))[0]
    if len(idx) == 0:
        # fallback: nearest
        j = int(np.argmin(np.abs(dif)))
        return r[j]
    j = idx[0]
    y1, y2 = y[j], y[j+1]
    if y2 == y1:
        return r[j]
    t = (y0 - y1) / (y2 - y1)
    return r[j] + t*(r[j+1]-r[j])

def radius_at_cumfrac(r, cum, frac):
    """Radius where cumulative flux reaches given fraction of total."""
    target = frac * cum[-1]
    if target <= 0:
        return 0.0
    j = np.searchsorted(cum, target)
    if j == 0:
        return r[0]
    f1, f2 = cum[j-1], cum[j]
    if f2 == f1:
        return r[j]
    t = (target - f1)/(f2 - f1)
    return r[j-1] + t*(r[j]-r[j-1])

def measure_radii(img, pixscale, dr, eta0=0.2):
    """Return dict of observed R50, Petrosian r_P, rP50, rP90 (within 2 r_P)."""
    prof = radial_profile(img, pixscale, dr)
    r = prof['r']; cum = prof['cum']; tot = prof['tot']
    R50  = radius_at_cumfrac(r, cum, 0.5)
    rP   = interp_radius_at_y(r, prof['eta'], eta0)
    cum_2rP = np.interp(2.0*rP, r, cum)
    rP50 = radius_at_cumfrac(r, cum, 0.5 * (cum_2rP / tot))
    rP90 = radius_at_cumfrac(r, cum, 0.9 * (cum_2rP / tot))
    return dict(R50=R50, rP=rP, rP50=rP50, rP90=rP90)

# ------------------------------ Configuration ------------------------------
# (Survey name, pixel scale ["/pix], FWHM PSF ["])
SURVEYS = [
    ("SDSS r",        0.396, 1.30),
    ("Pan-STARRS1 r", 0.258, 1.00),
    ("KiDS r",        0.214, 0.72),
    ("Euclid VIS",    0.101, 0.17),
]

SERSIC_N = [1, 2, 4]
RE_OVER_FWHM = np.geomspace(0.25, 5.0, 16)  # from PSF-limited to well-resolved
FINE_PSF_SAMPLING = 30.0  # fine rendering: pixels per FWHM
EXTENT_FACTOR = 5.0       # image size ~ EXTENT_FACTOR * max(Re, FWHM)

# ------------------------------- Main routine ------------------------------
def run():
    rows = []
    for survey_name, pixscale, fwhm in SURVEYS:
        for n in SERSIC_N:
            for f in RE_OVER_FWHM:
                Re = f * fwhm  # intrinsic effective radius [arcsec]
                hires = fwhm / FINE_PSF_SAMPLING

                # Intrinsic (no PSF) for reference
                img_intr, _ = make_image(Re, n, fwhm, hires_pix=hires, extent_factor=EXTENT_FACTOR)
                dr = min(hires, pixscale/2.0)
                m_intr = measure_radii(img_intr, hires, dr)

                # Observed with PSF
                img_obs = convolve_psf(img_intr, fwhm_pix=fwhm/hires)
                m_obs  = measure_radii(img_obs, hires, dr)

                rows.append(dict(
                    Survey=survey_name, n=n, FWHM=fwhm, pixscale=pixscale,
                    Re_intrinsic_arcsec=Re, Re_over_FWHM=f,
                    R50_obs=m_obs['R50'],   R50_intr=m_intr['R50'],
                    rP_obs=m_obs['rP'],
                    rP50_obs=m_obs['rP50'], rP50_intr=m_intr['rP50'],
                    rP90_obs=m_obs['rP90'], rP90_intr=m_intr['rP90'],
                ))

    df = pd.DataFrame(rows)
    outdir = Path("outputs")
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "radii_dependences_by_survey.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved results to {csv_path}")

    # --- Plotting ---
    def plot_metric(metric_key, title, ylabel):
        for n in sorted(df['n'].unique()):
            sub = df[df['n']==n]
            plt.figure()
            for survey_name, _, _ in SURVEYS:
                ss = sub[sub['Survey']==survey_name]
                x = ss['Re_intrinsic_arcsec'].values
                y = ss[metric_key].values
                plt.plot(x, y, label=survey_name)
            if metric_key == 'R50_obs':
                # 1:1 reference line
                xx = np.linspace(sub['Re_intrinsic_arcsec'].min(), sub['Re_intrinsic_arcsec'].max(), 200)
                plt.plot(xx, xx, linestyle='--', label='1:1 reference')
            plt.xlabel(r'Intrinsic $R_e$ [arcsec]')
            plt.ylabel(ylabel)
            plt.title(f"{title} (Sérsic n={n})")
            plt.legend()
            plt.grid(True)
            fig_path = outdir / f"{metric_key}_vs_Re_n{n}.png"
            plt.savefig(fig_path, dpi=160, bbox_inches='tight')
            print(f"Saved {fig_path}")

    plot_metric('R50_obs',  r'Observed half-light radius $R_{50}$ vs. intrinsic $R_e$', r'Observed $R_{50}$ [arcsec]')
    plot_metric('rP50_obs', r'Observed Petrosian50 vs. intrinsic $R_e$',               r'Observed $r_{P50}$ [arcsec]')
    plot_metric('rP90_obs', r'Observed Petrosian90 vs. intrinsic $R_e$',               r'Observed $r_{P90}$ [arcsec]')

if __name__ == "__main__":
    run()
