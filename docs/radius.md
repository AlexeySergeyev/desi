# Comparing Petrosian, Sérsic Effective, de Vaucouleurs, and Kron Radii

## Introduction

Defining a single “size” for a galaxy is non-trivial. Common choices include:
	•	Petrosian radius $R_P$: an observational, profile-driven metric radius from a surface-brightness ratio.
	•	Sérsic effective radius $R_e$: the (model) half-light radius from a PSF-convolved Sérsic fit.
	•	de Vaucouleurs radius: $R_e$ for a Sérsic profile with $n=4$ (the classic $R^{1/4}$ law).
	•	Kron radius $r_K$: a fixed multiple of the first-moment radius $r_1$.

Below are consistent definitions, equations (each with a brief description), analytic comparisons for idealized profiles, survey implementations (SDSS, Pan-STARRS, KiDS, Euclid), PSF-convolution effects (including the near-resolution limit), a quick-reference table, and rules of thumb.

⸻

## 1) Definitions

### 1.1 Sérsic profile (base model)

$$
I(r)=I_e \exp!\left[-,b_n!\left(\left(\frac{r}{R_e}\right)^{1/n}-1\right)\right],
\qquad
b_n \approx 2n-\frac13+\frac{4}{405n}+\cdots
$$
Description: Intensity profile for a Sérsic model (left), with an accurate asymptotic approximation for the structural constant $b_n$ (right) ensuring that $R_e$ encloses half the total model light.
	•	$R_e$: (model) half-light radius.
	•	$n$: Sérsic index ($n=1$ disk-like exponential; $n=4$ de Vaucouleurs).

Cumulative light fraction:
$$
f(<r)=\frac{\gamma!\big(2n,, b_n (r/R_e)^{1/n}\big)}{\Gamma(2n)}.
$$
Description: Fraction of total Sérsic light enclosed within $r$, expressed via the lower incomplete gamma function (normalised by $\Gamma(2n)$).

### 1.2 Petrosian radius

General ratio:
$$
\eta(r)\equiv \frac{I(r)}{\langle I\rangle(<r)}
= \frac{I(r)}{\displaystyle \frac{1}{\pi r^2}\int_0^r 2\pi t, I(t),dt}.
$$
Description: Definition of the Petrosian ratio as the local surface brightness divided by the mean interior surface brightness within $r$.

SDSS/PS1 operational ratio (mean in an annulus) and threshold:
$$
\eta(r)=\frac{\langle I\rangle_{[0.8r,,1.25r]}}{\langle I\rangle(<r)},
\qquad
R_P:\ \eta(R_P)=0.2 .
$$
Description: Practical SDSS/Pan-STARRS implementation: the local term is measured in an annulus; $R_P$ is the radius where the ratio equals 0.2.

Petrosian flux and radii:
$$
F_P=\int_0^{2R_P} 2\pi t,I(t),dt, \qquad
f(<R_{50})=\tfrac12 f(<2R_P),\quad
f(<R_{90})=0.9, f(<2R_P).
$$
Description: Petrosian flux is the light within $2R_P$; $R_{50}$ and $R_{90}$ enclose 50% and 90% of that Petrosian flux, respectively.

### 1.3 Sérsic effective radius $R_e$

$R_e$ is the radius enclosing 50% of the total Sérsic model light (i.e., extrapolated to infinity). In practice, surveys fit PSF-convolved models in 2D to infer intrinsic $R_e$ and $n$.

### 1.4 de Vaucouleurs radius

The de Vaucouleurs law is the $n=4$ Sérsic case. The “de Vaucouleurs radius” is $R_e$ from an $n=4$ fit. SDSS also reports best-fit $r_{\rm deV}$ (PSF-convolved $n=4$) and $r_{\rm exp}$ (PSF-convolved $n=1$).

### 1.5 Kron radius

First-moment and Kron radii:
$$
r_1=\frac{\displaystyle\int_0^\infty 2\pi r^2 I(r),dr}{\displaystyle\int_0^\infty 2\pi r, I(r),dr},
\qquad
r_K=k,r_1\quad (k\simeq 2\text{–}2.5;\ \text{SExtractor AUTO uses }k=2.5).
$$
Description: The first-moment radius $r_1$ measures the flux-weighted mean radius; the Kron radius $r_K$ scales this by a constant factor $k$ to define a photometric aperture.

Closed forms for Sérsic:
$$
\boxed{\frac{r_1}{R_e}=\frac{\Gamma(3n)}{\Gamma(2n)}, b_n^{-n}},
\qquad
\boxed{r_K=k,R_e,\frac{\Gamma(3n)}{\Gamma(2n)}, b_n^{-n}}.
$$
Description: Exact relations linking $r_1$ and $r_K$ to the Sérsic $R_e$ and $n$ through gamma functions; useful for predicting Kron apertures from fitted Sérsic parameters.

⸻

## 2) Analytic comparisons (idealized, no PSF)

Benchmarks for SDSS/PS1 Petrosian definition:

- **Exponential ($n=1$)**  
  - $b_1 \simeq 1.678$, $\displaystyle r_1/R_e \simeq 2/1.678 \approx 1.192 \Rightarrow r_K(2.5) \approx 2.98 R_e$  
  - Petrosian: $\displaystyle R_P \approx 2.1129 R_e$, $f(<2R_P) \approx 0.9933$  
  - Hence, $r_K/R_P \approx 1.41$

- **de Vaucouleurs ($n=4$)**  
  - $b_4 \simeq 7.669$, $\displaystyle r_1/R_e \approx 2.289 \Rightarrow r_K(2.5) \approx 5.72 R_e$  **
  - Petrosian: $\displaystyle R_P \approx 1.7140 R_e$, $f(<2R_P) \approx 0.8165$  
  - Hence, $r_K/R_P \approx 3.34$

**Trend with $n$**:  
- $R_P/R_e$ decreases monotonically with increasing $n$.  
- Kron ($k=2.5$) apertures are typically larger than $R_P$, especially at high $n$.  
- Flux within Kron(2.5):  
  - $f \approx 0.96$ ($n=1$)  
  - $f \approx 0.90$ ($n=4$ ).

⸻

## 3) Survey implementations

### Survey Implementations

**SDSS**
- Computes the Petrosian radius $R_P$ in the $r$-band using $\eta(R_P)=0.2$; measures Petrosian flux within $2R_P$.
- Reports Petrosian $R_{50}$ and $R_{90}$ (circular apertures).
- Provides model radii from PSF-convolved fits: $r_{\rm exp}$ ($n=1$) and $r_{\rm deV}$ ($n=4$).

**Pan-STARRS (PS1)**
- Adopts the SDSS Petrosian scheme ($\eta=0.2$, flux within $2R_P$).
- Implementation details (e.g., annulus spacing, sky subtraction, masking) differ, but the definitions are consistent with SDSS.

**KiDS**
- Uses PSF-convolved single-Sérsic 2D fits (e.g., 2DPHOT).
- Reports $R_{e,\rm maj}$ (major-axis half-light radius), Sérsic index $n$, axis ratio $q$, etc.
- Circularized size is given by $R_{e,\rm circ}=R_{e,\rm maj}\sqrt{q}$.

**Euclid**
- Employs PSF-convolved Sérsic fits (e.g., SourceXtractor++, GALFIT).
- Reports $R_e$, $n$, ellipticity, and magnitudes.
- Multi-component fits are planned for future data releases.

**Summary:**  
- **SDSS/PS1:** Empirical, Petrosian-based apertures.
- **KiDS/Euclid:** Model-based $R_e$ from PSF-aware Sérsic fits.

⸻

5) Quick reference (no PSF)

Table 1: Benchmarks for ideal Sérsic profiles without PSF, using SDSS/PS1 Petrosian ($\eta=0.2$) and Kron with $k=2.5$.

<div align="center">

| Quantity                  | $n=1$   | $n=4$   |
|--------------------------|---------|---------|
| $R_P/R_e$               | 2.1129  | 1.7140  |
| $f(<2R_P)$              | 0.9933  | 0.8165  |
| $R_{50}/R_e$ (Petro)    | 0.9936  | 0.7127  |
| $R_{90}/R_e$ (Petro)    | 2.2735  | 2.3881  |
| $R_{50}/R_P$            | 0.4703  | 0.4158  |
| $R_{90}/R_P$            | 1.0760  | 1.3933  |
| $C\equiv R_{90}/R_{50}$ | 2.288   | 3.351   |
| $r_1/R_e$               | 1.192   | 2.289   |
| $r_K(2.5)/R_e$          | 2.98    | 5.72    |
| $f(<r_K)$               | 0.960   | 0.904   |

</div>

Description: Summary of key ratios for $n=1$ (exponential) and $n=4$ (de Vaucouleurs) profiles. Shows how Petrosian, half-light, Kron, and concentration measures scale relative to $R_e$ and each other.

⸻

1) Rules of thumb
- Disks (low $n$): $R_P \sim 2R_e$; $R_{50}\approx R_e$; Kron(2.5) $\sim 3R_e$.
- Bulges (high $n$): $R_P \sim 1.7R_e$; $R_{50}<R_e$; Kron(2.5) $\sim 5$–$6R_e$.
- Seeing inflates all observed sizes; PSF-convolved model fits yield the cleanest intrinsic $R_e$ when the source is well resolved.
- Near the PSF limit, report a resolution metric (e.g., $R_{e,\rm obs}/{\rm FWHM}_{\rm PSF}$) and flag marginal cases.

⸻

7) Notes on survey choices
- Petrosian (SDSS/PS1): empirical, profile-driven apertures; robust to surface-brightness dimming and moderate seeing; can miss flux for high-$n$ halos.
- Sérsic $R_e$ (KiDS/Euclid): physically interpretable half-light sizes; depends on correct PSF modeling and the model family (single Sérsic vs multi-component).

⸻

8) References (indicative)

Petrosian (1976); SDSS technical descriptions (e.g., Strauss et al. 2002); Sérsic profile theory (e.g., Ciotti 1991) and $b_n$ approximations; fitting/systematics studies (e.g., Häußler et al. 2007; Biernaux et al. 2016); KiDS/Euclid morphology papers and challenges.