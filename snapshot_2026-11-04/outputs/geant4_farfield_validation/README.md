# Far-field experimental and Geant4 attenuation comparison

## Data and reproducibility

Run from any working directory on the same host:

```bash
MPLCONFIGDIR=/tmp/codex-mpl /opt/anaconda3/bin/python \
  /Users/chenqiyue/Documents/Codex/2026-10-29/r/outputs/geant4_farfield_validation/farfield_analysis.py
```

Requires NumPy, pandas, SciPy and Matplotlib. The script reads `Task16_2.csv`
from the Radlab folder and imports the existing verified geometry and staged
likelihood modules from `outputs/air_attenuation_comparison`. It reads existing
18 mm batches and `plan.json` in `outputs/geant4_aperture975_offset18`; it does
not modify or rerun them. The input SHA-256, batch seeds, generated histories,
bootstrap seed and numerical results are recorded in `summary.json`. The
script validates 32 shared distances, threshold, air/vacuum offset, source and
aperture metadata where present, positive probabilities and disjoint seeds.

The Monte Carlo `probability` is a weighted estimator of detected events per
original source history, not necessarily an unweighted count divided by a
generated count: the air calculation uses importance sampling and replay; the
vacuum calculation uses a source-exit library and prefilter. Scoring uses the
same assumed 10 keV energy-deposition threshold. There are eight independent air
libraries of 400,000 source histories and twelve independent vacuum libraries of
1,000,000 histories, each reused across distances. Whole source-library vectors
are resampled 20,000 times; air and vacuum batches use disjoint seeds and are
resampled independently. Do not count the 32 uses of a library as 32 new
independent simulations.

## Geometry and preselected region

For source radius 4.875 mm, square 10 mm detector and nominal 18 mm offset,
24-point radial Gauss–Legendre and 96-point angular integration gives a 1%
finite-geometry deviation at physical separation 64.982 mm, or rail distance
46.982 mm. A dense monotonicity check verifies the stated threshold holds over
the rest of the measured domain. The first available point is rail 50 mm
(nominal physical separation 68 mm), leaving 30/32 observations. The 0.5% and
2% tolerances give physical thresholds 92.222 and 45.624 mm, and first
available rail positions 75 and 50 mm, respectively. This criterion addresses
only the ideal analytical geometry and does not establish an air-free baseline
or validate the simulated detector response.

The task specifies 30 cm as the preferred reference where available. It is
available and used here (physical separation 31.8 cm); this differs from the
previous 35 cm display. Because the offset is common, `s-s_ref=d-d_ref`.
Normalisation introduces correlations and makes the reference error exactly
zero by construction, not by independent measurement.

## Experimental likelihood

The M3 likelihood uses original counts and times, with one Gaussian latent
rail-position shift per observation, a shared Gaussian offset constraint of
18.0 ± 1.5 mm, a single background nuisance constrained by the independent
8-count/60-s control, and fixed effective non-paralyzable tau = 1.97884 μs.
The 0.5 mm distance and 1.5 mm offset are treated as standard deviations.
The pre-dead-time B is obtained from the observed background nuisance β as
`B=β/(1−tau*β)`. This assumes the control applies to the scan. A genuine
non-paralyzable counter is not exactly Poisson; the likelihood is a working
approximation. The supplied timing cutoff has not been independently validated
as this model's dead time.

The verified original code uses `C_d exp(−mu*d)`. The requested physical
distance form is mathematically identical with
`C_physical=C_d exp(mu*offset)` and `s=d+offset`; both C conventions are saved,
and numerical equality is checked. No old model files were changed. The
experimental μ profile interval is the one-unit total-deviance profile with
background, individual distances and common offset profiled, conditional on
tau and model correctness. The rather poor fit and unmodelled repeatability
mean this is not a complete uncertainty budget.

The primary constrained-offset fit gives mu = 0.4718 m⁻¹, profile interval
0.4531–0.4905 m⁻¹, scan Poisson deviance 177.94 and augmented total deviance
204.99 for 28 nominal df. The corresponding full-range M3 result is 0.4408
m⁻¹. Fixing offset at 16.5, 18 and 19.5 mm instead gives 0.3761, 0.4113 and
0.4463 m⁻¹. These are sensitivity refits on the same 30 observations; fixed
offset does not include a Gaussian offset penalty. In the primary joint fit
the offset profiles to 20.60 mm, beyond the nominal +1σ but constrained by the
measurement. Do not interpret that fitted offset as an independently measured
value.

## MC slope, vacuum check, and interpretation

The vacuum probability divided by analytical geometry is normalised at 30 cm.
The ratio varies markedly at shorter distances even within the selected
geometry-only region. A constant ratio cannot be treated as validated. Such a
variation can bias interpretation of `P_air/P_vac` as the only multiplicative
correction to the analytical model.

The MC exponential summary fits `log(T/Tref)=-mu(d-dref)` with a covariance
estimated from whole-batch bootstrap vectors. Since the sample covariance is
rank-limited, a fixed 20% diagonal shrinkage is applied. The same covariance
weights are used for each bootstrap slope. This regularisation and the small
number of source libraries mean neither the resulting shape statistic nor its
centred-bootstrap tail fraction should be treated as a fully calibrated p-value.
The vacuum constant diagnostic uses the same approach. No fitted physical
parameter is introduced to force consistency.

The anchored MC slope is 0.5421 m⁻¹; the batch-bootstrap median is 0.5212 m⁻¹
(16–84%: 0.3153–0.7373 m⁻¹). The central MC-minus-experimental slope is
0.0703 m⁻¹. Conditional independent MC/experimental resampling gives a
difference median of 0.0495 m⁻¹ and 16–84% interval −0.1567 to 0.2659 m⁻¹;
this is not a full systematic-uncertainty statement or a calibrated hypothesis
test. Experimental μ draws use the local formal standard error; MC uses whole
batches. Fit residuals and uncertainties are recorded separately.

At rail 85 cm versus the common 30 cm reference, the simulated ratio is
0.537 (MC 16–84%: 0.446–0.698), while experimental M3 predicts 0.771.
The experimental model is outside the pointwise MC 16–84% band at 27.5 and
85 cm, but those intervals are not simultaneous confidence bands and cannot
alone establish significance. The slopes and the full curve therefore give
different impressions.

MC endpoint effective slopes differ: 30–55 cm gives 0.508 m⁻¹ (batch 16–84%:
−0.089–0.993), while 55–85 cm gives 1.651 m⁻¹ (1.111–2.032). The paired
bootstrap of their difference has median 1.144 m⁻¹, interval 0.538–1.704.
This indicates curvature in this assumed simulation; small-library and model
uncertainties limit the inference. A single exponential coefficient is an
interval- and weighting-dependent descriptive summary, not an exact transport
law. Alternative predeclared cutoffs and their recalculated reference values
are in `cutoff_sensitivity.csv`.

Important model discrepancies remain: the Geant4 emitting disk is assumed
22 mm across but analytical source is 9.75 mm across; the Geant4 exit aperture
is 9.75 mm. Geant4 also assumes a 10 keV threshold, encapsulation/foil/acrylic,
material and angular/energy distributions that are not independently fully
validated. Its MC bands include sampling uncertainty only. The experimental
exponential is phenomenological and does not isolate air. The data support
decreasing response qualitatively, but do not establish quantitative validation
or unique air causation.

`figure1_farfield_selection`, `figure2_transmission`, and
`figure3_log_diagnostic` are supplied in PDF and PNG; use the main comparison
selectively in a two-page report. Other CSVs contain the fitted values, selected
points, individual residuals, offset checks, interval slopes and MC diagnostics.
The optional simulation-forward model is not claimed or included.
