# Finite geometry and effective attenuation

## Reproduce

The included `data_used.csv` is a snapshot of all 32 measurement rows in
`Task16_2.csv`; `metadata.json` records the source path and SHA-256 digest.
Requires Python 3 with NumPy, pandas, SciPy and Matplotlib. For example:

```bash
/opt/anaconda3/bin/python analysis.py data_used.csv --output-dir results
```

Run from this directory. The script generates the comparison CSV, the two-panel
PDF/PNG, a pointwise residual CSV and provenance/geometry-check metadata.

## Statistical model

For measured rail position `d_i`, the latent true rail position `x_i` is
constrained by `(x_i-d_i)/d_unc_i ~ N(0,1)`. Counts follow independent Poisson
distributions with expectation `t_i R(x_i)`. Offset is either fixed at 16.5,
18.0 or 19.5 mm or treated as one shared nuisance parameter with an independent
normal constraint of 18.0 ± 1.5 mm. These uncertainties are interpreted as
standard deviations; this must be revisited if the reported ± values instead
represent hard bounds. Every fit uses the same observations.

The figure shows fixed-offset 18 mm fits. In panel (a), y-error bars propagate
both Poisson counting and independent distance-reading errors into `R d²`, while
x-error bars show the same reading error: x and y are correlated. Fitting uses
the original counts and positions, not those transformed error bars. Panel (b)
uses residuals divided by the local counting plus independent-distance error,
evaluated at the observed distance. Residuals are not independent after fitting;
common-offset uncertainty, fitted parameter covariance, and unmodelled errors
are not included in their denominators. They are diagnostic, not significance
tests for individual points.

`augmented_deviance` is the Poisson likelihood-ratio deviance plus the squared
Gaussian constraint residuals. It is asymptotically chi-square distributed if
the model and error assumptions hold. The quoted p-values and reduced deviance
are therefore nominal approximations. The information-criterion differences
are approximate diagnostics, not independently validated evidence: every point
has a latent distance, and model misspecification and overdispersion undermine
the usual calibration. ΔAIC and ΔBIC are relative to M1 for the same treatment
of offset. They do not establish a unique cause.

Uncertainties reported for C and μ from the local curvature are formal 1σ
estimates under the stated likelihood. Profile intervals in the table use a
one-unit deviance increase and are provided for the nominal fixed-offset and
shared-offset fits. Fixed-offset results are conditional on that offset. The
shared-offset profile includes its correlated measurement uncertainty, but not
unknown setup errors or model inadequacy. Since the fits are poor, none of
these intervals represents a complete experimental uncertainty budget.

There is no independently verified background, no dead-time correction, and
no Geant4 input in this analysis. C incorporates emission and detection
efficiencies and is not a calibrated source activity. The fitted μ describes
effective distance dependence, not a direct measurement of the air attenuation
coefficient. The exponential multiplies exp(-μ d) as specified, while the
geometric separation is d + offset.

## Checks

The script independently compares source integration with a finer quadrature,
checks its analytic distance derivative against finite differences, and compares
the solid-angle integral to a separate source-and-detector area integral. It
aborts if any relative disagreement exceeds 2 × 10⁻⁶. A separate numerical
Jacobian check of the complete constrained fit was also performed.
