# M1/M2/M3 comparison: methods and limitations

Run from the parent `air_attenuation_comparison` directory using Python with
NumPy, pandas, SciPy and Matplotlib:

```bash
python analysis_m3.py results_m3/data_used.csv --output-dir results_m3_reproduced
```

`analysis.py` must be beside `analysis_m3.py`: it supplies the unchanged,
verified finite-geometry integration and input checks. The original M1/M2
analysis and its fixed-offset figures remain in the parent folder. All 32 rows
are fitted in every model. The CSV is a snapshot of `Task16_2.csv`; its source
hash and geometry validation results are in `metadata.json`.

The data are modelled as Poisson counts with means t_i R(d_i). Measured rail
distances have individual Gaussian constraints of `d_unc`; the offset has one
shared Gaussian constraint of 18.0 ± 1.5 mm. Both ± values are assumed standard
deviations. All model parameters, including C and μ in M3, are refitted.

For M3, the eight background counts in 60 s are treated as a separate Poisson
measurement of the observed background rate β. β is one positive, common
profiled nuisance parameter, constrained by that control measurement and by
the scan data. The pre-dead-time rate is B=β/(1−τβ). M1, M2 and the dead-time
only diagnostic do not use background in their scan predictions, but the same
independent control is included at its saturated likelihood; hence the same
data and control are used for model comparison. Its nuisance parameter and
control observation cancel in the approximate degrees-of-freedom comparison.
The background measurement is assumed relevant to the scan configuration.

The nominal τ is held fixed at 1.97884 μs; its 1.95758 and 2.00010 μs
alternatives are refitted for each offset scenario. The two offset-and-tau
scenarios also include a constrained-offset fit at each τ. Formal μ standard
errors and profile intervals account for the shared offset and, for M3,
background control when fitted, but are conditional on τ. Tau sensitivity is
reported separately; it cannot represent model uncertainty. In particular, the
histogram cutoff has not been established as a non-paralyzable dead time.

Total deviance is Poisson count deviance plus per-point distance constraints,
the common offset constraint and, where applicable, the Poisson background
control deviance. Degrees of freedom, p-values and AIC/BIC differences are
nominal approximations; the many latent distances and clear lack of fit make
their usual calibration unreliable. The background-only and dead-only rows
are exploratory refits, not new unconstrained corrections. They reveal whether
each prescribed correction increases or decreases the objective, but objective
changes cannot be uniquely attributed by simple subtraction because parameters
are refitted and effects interact.

`effect_sizes.csv` isolates the direct change to the predicted rate with C, μ,
offset and the profiled M3 background held fixed. It compares S+B, S/(1+τS),
and (S+B)/(1+τ(S+B)) to S at the 32 observed distances. These are not fit
improvements, and maxima can change if other distances or parameters are used.

The plot shows the common-offset-constrained fits, and the marginal vertical
error bars include Poisson counting and propagated distance uncertainty. The
same reading uncertainty contributes to the horizontal bars, so these are
correlated. Residuals use local count and independent-distance uncertainties;
they omit the shared offset, background, tau, fitted-parameter covariance and
unmodelled errors. M2 and M3 curves may nearly overlap.

The imposed Poisson count likelihood is a working approximation: a genuine
non-paralyzable dead-time process is not exactly Poisson and can have reduced
count variance. Repeatability, background drift and other systematic effects
are also absent. Since reduced deviance is much larger than one, none of the
formal intervals should be presented as the complete measurement uncertainty.
