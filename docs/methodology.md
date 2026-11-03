# Current methodology: common-response M0–M2

## Physical hierarchy and fixed inputs

For latent rail distance $x_i$, use $S_0=K/x_i^2$, $S_1=C G(x_i;\delta)$ and $S_2=C G(x_i;\delta)e^{-\mu x_i}$. Attenuation uses rail distance as in the established analysis. For fixed μ and free normalization, using $e^{-\mu(x+\delta)}$ would simply reparameterize C, but does not establish a physical propagation path.

G is dimensionless: the average of detector solid angle divided by $4\pi$ over a uniform emitting disc of radius 4.875 mm. The square detector is the established 10 × 10 mm² rectangle. Gauss–Legendre integration uses 24 nodes in squared radius and 96 uniformly spaced angles. The existing solid-angle formula and derivative are unchanged. Verification compares doubled quadrature, an independent detector-plane area integral and finite-difference derivatives; the threshold is 2e-6 relative discrepancy.

M0→M1 changes both the distance reference and finite geometry. It is not a controlled test of finite source size alone. P1/P2 preserve the previous approximation $A/[4\pi(x+\delta)^2]$, without/with attenuation, using the same response and nuisance treatment as the corresponding finite models. This is the old point-source/small-detector approximation, not an exact point-to-square solid angle.

## Common measurement response and joint likelihood

Every core and sensitivity model uses

$$R_i=\frac{S_i+B}{1+\tau(S_i+B)},\quad B=\frac{\beta}{1-\tau\beta},\quad \tau=1.97884\,\mu\mathrm{s}.$$

The independent background observation is $n_b=8$ in $t_b=60$ s, so its expected counts are $t_b\beta$. This parameterization distinguishes the observed background β from the physical pre-dead-time B. One positive β is shared between the scan and control. It is fitted rather than subtracting a noisy control from the scan and then imposing a Poisson model on corrected counts.

The scan means are $\lambda_i=t_iR_i$. Measured distances $d_i$ constrain $x_i=d_i+s_i z_i$ with independent standard-normal residuals $z_i$, where $s_i$ is the original `d_unc`. M1/M2 have one common $\delta\sim N(18.0,1.5^2)$ mm; M0 has no δ and no dummy offset parameter. This common parameter induces correlations across all scan predictions. Gaussian constraints interpret the supplied ± values as standard deviations, not hard limits.

The augmented deviance is

$$D_{\rm total}=D_{\rm scan}+D_{\rm bg}+\sum_i z_i^2+[(\delta-0.018)/0.0015]^2,$$

with the final term present only where δ is fitted. Each Poisson component uses $2\sum[\lambda-n+n\log(n/\lambda)]$, with the zero-count limit. These terms are exported separately; the sum is checked numerically against the optimizer objective. Fixed-offset scenarios omit the offset constraint rather than adding a penalty for a value that is not being estimated.

The same latent-distance domain is used for every model: positive x, with no arbitrary ±12σ clipping. Broad established bounds on log normalization, log β, nonnegative μ (0–20 m⁻¹), and shared offset pull (±8σ) are retained for numerical optimization. None of the core best fits touches a bound. No constants, source dimensions or uncertainties are tuned to improve the result. M0 has approximately a 30σ maximum distance pull, making its severe inadequacy explicit.

There are 32 scan counts, 32 distance observations, one background control and (where applicable) one offset observation. M0 fits 34 parameters (32 latent distances, K, β) against 65 observations. M1 fits 35 against 66; M2 fits 36 against 66. Nominal degrees of freedom are therefore 31, 31 and 30. Extra offset/background nuisance parameters are balanced by their independent constraints. These are nominal augmented-likelihood degrees of freedom; the counting-only term is not assigned this same calibration by itself.

## Uncertainties and model comparisons

Formal covariance is $(J^TJ)^{-1}$ (numerical pseudoinverse) for signed deviance residuals, without multiplying by reduced deviance. Normalization, background and offset uncertainties are transformed from the fitting coordinates. Each covariance profiles the remaining fitted nuisance parameters locally. The M2 μ profile interval refits normalization, β, δ and all latent distances at fixed μ until total deviance increases by one. All such uncertainties are conditional on fixed τ, correct likelihood/constraints and the assumed physical model. They do not include model inadequacy, source dimensions, alignment, threshold response or unknown systematics.

The exported improvement is the absolute and percentage reduction in **total** deviance from the preceding core model. M0/M1 are not a simple nested single-mechanism test. M1/M2 differ by nonnegative μ, with the null at its boundary, so even a conventional one-parameter likelihood-ratio significance requires care. No deviance reduction is used as proof of absorption or another causal mechanism. Nominal chi-square goodness-of-fit p-values document poor fit but are not calibrated under unmodelled variance and systematics. No AIC/BIC-based physical claim is made.

## Residual and plotting definitions

Panel (a) plots observed $R_i^{\rm data}=n_i/t_i$ times the measured $d_i^2$, with

$$s_{Rd^2}=\sqrt{(d_i^2\sqrt{n_i}/t_i)^2+(2R_i^{\rm data}d_i s_i)^2}.$$

Horizontal bars are $s_i$. The same distance reading contributes to x and y, so the displayed errors are correlated. Curves use fitted observed rates at nominal rail d, rather than replacing x-axis points by fitted latent distances. Fitting uses raw counts, not transformed y or its propagated errors.

Panel (b) reports

$$r_i=\frac{n_i/t_i-\widehat R(d_i)}{\sqrt{\widehat R(d_i)/t_i+[s_i\widehat R'(d_i)]^2}}.$$

It compares M1 and M2 using predicted counting variance and local independent reading error. The denominator omits fitted-parameter covariance, common-offset/background correlations, dead-time uncertainty and model inadequacy. Residuals share fitted parameters and are not independent significance tests. The ±2 lines are visual guides, not discovery thresholds. These diagnostic residuals also differ from the signed Poisson deviance residuals at fitted latent distances (exported separately).

Panel (c) refits M2 with fixed δ=16.5/18/19.5 mm and nominal τ. Error bars are local conditional statistical 1σ, with β, C and latent distances fitted, but **no** uncertainty in the fixed δ or τ. The three scenarios do not have an assigned probability distribution and their span is not a confidence interval. The primary shared-offset fit has a different, larger formal μ uncertainty because δ is allowed to vary and correlates with μ.

The final figure is 178 × 130 mm: a full-width `figure*` for a two-column A4 report. It has vector lines/text and embedded fonts in PDF, and a 300 dpi PNG. Labels use 7–8 pt at the intended width, with a grey/blue/vermillion colour-blind-friendly palette and distinct line/marker shapes. It should not be shrunk to a single 85 mm column.

## Additional sensitivities and unresolved assumptions

`offset_deadtime_sensitivity.csv` refits M2 at τ=1.95758/1.97884/2.00010 µs for each fixed offset and for the shared constrained offset. Varying τ this way is a sensitivity calculation, not a marginalized combined uncertainty. `point_source_sensitivity.csv` is separate from the core hierarchy.

The retained histogram analysis fits Gaussian recovery delay plus exponential waiting times, includes overflow and checks sharp-cutoff/free-amplitude alternatives. Its approximately 1.998734 µs soft delay is not identical to the prescribed 1.97884 µs scan response constant; no local derivation of the latter was found. The 31.25 µs display range and final-bin overflow interpretation are assumptions not encoded in the file. No recovery parameter is a calibrated hard non-paralysable dead time. A true counting-loss process need not generate Poisson counts. Background stationarity during the later scan is assumed, not independently established.

A broad beta energy spectrum, scattering, detector efficiency and threshold, source encapsulation/apertures and alignment can all produce departures from a single exponential. μ is phenomenological and cannot be interpreted as the true beta absorption coefficient of air. The poor fit and structured residuals rule out presenting the formal uncertainty as a complete experimental uncertainty budget. The original report and complete earlier M0–M3 results remain explicitly historical.
