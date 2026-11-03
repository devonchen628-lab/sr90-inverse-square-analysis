# Quantitative interpretation of the current three-model analysis

The core models use identical raw measurements and common measurement response. Original constants remain unchanged; the entire earlier analysis is archived.

| model | scan_poisson_deviance | background_control_deviance | distance_constraint_penalty | offset_constraint_penalty | total_deviance | dof | reduced_total_deviance |
|---|---|---|---|---|---|---|---|
| M0 | 772.662709 | 0.116967 | 2012.330807 | 0.000000 | 2785.110482 | 31.000000 | 89.842274 |
| M1 | 851.940718 | 8.223016 | 142.061700 | 20.564779 | 1022.790214 | 31.000000 | 32.993233 |
| M2 | 190.912286 | 0.408258 | 24.799884 | 0.484970 | 216.605398 | 30.000000 | 7.220180 |

These are distinct deviance components. M0→M1 reduces total D by 1762.3203 (63.2765%), but the scan counting term rises from 772.6627 to 851.9407. The distance penalty drops from 2012.3308 to 142.0617. Thus the improvement mainly relieves extreme fitted distance adjustments; it must not be described as a 63% improvement in counting-only agreement or assigned entirely to finite source size.

M1 adds the established finite disc/square acceptance and one shared distance-reference offset. Its fitted δ=11.1977 ± 0.2717 mm (formal local 1σ) is about 4.53 prior standard deviations below 18 mm, and its maximum individual distance pull is about 8.60σ. The reduced total deviance is still 32.99. Fitted geometry-dependent parameters must not be interpreted as an independent measurement of the hardware.

M2 adds exp(−μx). It reduces total D by 806.1848 (78.8221%) relative to M1, and by approximately 92.22% relative to the newly fitted M0. Its fitted normalization is C=(2.071569 ± 0.017983)×10⁶ s⁻¹, μ=0.440778 ± 0.015888 m⁻¹, δ=19.0446 ± 0.4017 mm, and observed background β=0.105437 ± 0.040171 s⁻¹. M0 has K=11.90854 ± 0.02893 m² s⁻¹ and β≈0.117854 ± 0.041964 s⁻¹; M1 has C=(1.650568 ± 0.005399)×10⁶ s⁻¹ and β≈0.039436 ± 0.019981 s⁻¹. These are conditional local errors at fixed τ; all fitted nuisance correlations are included locally, but model/systematic uncertainty is not. The primary M2 μ–δ correlation is approximately +0.702. The M2 profile interval ΔD=1 is 0.425031–0.456552 m⁻¹.

## Offset sensitivity

| Fixed δ (mm) | μ (m⁻¹) | Conditional local 1σ (m⁻¹) | Total D |
|---|---:|---:|---:|
| 16.5 | 0.370450 | 0.011317 | 256.5132 |
| 18.0 | 0.411925 | 0.011316 | 223.4553 |
| 19.5 | 0.453350 | 0.011318 | 216.8990 |

The scenario span is 0.082900 m⁻¹ over 3 mm, or about 0.02763 m⁻¹ per mm. This is about 7.3 times an individual conditional standard error and illustrates sensitivity to distance reference. It is not a confidence interval or a calibrated systematic-error distribution. A fixed 18 mm offset gives a different μ than the primary shared-offset fit, which selects 19.0446 mm. The constrained-offset τ sensitivity changes μ from about 0.440847 to 0.440708 m⁻¹ across the prescribed ±0.02126 µs range; this is much smaller than the selected offset-scenario span but does not bound dead-time model inadequacy.

## Is finite geometry independently established?

The preserved shifted point-source/small-detector check P1 gives D=1052.8962, versus 1022.7902 for finite geometry M1, a 30.1060 (2.8594%) decrease at the same nominal degrees of freedom. This is far smaller than the combined M0→M1 change, so that change must not be attributed to finite acceptance alone.

With attenuation allowed, P2 gives μ=0.448585 m⁻¹ and D=212.7110, slightly **lower** than the finite model's 216.6054; finite geometry shifts μ by approximately −0.00781 m⁻¹. The assumed physical dimensions motivate the finite model, but this dataset does not uniquely validate those dimensions or prefer the finite acceptance in every comparison. No radius/area was retuned to reverse that finding.

## Remaining residuals and evidential limits

M1 retains a broad positive middle-distance band and a strong terminal deficit: its nominal-position standardized residual at 85 cm is about −12.62. M2 reduces these large-scale trends, but keeps a positive band around 35–52.5 cm (for example +5.45 at 47.5 cm), and negative terminal departures (−3.76 at 80 cm, −4.98 at 85 cm). These labels refer only to the defined diagnostic residual scale, not independent sigma-level discoveries. The common fitted parameters correlate residuals, and their denominators omit shared systematic uncertainty.

M2's D/ν=216.6054/30=7.22 and nominal chi-square p≈3.7×10⁻³⁰ demonstrate severe lack of fit under the working likelihood, not an accurately calibrated probability under unknown systematics. A one-exponential distance dependence, simple uniform geometry and fixed response do not fully describe the observations. Energy-dependent transport, scattering/threshold response, misalignment, source structure, drift and extra repeat variation are plausible unresolved effects; these data and fits do not establish which causes the remaining pattern. Adding arbitrary flexible fit functions or inflating individual measurement errors to hide it is unwarranted.
