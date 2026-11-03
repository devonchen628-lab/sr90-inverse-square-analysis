# Verification of the current M0–M2 revision

Both the staged project and an independent input-only copy were run with /private/tmp as the current directory. The copy contained source, immutable data, documentation and historical references, but no current results directory. Both complete reproduce.py runs passed and generated all current tables and figures.

## Numerical checks

- Original raw CSV SHA-256 values match provenance and remain unchanged after each run.
- All 32 scan rows are fitted, with a single shared background parameter and fixed τ for every core model. M0 is newly fitted with that response and has no dummy offset parameter.
- Geometry refinement, independent detector-plane integration and finite-difference derivative checks all passed the existing 2e-6 threshold. Relative discrepancies were 0, approximately 8.9e-16 and 2.1e-10 respectively at the test separations.
- Individual component sums agree with total optimizer deviance. Nominal likelihood-observation and parameter counts independently give degrees of freedom 31/31/30.
- Alternate starts (1.5× normalization, initial μ=0.8 m⁻¹ where fitted) give total-deviance differences below 1e-8. Analytic residual Jacobians agree with independent three-point numerical derivatives, with maximum relative column errors approximately 1.0e-8 (M0), 1.6e-5 (M1), 9.9e-7 (M2).
- No core fit has an active optimization bound. All latent distances share the same positivity domain; the prior arbitrary model-dependent sigma clipping is removed.
- Current M1/M2 totals and M2 μ agree with historical finite-geometry/attenuation results; these are the unchanged physical predictions under new labels. Current M0 is explicitly checked to differ from the old no-response M0. Its total is 2785.1104823, not 2894.5347284.
- All three fixed-offset M2 scenarios match the original sensitivity calculations; their full component deviances and conditional errors are regenerated. The point-source checks and all twelve offset/dead-time scenarios are also regenerated.
- Every CSV from the independent fresh copy matches the staged output (rtol=1e-9, atol=1e-8 for numeric columns; exact nonnumeric content). Histogram and numerical-verification JSON results match exactly.

The current tables checked were: latent_distances.csv, model_comparison.csv, model_predictions.csv, offset_deadtime_sensitivity.csv, offset_sensitivity.csv, point_source_sensitivity.csv, residuals.csv. Package versions remain those pinned in requirements.txt (Python 3.13.5 for these runs). A new dependency installation was not separately tested.

## Figure and preservation checks

The final vector PDF was rendered with Poppler and visually inspected. It is one page, 178 × 130 mm; it embeds TrueType font data and contains no raster image objects. The PNG dimension and 300 dpi metadata checks pass. A separate two-column A4 LaTeX print proof using docs/figure_usage.tex compiled successfully and was visually reviewed for text, units, panel labels, caption and clipping. Font sizes are intended for the full 178 mm width, not one column. The proof is an internal validation artifact rather than a replacement historical report.

Twenty-five historical code/result/document files were checked byte-for-byte against the previous Git commit after archival. Original raw CSV files, historical report source/class and both report figures were also checked unchanged. All prior result files remain available in archive/2026-11-03_M0-M3/.

## Scientific limits

Passing numerical checks verifies implementation and reproducibility, not the physical adequacy of the models. The M2 total/dof remains 7.22, shared-offset and μ estimates correlate, and diagnostic residuals retain structure. The histogram soft-delay discrepancy and the unvalidated fixed non-paralysable response are preserved and documented. No artificial error rescaling or geometry adjustment was introduced to improve fit.
