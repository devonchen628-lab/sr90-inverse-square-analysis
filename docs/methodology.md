# Statistical and implementation details

Geometry integrates over squared radius using Gauss–Legendre nodes and uniform azimuth; the rectangular solid angle is analytic. Validation compares doubled quadrature resolution, an independent detector-plane surface integral, and finite-difference derivatives, requiring relative error below 2e-6.

The optimizer minimizes squared signed Poisson deviance residuals plus Gaussian constraint residuals. The counting deviance is 2 Σ[λ−n+n log(n/λ)]. Background contributes the same expression for the independent 8-count control. Total adds Σz_i² and ((δ−0.018)/0.0015)². M0 has no physical dependence on δ; its shared constraint is consequently zero. A saturated independent background control would also add zero for M0.

M0 allows large latent-distance shifts to avoid the earlier arbitrary 12σ clipping. Other models retain the established ±12σ latent-distance and ±8σ offset bounds. M0's distance penalty is 2129.20: model inadequacy drives implausibly large adjustments, and its counting-only deviance must not be substituted for the total. M1 and M2 pull the common offset to about 11.69 and 11.20 mm; these are model-dependent fits, not hardware measurements.

The approximate degrees of freedom are 31 for M0–M2 and 30 for M3. Counting, Gaussian constraints, latent nuisance parameters and shared background must be considered together. AIC/BIC in the generated comparison are heuristic, not calibrated evidence of a physical cause.

Plot error bars propagate Poisson counting and independent reading errors into R d²; horizontal and vertical errors therefore correlate. Diagnostic residuals use predicted rate variance plus the local distance-slope term, evaluated at observed rail distances. They omit common-offset and fitted-parameter covariance, background and dead-time uncertainty and unmodelled systematics; they are not independent significance tests.

The M3 profile interval refits background, normalization, offset and latent distances at fixed μ and uses ΔD=1. It is conditional on τ. Offset/dead-time sensitivity covers fixed offsets 16.5/18/19.5 mm and a shared constrained offset, each at τ=1.95758/1.97884/2.00010 µs.

The retained deadtime.py fits Gaussian recovery delay convolved with exponential waiting times, including the overflow bin. It also compares a sharp cutoff and free finite-range amplitude. This descriptive recovery fit does not establish a hard non-paralysable detector dead time. The 1.97884 µs scan constant is a preserved supplied value, not silently replaced by the approximately 1.99873 µs local fit.

The report's assumed 21 mm offset, illustrative 22 mm source, Geant4 transport and older statistics are historical and differ from the final analysis. Its source and referenced images are preserved unedited. Current analysis contains no Geant4 dependencies. Rebuilding its archived plots would require the earlier, separate research pipeline, deliberately outside the final model scope.
