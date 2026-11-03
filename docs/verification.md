# Verification record — 3 November 2026

The full pipeline ran with /private/tmp as working directory and newly created temporary fit-output directories. All 32 rows were used. Total deviance checks passed at rtol=1e-7, atol=1e-4; μ matched within 1e-6 m⁻¹. All 12 fixed-offset M0–M3 fits matched the original sensitivity table within rtol=2e-6, atol=1e-6. Twelve additional M3 offset/dead-time combinations were regenerated. Geometry integration, independent surface integral and derivative checks passed the 2e-6 threshold.

The final M0–M3 curve/residual PDF and PNG and the dead-time diagnostic PNG were generated. Source raw files were not modified. Pinned package versions were those used for this run; dependency installation into a new virtual environment was not separately tested.

The preserved LaTeX report compiled twice in a scratch copy with pdflatex to two pages. It had no missing image dependencies; its existing class-name and underfull-box warnings remain. Its scientific content predates the final models.

The histogram soft delay is about 1.998734 µs, versus 1.97884 µs prescribed for the scan. This discrepancy is documented rather than fitted away. M3 remains a poor fit despite its improvement; its total deviance is not its counting-only deviance.

Machine-readable values: results/tables/verification.json. Exact component breakdown: M0_M3_comparison.csv. Offset/dead-time scenarios: offset_deadtime_sensitivity.csv.

Source-directory integrity note: all operations performed by this consolidation on the original directory were read-only. A size/mtime audit observed unrelated scan_dense32 and scan_precision64 simulation CSV/log files changing during inspection, consistent with concurrent external simulation activity. The selected copied inputs and their provenance hashes are fixed; reproducibility does not depend on those changing simulations.
