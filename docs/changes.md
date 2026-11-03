# Files changed in the three-model revision

- src/models.py: reuse the original likelihood/analytic Jacobian; define core M0/M1/M2 plus P1/P2 sensitivity models; unify background/dead-time response and distance bounds; return conditional background/offset uncertainties and μ–δ correlation; expose rate derivative.
- src/geometry.py: preserve geometry and input-validation algorithms; remove unused imports.
- src/comparison.py: regenerate component deviance and parameter tables, offset/tau/point-source checks, latent distances, nominal-position residuals, prediction grids and numerical verification.
- src/plotting.py (new): publication-size three-panel vector PDF and 300 dpi PNG with accessible styles and conditional offset-sensitivity errors.
- reproduce.py: run the current analysis, figure and histogram diagnostic; validate raw hashes and image export dimensions.
- results/tables/: current model_comparison.csv, residuals.csv, offset_sensitivity.csv, offset_deadtime_sensitivity.csv, point_source_sensitivity.csv, latent_distances.csv, model_predictions.csv and verification.json; the retained deadtime_fit.json is regenerated too.
- results/figures/: final_three_panel.pdf and final_three_panel.png; the retained deadtime_recovery.png is regenerated too.
- README.md and docs/methodology.md: current model/response definitions, comparable statistics and limitations.
- docs/results.md, docs/figure_usage.tex, docs/verification.md: quantitative interpretation, report insertion/caption and verification evidence.
- docs/cleanup.md and this file: preservation and revision record.
- report/README.md: explicitly mark the preserved report as historical. Historical .tex/.cls/images are unchanged.
- archive/2026-11-03_M0-M3/: preserve previous scripts, all numerical/figure outputs and documentation before replacement of the active hierarchy.

Raw data, experimental constants, requirements.txt and historical report assets are unchanged. No writes were made to the original 2026-10-29/r project.
