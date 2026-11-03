# Cleanup and source selection

The original project was read only. A separate allowlist staging directory was created; no original files were deleted or edited. Content inspection, source imports, figure references and timestamps were used together.

Retained: byte-identical original Radlab Task16_2.csv (32 rows; hash matches original metadata), deadtime.csv from the archived research inputs, final analysis_staged.py model implementation and final_comparison_table.py logic. Geometry and input-validation functions were extracted from analysis.py because these are active dependencies; its older Fit class and main routine were omitted. Renamed modules use portable imports. M0 was added to the final plot using the unchanged model; all M0–M3 statistics remain unaltered. The misleading inherited background comment was corrected.

The final comparison script (Nov 3 15:30) imports analysis_staged.py (15:31), and its reference table matches the requested results. analysis.py (14:48), analysis_m3.py (15:17) and analysis_merged.py (15:27) are earlier model variants; their full routines and copied results folders are excluded because their model/background definitions differ. Earlier documentation also describes outdated model choices and was replaced with implementation-based documentation.

Excluded: root analysis outputs, results_m3/results_merged/results_shared_background/results_staged duplicates, redundant data_used.csv snapshots, intermediate CSV/PNG/PDF copies, cache files, .DS_Store, TeX build products, ZIP archives, debugging work, repository reviews, unrelated lab tasks and example reports. Final outputs are regenerated once under results/. A single original sensitivity CSV remains under docs/ as an explicit regression reference, rather than an unexplained duplicate result.

The most recent local two-column report source found is outputs/sr90_research/tex/main.tex; it references exactly report_data.png and report_geometry.png. Its custom class and both images were retained. Older radioactivity_report/theme_report/reframed versions were excluded as superseded. No report references point to excluded assets. Its archived PDF and duplicate PDFs were omitted; source builds independently.

src/deadtime.py is adapted from work/radioactivity/fit_new_deadtime.py with portable input/output paths; the obsolete sr90_research/02_timing.py duplicate is excluded. Histogram values are not used to overwrite the established scan constant.

Provenance records original source hashes and timestamps for copied assets. Raw data are never rewritten. No secrets, unrelated personal documents or external credentials are included. Original inventory and build logs remain outside the repository in the task work directory.
