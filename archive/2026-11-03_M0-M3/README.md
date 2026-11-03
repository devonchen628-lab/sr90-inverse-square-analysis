# Sr-90/Y-90 inverse-square analysis

Imperial College London Year 2 radioactivity laboratory project. This analysis tests how finite geometry, counting losses and effective attenuation describe deviations from inverse-square behaviour in a sealed Sr-90/Y-90 source distance scan. Intended repository visibility: **PRIVATE**.

A circular emitting source (radius 4.875 mm) faces a 10 × 10 mm² silicon detector on a rail. The active-plane separation is the measured rail distance plus one common offset, 18.0 ± 1.5 mm. Source emission and detector efficiency are absorbed into a fitted normalization; this is not an absolute activity calibration.

## Data and physical models

`data/raw/Task16_2.csv` is the byte-preserved original 32-position dataset (1–85 cm), including counts `n`, exposure `t` in seconds, rail distance `d` and reading standard uncertainty `d_unc` in metres. Derived spreadsheet columns are retained but not fitted. `deadtime.csv` is the original 500-bin interval histogram. Background is a separate observation of 8 counts in 60 s.

For latent rail distance x, offset δ and source rate S:

- **M0:** S = C/x²; no offset correction, background or dead time.
- **M1:** S = C A/[4π(x+δ)²], with background and dead time.
- **M2:** S = C ε(x+δ), with background and dead time.
- **M3:** S = C ε(x+δ) exp(−μx), with background and dead time.

ε averages the analytic rectangular-detector solid angle over a uniform circular source, using 24 radial Gauss–Legendre × 96 angular samples. M1–M3 predict R = (S+B)/[1+τ(S+B)], with B = β/(1−τβ). One observed background rate β is shared between the scan and its Poisson control. τ is fixed at 1.97884 µs, with separate sensitivity fits at ±0.02126 µs.

## Statistics and findings

Independent Poisson scan counts are fitted jointly with Gaussian constraints on each latent distance and one Gaussian constraint on the shared offset. Offset uncertainty is correlated across all observations, not added independently to each point. All stated distance uncertainties are interpreted as standard deviations.

| Model | Scan counting deviance | Total deviance |
|---|---:|---:|
| M0 | 765.3310 | 2894.5347 |
| M1 | 867.1804 | 1052.8962 |
| M2 | 851.9419 | 1022.7902 |
| M3 | 190.9123 | 216.6054 |

Total deviance includes the scan counting term, background-control deviance where applicable, distance penalties and the common-offset penalty. These are distinct statistics. M3 gives μ = 0.440778 m⁻¹ (formal standard error 0.015888 m⁻¹), a one-unit profile interval of approximately 0.42503–0.45655 m⁻¹, and fitted offset 19.0446 mm. Its nominal total deviance per degree of freedom is 7.22 (30 degrees), indicating substantial unexplained variation.

At fixed offsets 16.5, 18.0 and 19.5 mm, M3 μ is 0.370450, 0.411925 and 0.453350 m⁻¹ respectively. These are conditional sensitivity scenarios, not independent measurements.

## Reproduce

Use Python 3.11+ (verified with Anaconda Python and the pinned packages):

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

The entry point works from any current directory, regenerates every final analysis figure and numerical table, checks total deviances and offset sensitivity against the original results, and writes `results/tables/verification.json`. Typical runtime is a few minutes. No external data, Geant4 or network access is required. For the historical report, install a TeX distribution and run `cd report && pdflatex main.tex && pdflatex main.tex`. Its two supplied figures are archived assets, not outputs of the current M0–M3 pipeline.

## Limitations

μ describes effective distance dependence, not a calibrated air absorption coefficient. The non-paralysable response is a working approximation; a dead-time process is not strictly Poisson. Background stationarity, alignment, aperture, energy response and additional repeat variation are not established. Formal intervals are conditional on fixed τ and incomplete error models. Nominal p-values and information criteria are heuristic with many latent distances and severe lack of fit.

The available histogram fit gives 1.998734 ± 0.020807 µs, rather than the prescribed scan constant 1.97884 ± 0.02126 µs. The latter is deliberately preserved; its exact derivation was not found locally. The histogram time axis (31.25 µs) and last-bin overflow interpretation are assumptions not encoded in the CSV. See `docs/methodology.md`.

The preserved two-column LaTeX report predates the final analysis and describes different geometry and transport models. It is a historical report, not a report of the M0–M3 results above. No matching final two-page M0–M3 LaTeX report was found.

## Layout

```text
data/raw/          immutable original scan and interval data
src/               geometry, models, comparison and dead-time analysis
results/figures/   regenerated model/residual and histogram figures
results/tables/    numerical results, sensitivity and verification
report/            preserved historical LaTeX source, class and images
docs/              methods, cleanup record, provenance and reference sensitivity
reproduce.py       complete analysis entry point and numerical verification
```
