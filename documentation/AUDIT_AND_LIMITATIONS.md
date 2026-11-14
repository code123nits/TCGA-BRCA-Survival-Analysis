# Completion and integrity audit

**Completed in this handoff:** consolidated manuscript, 3 original result figures, reproducibility script copied from prior archived package, syntax compilation, archived predictions and model-metric JSON, supplementary feature/CV tables, and ZIP packaging.

**Not claimed:** fresh model fitting in this session, full raw-to-result rerun, a truly untouched held-out test cohort, independent external validation, fully verified literature metadata, or clinically deployable predictions.

**Source conflicts resolved:** older Stage 14E results (C-indices ~0.6813/0.6835) describe a different anchored/grouped-stage model specification and MUST NOT be combined with frozen 13-feature fixed-alpha results (0.6712/0.6175). Separate training-CV-selected alpha results (clinical ~0.5832, combined ~0.5854) are not the fixed-alpha primary model results. Older calibration/stability plots are excluded.

**Interpretation:** intervals include zero; no established superiority, no claim genes are biologically irrelevant. SCAN-B readiness is not external validation.


---

## Update 2026-11-13: checks added in the package revision

**Verified from archived outputs (no model refitting)**
- Harrell C-indices recomputed from `results/test_predictions.csv` equal `results/model_metrics.json` exactly
  (clinical 0.6712148, combined 0.6175176, combined sensitivity 0.5853873).
- `code/bootstrap_cindex.py` with seed 42 and 2,000 resamples reproduces every manuscript interval:
  clinical 0.5243-0.8040; combined (alpha 0.03) 0.4848-0.7420; combined (alpha 0.05) 0.4380-0.7363;
  differences -0.2126 to 0.0972 and -0.2379 to 0.0481.

**Tested only for plumbing (not on real data)**
- `code/frozen_refit.py` and `code/repeated_cv.py` were run once against synthetic data with a throwaway
  stand-in for scikit-survival, to check file handling, output names and fold logic. No real model was fitted,
  so these scripts have **not** been shown to reproduce the archived models. Run `frozen_refit.py`, then
  `compare_to_archive.py`, before claiming reproduction.

**Penalty sensitivity (exploratory, from the 3-fold CV files)**
| Model | Frozen alpha | Test C-index | CV-selected alpha | Test C-index |
|---|---|---|---|---|
| Clinical | 0.001 | 0.6712 | 0.01 | 0.5832 |
| Clinical + RNA | 0.05 (sensitivity model) | 0.5854 | 0.05 | 0.5854 |

With CV-selected penalties the two models are essentially tied on the held-out patients, so the ordering in the
primary analysis depends on the penalty chosen and on only 30 test events. Caveats on the 3-fold run: three
folds only; the best clinical alpha (0.01) is the smallest value in that grid; the grid did not include the frozen
clinical alpha (0.001); alpha = 0.2 zeroed all coefficients. `repeated_cv.py` addresses these but has not yet been
run on the real data.

**Still outstanding**
Raw-to-matrix build script and GDC manifest; alpha-selection provenance for the frozen models; a harmonized
external cohort; calibration and time-dependent evaluation not touching the reused test set; proportional-hazards
check; full reference verification; author declarations.
