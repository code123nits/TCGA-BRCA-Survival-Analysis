## Public GitHub draft cleanup (November 2026)
- Removed bytecode cache, superseded manuscript, internal editing notes, and publication checklist.
- Replaced TCGA case IDs in the saved prediction table with sequential sample labels; predictions and outcomes are unchanged.
- Updated README and manuscript to refer to the current draft.
- No new model fitting or external validation was performed.

# Changelog

## 2026-11-13 (package revision)
**Added**
- `code/cindex_utils.py`: Harrell's C-index following scikit-survival conventions. Reproduces the archived test
  C-indices exactly (abs. difference 0.0).
- `code/bootstrap_cindex.py`: paired patient-level bootstrap with fixed seed 42. Output
  `results/bootstrap_results.json` reproduces every interval reported in the manuscript.
- `code/compare_to_archive.py`: checks a fresh rerun against `results/` with explicit tolerances.
- `code/repeated_cv.py`: repeated stratified CV on the training cohort only, wider alpha grid, preprocessing
  and gene selection refit inside each fold. Not yet run on the real data.
- `code/record_environment.py`, `requirements.txt`, `.gitignore`, `LICENSE` (MIT, author named),
  `documentation/DATA_ACCESS.md`, `manuscript/SUGGESTED_EDITS.md`.

**Changed**
- `code/frozen_refit.py`: output names and metrics layout now match the archive (`clinical`, `combined`,
  `combined_sensitivity`; `*_risk` columns; flat `model_metrics.json`); logs library versions. Modelling logic
  is unchanged.
- Renamed result files for clarity (contents unchanged): `cv_cv_alpha_summary.csv` ->
  `exploratory_3fold_cv_alpha_summary.csv`, `cv_cv_fold_results.csv` -> `exploratory_3fold_cv_fold_results.csv`,
  `cv_nested_test_metrics.csv` -> `exploratory_3fold_cv_selected_alpha_test_metrics.csv`.
- README and audit rewritten to report the penalty-sensitivity finding and what has and has not been run.

**Removed**
- `code/__pycache__/` (compiled bytecode should not be committed).

**Not changed in the original archive**
- Archived results and figures; manuscript v4 remains archived.

## 2026-11-13 (submission-ready documentation revision)
- Reconciled README, limitations, manuscript source and execution status without changing archived numerical results.
- Added an editable manuscript and PDF labeled as a research draft.
- Corrected stale references to the input hash manifest and license.
- Kept unresolved provenance, tuning and validation limitations explicit.

## 2026-11-13 (pre-publication update)
- Replaced the manuscript with the v4 PDF (`manuscript/TCGA_BRCA_Manuscript_v4.pdf`); removed the superseded
  v3 Word and PDF files so the public package shows one version.
- Added `documentation/EXECUTION_STATUS.md` and `documentation/RAW_GDC_BUILD_STATUS.md` (referenced by the manuscript).
- LICENSE copyright holder set to the manuscript author.
- Rewrote `manuscript/SUGGESTED_EDITS.md` to the edits that remain for v4.
