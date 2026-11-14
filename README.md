# TCGA-BRCA exploratory survival prediction: does RNA add to age and stage?

**Status:** exploratory analysis of public data. Not peer reviewed, not journal-ready, not externally validated.
The held-out test patients have been examined repeatedly and should not be treated as an untouched test set.

## Question and design
Do 500 variance-selected RNA-expression features improve held-out overall-survival discrimination beyond age and
pathological stage? Elastic-net Cox models (scikit-survival `CoxnetSurvivalAnalysis`, `l1_ratio = 0.5`) were fitted on
857 training patients and scored on 215 held-out patients (30 deaths) from 1,072 matched TCGA-BRCA patients
(149 deaths). Gene selection and all scaling use the training cohort only.

## Main results (archived frozen models, fixed penalties)
| Model | alpha | Test C-index | 95% bootstrap CI |
|---|---|---|---|
| Clinical (age + stage) | 0.001 | 0.6712 | 0.5243 to 0.8040 |
| Clinical + 500 genes | 0.03 | 0.6175 | 0.4848 to 0.7420 |
| Clinical + 500 genes | 0.05 | 0.5854 | 0.4380 to 0.7363 |

Paired differences from the clinical model: -0.0537 (-0.2126 to 0.0972) and -0.0858 (-0.2379 to 0.0481). Both
intervals include zero. Intervals resample fixed predictions and exclude training, gene-selection and penalty
uncertainty. Interpretation: the archived comparisons do not establish an improvement from adding RNA features; the wide intervals and penalty sensitivity prevent a firm conclusion about incremental value.

## Penalty sensitivity (exploratory)
When penalties were chosen by 3-fold cross-validation inside the training cohort instead, the clinical model's test
C-index was 0.5832 (alpha 0.01) and the combined model's was 0.5854 (alpha 0.05), a tie. The ordering in the table above
therefore depends on the penalty chosen and on only 30 test events. The 3-fold run used a coarse grid (best clinical
alpha at the grid edge, frozen alpha 0.001 not included); `code/repeated_cv.py` is provided to redo it properly and
has not yet been run on the real data.

## Repository layout
- `code/frozen_refit.py`: refit the three fixed-penalty models (needs scikit-survival and the input files).
- `code/compare_to_archive.py`: check a rerun against `results/` with explicit tolerances.
- `code/bootstrap_cindex.py`: paired bootstrap from saved predictions (seed 42). Reproduces the manuscript intervals
  exactly; output in `results/bootstrap_results.json`.
- `code/repeated_cv.py`: repeated stratified CV on training data only (not yet run on real data).
- `code/cindex_utils.py`: Harrell's C-index (scikit-survival conventions).
- `code/record_environment.py`: writes exact library versions.
- `results/`: archived metrics, saved test predictions with pseudonymous row labels, selected genes, clinical feature list, exploratory 3-fold CV outputs.
- `figures/`, `manuscript/`: three figures and the current research manuscript draft (PDF, DOCX, Markdown).
- `documentation/`: `DATA_ACCESS.md` (input files and sources) and `AUDIT_AND_LIMITATIONS.md`.

## Before running
The processed inputs are not supplied in this repository. The raw-data construction script and GDC manifest are also absent. A complete rerun therefore requires the original processed inputs, and an end-to-end reconstruction is not possible from this archive alone. Do not download an untrusted pickle: only load a file you generated or trust.

## Reproduce
```bash
pip install -r requirements.txt          # then record exact versions: python code/record_environment.py
# 1. check numbers from saved predictions (no data or scikit-survival needed)
python code/bootstrap_cindex.py --predictions results/test_predictions.csv
# 2. refit the frozen models (needs the input files described in documentation/DATA_ACCESS.md)
python code/frozen_refit.py --data /path/to/data --output rerun_results
python code/compare_to_archive.py --rerun rerun_results --archive results
# 3. cross-validated penalty selection on training data only
python code/repeated_cv.py --data /path/to/data --output results/repeated_cv
```
Do not describe the frozen models as reproduced until step 2 reports `ALL CHECKS PASSED`.

## What has and has not been verified
- Verified: the C-indices in `results/model_metrics.json` and every manuscript interval can be regenerated from
  `results/test_predictions.csv`.
- Not verified: a fresh model fit (never run in this package), the raw-to-matrix build, and the original basis for
  the three penalty values.

## Interpretation and scope
The saved predictions allow an audit of reported discrimination and conditional bootstrap uncertainty, not an independent check of preprocessing, fitting, or test-set independence. The repeated-CV script is a proposed sensitivity analysis, not a completed result. Test-set inspection in earlier work means the label "held-out" describes the original split, not a pristine confirmatory evaluation.

## Not claimed
External validation, calibration, time-dependent performance, clinical utility, treatment effects, a proportional-hazards
check, or a pristine held-out test set. A SCAN-B (GSE96058) feasibility check was done; no external scores were computed.

## Data and license
Inputs are public (NCI GDC TCGA-BRCA; GEO GSE96058) and are not redistributed; see `documentation/DATA_ACCESS.md`.
Code is under the MIT License (copyright holder set in `LICENSE`).

## Processed input integrity
The `documentation/input_sha256_manifest.csv` records SHA-256 hashes of the three supplied processed input files (not distributed in this repository). Check your local copies with `python code/verify_input_integrity.py --data /path/to/inputs`. This is not an end-to-end raw GDC provenance verification.

## Public prediction file
`results/test_predictions.csv` contains pseudonymous sequential sample labels instead of TCGA case identifiers. Survival outcomes and risk predictions are unchanged; this is not a newly fitted model or a fully de-identified clinical dataset. The row-ID mapping is not distributed.
