# Execution status (13 Nov 2026)

## Run and verified
- Harrell C-indices recomputed from `results/test_predictions.csv` match `results/model_metrics.json` exactly.
- `code/bootstrap_cindex.py` (seed 42, 2,000 resamples) reproduces every interval in the manuscript;
  output saved as `results/bootstrap_results.json`.
- `code/compare_to_archive.py` passes on an identical copy of the archive and fails on a deliberately altered one.
- All scripts compile.

## Plumbing test only (no real model fitted)
- `code/frozen_refit.py` and `code/repeated_cv.py` were run once on synthetic data with a stand-in for
  scikit-survival, to check file handling and fold logic. This says nothing about reproducing the archived models.

## NOT executed
- A real refit of the three frozen models (scikit-survival was not installed in the environment used) and the
  comparison of that refit with `results/`.
- Repeated 5x5 training-only cross-validation (`code/repeated_cv.py`) on the real data.
- Raw GDC-to-processed-matrix reconstruction (see `RAW_GDC_BUILD_STATUS.md`).

## Note
The manuscript also mentions an input SHA-256 manifest and an input-integrity audit script. They have now been added to this package: `documentation/input_sha256_manifest.csv` and
`code/verify_input_integrity.py`. These hashes apply to the three uploaded processed files,
not the raw GDC source files.

## To complete
1. `pip install -r requirements.txt`, then `python code/record_environment.py`.
2. `python code/frozen_refit.py --data <data_dir> --output rerun_results`
3. `python code/compare_to_archive.py --rerun rerun_results --archive results`
4. Only if step 3 prints ALL CHECKS PASSED, update the manuscript and README to say the frozen models were reproduced.

## Draft revision note
The revised manuscript and documentation were prepared from the archived files. Editing text does not constitute a model refit, a raw-data reconstruction, or new validation.
