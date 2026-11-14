"""Compare a fresh frozen-model rerun with the archived results.

Usage:
  python code/compare_to_archive.py --rerun rerun_results --archive results

Checks (default tolerances are deliberately small but allow tiny floating-point
differences between library versions):
  * the same 215 test patients, in the same order
  * risk scores for each model: max absolute difference <= --risk-tol
  * C-index for each model: absolute difference <= --cindex-tol
  * number of nonzero coefficients for each model: must be equal
Exit code 0 = all checks pass, 1 = at least one check failed.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODELS = ['clinical', 'combined', 'combined_sensitivity']


def main(rerun, archive, risk_tol, cindex_tol):
    rerun, archive = Path(rerun), Path(archive)
    ok = True

    pr = pd.read_csv(rerun / 'test_predictions.csv')
    pa = pd.read_csv(archive / 'test_predictions.csv')
    if len(pr) != len(pa) or not (pr['patient_id'].astype(str).values == pa['patient_id'].astype(str).values).all():
        print('FAIL  test patients/order differ from archive')
        ok = False
    else:
        print(f'PASS  same {len(pa)} test patients in the same order')
        for m in MODELS:
            col = f'{m}_risk'
            d = float(np.max(np.abs(pr[col].to_numpy() - pa[col].to_numpy())))
            good = d <= risk_tol
            ok &= good
            print(f"{'PASS' if good else 'FAIL'}  {col}: max |difference| = {d:.3e} (tolerance {risk_tol:g})")

    mr = json.loads((rerun / 'model_metrics.json').read_text())
    ma = json.loads((archive / 'model_metrics.json').read_text())
    for m in MODELS:
        d = abs(mr[f'{m}_cindex'] - ma[f'{m}_cindex'])
        good = d <= cindex_tol
        ok &= good
        print(f"{'PASS' if good else 'FAIL'}  {m} C-index: rerun {mr[f'{m}_cindex']:.6f} vs archive "
              f"{ma[f'{m}_cindex']:.6f} (diff {d:.2e})")
        nr, na = mr[f'{m}_nonzero_coefficients'], ma[f'{m}_nonzero_coefficients']
        good = nr == na
        ok &= good
        print(f"{'PASS' if good else 'FAIL'}  {m} nonzero coefficients: rerun {nr} vs archive {na}")

    print('\nALL CHECKS PASSED' if ok else '\nONE OR MORE CHECKS FAILED: do not describe the frozen models as reproduced.')
    return 0 if ok else 1


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--rerun', required=True)
    p.add_argument('--archive', default='results')
    p.add_argument('--risk-tol', type=float, default=1e-6)
    p.add_argument('--cindex-tol', type=float, default=1e-6)
    a = p.parse_args()
    sys.exit(main(a.rerun, a.archive, a.risk_tol, a.cindex_tol))
