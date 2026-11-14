"""Repeated stratified cross-validation on the TRAINING cohort only.

Purpose: choose the elastic-net penalty (alpha) for each model type without
touching the held-out test patients, and quantify how much CV C-index varies
across folds and repeats.

Usage:
  python code/repeated_cv.py --data /path/to/data --output results/repeated_cv \
      --n-splits 5 --n-repeats 5 --seed 42

What it does
  * Re-creates the original 80:20 split (random_state=42) and asserts that the
    test IDs match the saved ones. Test patients are never used here.
  * Inside EVERY fold it re-fits the age imputer/scaler, the stage encoder, the
    500-gene variance selection and the gene scaler on that fold's training
    part only, so nothing leaks from the validation part.
  * Evaluates a wider alpha grid than the earlier 3-fold exploratory run,
    including the frozen values (0.001, 0.03, 0.05) and smaller values.
  * Reports mean, SD and number of valid folds per model and alpha, plus the
    alpha with the best mean and the largest alpha within one standard error
    of the best (a more conservative choice).

Status: written to run with scikit-survival, but NOT yet executed on the real
TCGA-BRCA data in this package. Run it once and inspect the output before
using any number in a manuscript.
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.model_selection import RepeatedStratifiedKFold, train_test_split
from sklearn.preprocessing import OneHotEncoder, StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cindex_utils import harrell_c  # noqa: E402

STAGES = ['Stage I', 'Stage IA', 'Stage IB', 'Stage II', 'Stage IIA', 'Stage IIB',
          'Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC', 'Stage IV',
          'Stage X', 'Unknown']
ALPHAS = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.05, 0.1]


def build_features(df, exp, fit_ids, apply_ids, use_rna):
    """Fit all transforms on fit_ids only; return design matrices for both sets."""
    imp = SimpleImputer(strategy='median').fit(df.loc[fit_ids, ['age_at_diagnosis']])
    sc = StandardScaler().fit(imp.transform(df.loc[fit_ids, ['age_at_diagnosis']]))
    enc = OneHotEncoder(categories=[STAGES], drop='first', handle_unknown='ignore', sparse_output=False)
    enc.fit(df.loc[fit_ids, ['pathologic_stage']].fillna('Unknown'))

    def clin(ids):
        a = sc.transform(imp.transform(df.loc[ids, ['age_at_diagnosis']]))
        s = enc.transform(df.loc[ids, ['pathologic_stage']].fillna('Unknown'))
        return np.column_stack([a, s])

    Xa, Xb = clin(fit_ids), clin(apply_ids)
    if use_rna:
        genes = exp.loc[fit_ids].var().nlargest(500).index
        gs = StandardScaler().fit(exp.loc[fit_ids, genes])
        Xa = np.column_stack([Xa, gs.transform(exp.loc[fit_ids, genes])])
        Xb = np.column_stack([Xb, gs.transform(exp.loc[apply_ids, genes])])
    return Xa, Xb


def main(data, output, n_splits, n_repeats, seed):
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.util import Surv

    data, output = Path(data), Path(output)
    output.mkdir(parents=True, exist_ok=True)

    df = (pd.read_csv(data / 'tcga_brca_survival_dataset.csv')
          .rename(columns={'Unnamed: 0': 'patient_id'}).set_index('patient_id'))
    exp = pd.read_pickle(data / 'tcga_brca_matched_expression.pkl')
    common = df.index.intersection(exp.index)
    df, exp = df.loc[common], exp.loc[common]

    tr, te = train_test_split(df.index.to_numpy(), test_size=.2, stratify=df.event, random_state=42)
    saved = set(pd.read_csv(data / 'test_patient_ids.csv').iloc[:, 0].astype(str))
    assert set(te) == saved, 'Saved test IDs differ; stop.'

    tr_df = df.loc[tr]
    event = tr_df['event'].to_numpy().astype(bool)
    time = tr_df['survival_days'].to_numpy().astype(float)
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)

    rows = []
    for fold_id, (a, b) in enumerate(rskf.split(tr, event)):
        repeat, fold = divmod(fold_id, n_splits)
        fit_ids, val_ids = tr[a], tr[b]
        y_fit = Surv.from_arrays(df.loc[fit_ids, 'event'].astype(bool),
                                 df.loc[fit_ids, 'survival_days'].astype(float))
        ev_val = df.loc[val_ids, 'event'].to_numpy().astype(bool)
        t_val = df.loc[val_ids, 'survival_days'].to_numpy().astype(float)
        for model_name, use_rna in [('clinical', False), ('combined', True)]:
            Xa, Xb = build_features(df, exp, fit_ids, val_ids, use_rna)
            for alpha in ALPHAS:
                status, c, nz = 'ok', np.nan, np.nan
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    try:
                        m = CoxnetSurvivalAnalysis(l1_ratio=.5, alphas=[alpha], max_iter=200000)
                        m.fit(Xa, y_fit)
                        coefs = np.asarray(m.coef_).reshape(-1)
                        nz = int(np.count_nonzero(coefs))
                        if nz == 0:
                            status = 'all_coefficients_zero'   # C-index would be uninformative
                        else:
                            c = harrell_c(ev_val, t_val, m.predict(Xb))
                    except Exception as e:  # keep going; record the failure
                        status = f'error: {type(e).__name__}'
                rows.append(dict(model=model_name, alpha=alpha, repeat=repeat, fold=fold,
                                 cindex=c, nonzero=nz, status=status))
        print(f'finished repeat {repeat + 1}/{n_repeats}, fold {fold + 1}/{n_splits}', flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(output / 'repeated_cv_fold_results.csv', index=False)

    ok = res[res.status == 'ok']
    summ = (ok.groupby(['model', 'alpha'])['cindex']
            .agg(mean_cindex='mean', sd_cindex='std', valid_folds='count').reset_index())
    summ['se_cindex'] = summ['sd_cindex'] / np.sqrt(summ['valid_folds'])
    summ.to_csv(output / 'repeated_cv_alpha_summary.csv', index=False)

    picks = []
    for model_name, g in summ.groupby('model'):
        best = g.loc[g['mean_cindex'].idxmax()]
        within = g[g['mean_cindex'] >= best['mean_cindex'] - best['se_cindex']]
        one_se = within.loc[within['alpha'].idxmax()]
        picks.append(dict(model=model_name, best_alpha=best['alpha'], best_mean_cindex=best['mean_cindex'],
                          one_se_alpha=one_se['alpha'], one_se_mean_cindex=one_se['mean_cindex'],
                          best_alpha_is_grid_edge=bool(best['alpha'] in (min(ALPHAS), max(ALPHAS)))))
    pd.DataFrame(picks).to_csv(output / 'repeated_cv_selected_alpha.csv', index=False)
    print(summ.to_string(index=False))
    print(pd.DataFrame(picks).to_string(index=False))
    print('If best_alpha_is_grid_edge is True, extend the alpha grid before trusting the choice.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--data', required=True)
    p.add_argument('--output', default='results/repeated_cv')
    p.add_argument('--n-splits', type=int, default=5)
    p.add_argument('--n-repeats', type=int, default=5)
    p.add_argument('--seed', type=int, default=42)
    a = p.parse_args()
    main(a.data, a.output, a.n_splits, a.n_repeats, a.seed)
