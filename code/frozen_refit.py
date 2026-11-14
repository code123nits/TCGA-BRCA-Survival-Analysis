"""Frozen 13-feature TCGA-BRCA exploratory refit (fixed-alpha specifications).

Usage:
  python code/frozen_refit.py --data /path/to/data --output rerun_results

Requires scikit-survival. Run on the same processed expression matrix (see
documentation/DATA_ACCESS.md), NOT on raw GDC RNA files.

Model and output names match the archived files in results/:
  clinical               alpha = 0.001   -> clinical_risk
  combined               alpha = 0.03    -> combined_risk
  combined_sensitivity   alpha = 0.05    -> combined_sensitivity_risk

The modelling logic is unchanged from the original frozen run. Only output
naming, version logging and a metrics file in the archive's layout were added.
The archived survival_3y / survival_5y columns are NOT regenerated here.
After running, check agreement with:
  python code/compare_to_archive.py --rerun rerun_results --archive results
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder, StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cindex_utils import harrell_c  # noqa: E402

STAGES = ['Stage I', 'Stage IA', 'Stage IB', 'Stage II', 'Stage IIA', 'Stage IIB',
          'Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC', 'Stage IV',
          'Stage X', 'Unknown']

# name -> (use RNA features?, alpha)
MODEL_SPECS = {
    'clinical': (False, 0.001),
    'combined': (True, 0.03),
    'combined_sensitivity': (True, 0.05),
}


def main(data, output):
    import sklearn
    import sksurv
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.util import Surv

    data, output = Path(data), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    print(f"python {sys.version.split()[0]} | numpy {np.__version__} | pandas {pd.__version__} "
          f"| scikit-learn {sklearn.__version__} | scikit-survival {sksurv.__version__}")

    df = (pd.read_csv(data / 'tcga_brca_survival_dataset.csv')
          .rename(columns={'Unnamed: 0': 'patient_id'}).set_index('patient_id'))
    exp = pd.read_pickle(data / 'tcga_brca_matched_expression.pkl')
    common = df.index.intersection(exp.index)
    df, exp = df.loc[common], exp.loc[common]
    assert len(df) == 1072 and df.event.sum() == 149, 'Cohort size/events differ from the frozen run.'

    tr, te = train_test_split(df.index.to_numpy(), test_size=.2, stratify=df.event, random_state=42)
    saved = set(pd.read_csv(data / 'test_patient_ids.csv').iloc[:, 0].astype(str))
    assert set(te) == saved, 'Saved test IDs differ; stop.'

    ages = SimpleImputer(strategy='median').fit(df.loc[tr, ['age_at_diagnosis']])
    scaler = StandardScaler().fit(ages.transform(df.loc[tr, ['age_at_diagnosis']]))
    enc = OneHotEncoder(categories=[STAGES], drop='first', handle_unknown='ignore', sparse_output=False)
    enc.fit(df.loc[tr, ['pathologic_stage']].fillna('Unknown'))

    def clinical(ids):
        a = scaler.transform(ages.transform(df.loc[ids, ['age_at_diagnosis']]))
        s = enc.transform(df.loc[ids, ['pathologic_stage']].fillna('Unknown'))
        return np.column_stack([a, s])

    Ctr, Cte = clinical(tr), clinical(te)
    assert Ctr.shape[1] == 13

    # Gene selection and scaling use the training cohort only.
    genes = exp.loc[tr].var().nlargest(500).index
    gs = StandardScaler().fit(exp.loc[tr, genes])
    Xtr = np.column_stack([Ctr, gs.transform(exp.loc[tr, genes])])
    Xte = np.column_stack([Cte, gs.transform(exp.loc[te, genes])])

    ytr = Surv.from_arrays(df.loc[tr, 'event'].astype(bool), df.loc[tr, 'survival_days'].astype(float))
    yte_event = df.loc[te, 'event'].to_numpy().astype(bool)
    yte_time = df.loc[te, 'survival_days'].to_numpy().astype(float)

    clin_names = ['age_at_diagnosis'] + ['pathologic_stage_' + s for s in STAGES[1:]]
    predictions = pd.DataFrame({'patient_id': te, 'event': df.loc[te, 'event'].to_numpy(),
                                'survival_days': df.loc[te, 'survival_days'].to_numpy()})
    flat = {'n_test': int(len(te)), 'n_test_deaths': int(yte_event.sum())}

    for name, (use_rna, alpha) in MODEL_SPECS.items():
        xt, xv = (Xtr, Xte) if use_rna else (Ctr, Cte)
        model = CoxnetSurvivalAnalysis(l1_ratio=.5, alphas=[alpha], max_iter=200000)
        model.fit(xt, ytr)
        risk = model.predict(xv)
        coefs = np.asarray(model.coef_).reshape(-1)
        flat[f'{name}_cindex'] = float(harrell_c(yte_event, yte_time, risk))
        flat[f'{name}_nonzero_coefficients'] = int(np.count_nonzero(coefs))
        predictions[f'{name}_risk'] = risk
        feature_names = clin_names + (list(genes) if use_rna else [])
        pd.DataFrame({'feature': feature_names, 'coefficient': coefs}).to_csv(
            output / f'{name}_coefficients.csv', index=False)

    flat['combined_minus_clinical_cindex'] = flat['combined_cindex'] - flat['clinical_cindex']
    predictions.to_csv(output / 'test_predictions.csv', index=False)
    (output / 'model_metrics.json').write_text(json.dumps(flat, indent=2))
    print(json.dumps(flat, indent=2))
    print('Fresh frozen-model fit. Compare with the archive using code/compare_to_archive.py; '
          'do not assume exact agreement across library versions.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--data', required=True)
    p.add_argument('--output', default='refit_results')
    a = p.parse_args()
    main(a.data, a.output)
