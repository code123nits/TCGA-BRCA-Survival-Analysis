"""Paired patient-level bootstrap of held-out Harrell C-index.

Reads saved test-set risk scores (no model refitting) and writes percentile
95% intervals for each model's C-index and for each combined-minus-clinical
difference. The same resampled patients are used for every model in a
replicate. The random seed is fixed so the output is reproducible.

Usage:
  python code/bootstrap_cindex.py \
      --predictions results/test_predictions.csv \
      --output results/bootstrap_results.json \
      --n-boot 2000 --seed 42

Note: this resamples fixed predictions only. It does not capture uncertainty
from training, gene selection or penalty choice.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cindex_utils import harrell_c_numpy

MODELS = {
    "clinical": "clinical_risk",
    "combined_alpha_0.03": "combined_risk",
    "combined_alpha_0.05": "combined_sensitivity_risk",
}


def main(predictions, output, n_boot, seed):
    df = pd.read_csv(predictions)
    event = df["event"].to_numpy().astype(bool)
    time = df["survival_days"].to_numpy().astype(float)
    risks = {name: df[col].to_numpy() for name, col in MODELS.items()}
    n = len(df)

    point = {k: harrell_c_numpy(event, time, v) for k, v in risks.items()}
    rng = np.random.default_rng(seed)
    boots = {k: [] for k in risks}
    diffs = {k: [] for k in risks if k != "clinical"}
    valid = 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if event[idx].sum() < 2:
            continue
        c = {k: harrell_c_numpy(event[idx], time[idx], v[idx]) for k, v in risks.items()}
        if any(np.isnan(x) for x in c.values()):
            continue
        valid += 1
        for k in risks:
            boots[k].append(c[k])
        for k in diffs:
            diffs[k].append(c[k] - c["clinical"])

    def ci(x):
        lo, hi = np.percentile(x, [2.5, 97.5])
        return [float(lo), float(hi)]

    out = {
        "seed": seed,
        "n_boot_requested": n_boot,
        "n_boot_valid": valid,
        "n_test": int(n),
        "n_test_deaths": int(event.sum()),
        "cindex": {k: {"point": float(point[k]), "ci95": ci(boots[k])} for k in risks},
        "difference_vs_clinical": {
            k: {
                "point": float(point[k] - point["clinical"]),
                "ci95": ci(diffs[k]),
                "bootstrap_share_below_zero": float(np.mean(np.array(diffs[k]) < 0)),
            }
            for k in diffs
        },
        "note": "Resamples saved predictions only; excludes training, gene-selection "
                "and penalty-selection uncertainty.",
    }
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--predictions", default="results/test_predictions.csv")
    p.add_argument("--output", default="results/bootstrap_results.json")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=42)
    a = p.parse_args()
    main(a.predictions, a.output, a.n_boot, a.seed)
