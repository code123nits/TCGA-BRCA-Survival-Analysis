"""Harrell's concordance index, written to follow scikit-survival's
``concordance_index_censored`` conventions, so results can be recomputed from
saved predictions without installing scikit-survival.

Conventions (as in scikit-survival):
  * A pair (i, j) is comparable when i has an event and either time_j > time_i,
    or time_j == time_i and j is censored. Pairs where both samples have an
    event at the same time are not comparable.
  * Higher risk score should mean shorter survival.
  * Tied risk scores count as 0.5.

If scikit-survival is installed, ``harrell_c`` uses it. Otherwise it uses the
NumPy implementation below. For the archived test predictions, the NumPy
implementation reproduces the C-indices in results/model_metrics.json exactly
(absolute difference 0.0 for all three models).
"""
import numpy as np


def harrell_c_numpy(event, time, risk):
    event = np.asarray(event, dtype=bool)
    time = np.asarray(time, dtype=float)
    risk = np.asarray(risk, dtype=float)
    num = 0.0
    den = 0.0
    for i in np.flatnonzero(event):
        comparable = (time > time[i]) | ((time == time[i]) & ~event)
        comparable[i] = False
        n_comp = comparable.sum()
        if n_comp == 0:
            continue
        r = risk[comparable]
        num += (r < risk[i]).sum() + 0.5 * (r == risk[i]).sum()
        den += n_comp
    if den == 0:
        return float("nan")
    return num / den


def harrell_c(event, time, risk):
    try:
        from sksurv.metrics import concordance_index_censored
        return float(concordance_index_censored(
            np.asarray(event, dtype=bool), np.asarray(time, dtype=float),
            np.asarray(risk, dtype=float))[0])
    except ImportError:
        return float(harrell_c_numpy(event, time, risk))
