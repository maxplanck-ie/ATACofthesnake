import warnings

import pandas as pd
from joblib import Parallel, delayed
from sklearn.exceptions import ConvergenceWarning
from statsmodels.stats.multitest import multipletests


def run_parallel_gp(
    count_matrix_norm,
    time,
    time_grid,
    cov_encoded,
    fit_func,
    n_jobs,
    permutations=1000,
    extra_metric_names=None,
):
    def gp_row_test(row):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            fit_out = fit_func(
                time,
                row,
                time_grid,
                cov_encoded,
            )

        return fit_out

    results = Parallel(n_jobs=n_jobs)(
        delayed(gp_row_test)(row) for row in count_matrix_norm.values
    )
    unpacked = list(zip(*results))
    lrs, pvals, ypreds, ystds = unpacked[:4]

    extras = {}
    if extra_metric_names:
        for name, values in zip(extra_metric_names, unpacked[4:]):
            extras[name] = values

    return lrs, pvals, ypreds, ystds, extras


def build_results_df(peaks, lrs, pvals, ypreds, ystds):
    results_df = pd.DataFrame(
        {
            "LR": lrs,
            "pvalue": pvals,
            "y_pred": ypreds,
            "y_std": ystds,
        },
        index=peaks,
    )
    _, results_df["FDR"], _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
    return results_df
