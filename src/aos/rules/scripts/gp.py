import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
import warnings
from sklearn.exceptions import ConvergenceWarning
from sklearn.cluster import KMeans
from aos.helper import get_elbow
import os
# Avoid overthreading.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"


def fit_peak_gp(time, y, time_grid_pred, covariates=None, perms=100):
    '''
    Fit a gaussian process for a single peak.
    Performs a likelihood ratio against a flat kernel to get an idea on the effect size,
    and a permutation test to obtain a p-value (since distribution of LR statistic is ill-defined).

    time - array-like time points
    y - array-like counts.
    time_grid_pred - 
    '''
    time = np.asarray(time)[:, None]

    if covariates is not None:
        X = np.hstack([time, covariates.values])
        X_null = covariates.values
    else:
        X = time
        X_null = np.ones((len(time), 1))

    kernel = C(1.0, (1e-5, 1e5)) * RBF(
        length_scale=np.ones(X.shape[1]),
        length_scale_bounds=(1e-2, 1e2)
    )

    cte_kernel = C(1.0, (1e-3, 1e3))
    gp_cte = GaussianProcessRegressor(
        kernel=cte_kernel,
        alpha=0.1,
        optimizer=None
    )
    gp_cte.fit(X_null, y)
    lml_cte = gp_cte.log_marginal_likelihood()

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.1,
        n_restarts_optimizer=1
    )
    gp.fit(X, y)
    lml_gp = gp.log_marginal_likelihood()
    
    LR_obs = 2 * (lml_gp - lml_cte)

    kernel_fixed = gp.kernel_
    LR_perm = []

    # Smoothed prediction for clustering later
    if covariates is not None:
        # Reference levels for covariates
        cov_ref = covariates.mean(axis=0).values
        X_test = np.hstack([
            time_grid_pred,
            np.tile(cov_ref, (len(time_grid_pred), 1))
        ])
        
        y_pred, y_std = gp.predict(X_test, return_std=True)
    else:
        y_pred, y_std = gp.predict(time_grid_pred, return_std=True)
    
    rng = np.random.default_rng()
    for _ in range(perms):
        perm_idx = rng.permutation(len(time))

        time_perm = time[perm_idx]

        if covariates is not None:
            cov_perm = covariates.values[perm_idx]
            X_perm = np.hstack([time_perm, cov_perm])
            X_null_perm = cov_perm
        else:
            X_perm = time_perm
            X_null_perm = np.ones((len(time_perm), 1))


        gp_perm = GaussianProcessRegressor(
            kernel=kernel_fixed,
            alpha=0.1,
            optimizer=None,
        )

        gp_cte_perm = GaussianProcessRegressor(
            kernel=cte_kernel,
            alpha=0.1,
            optimizer=None,
        )

        gp_perm.fit(X_perm, y)
        gp_cte_perm.fit(X_null_perm, y)

        lml_gp_p = gp_perm.log_marginal_likelihood()
        lml_cte_p = gp_cte_perm.log_marginal_likelihood()

        LR_perm.append(2 * (lml_gp_p - lml_cte_p))

    LR_perm = np.asarray(LR_perm)
    p_value = (np.sum(LR_perm >= LR_obs) + 1) / (len(LR_perm) + 1)
    return LR_obs, p_value, y_pred, y_std




count_matrix = pd.read_csv(snakemake.input.mat, sep='\t')
peaks = count_matrix[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
samplesheet = pd.read_csv(snakemake.params.samplesheet, sep='\t', index_col=0)
count_matrix = count_matrix[list(samplesheet.index)]
# Normalize
count_matrix_norm = np.log1p( count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6 )

time = samplesheet[snakemake.params.comparison['time']].values
time_grid = np.linspace(time.min(), time.max(), 100)[:, None]
covars = samplesheet.drop(columns=[snakemake.params.comparison['time']])
cov_encoded = pd.get_dummies(covars, drop_first=False)

assert len(time) == count_matrix_norm.shape[1]
assert cov_encoded.shape[0] == count_matrix_norm.shape[1]

def gp_row_test(row):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        LR_obs, p_value, y_pred, y_std = fit_peak_gp(
            time,
            row,
            time_grid,
            cov_encoded
        )

    return pd.Series({
        "LR": LR_obs,
        "pvalue": p_value,
        "y_pred": y_pred,
        "y_std": y_std
    })

results = Parallel(n_jobs=snakemake.threads)(
    delayed(gp_row_test)(row)
    for row in count_matrix_norm.values
)
LRs, pvals, ypreds, ystds = zip(*results)
results_df = pd.DataFrame({
    "LR": LRs,
    "pvalue": pvals,
    "y_pred": ypreds,
    "y_std": ystds
}, index=peaks)

reject, results_df['FDR'], _, _ = multipletests(results_df['pvalue'], method="fdr_bh")
patterns = np.vstack(results_df.loc[results_df["FDR"] < 0.05, "y_pred"])
patterns_scaled = patterns - patterns.mean(axis=1, keepdims=True)
patterns_scaled /= patterns.std(axis=1, keepdims=True) + 1e-8

K_range = range(2, 20)
inertias = []
for k in K_range:
    km = KMeans(n_clusters=k, n_init=20, random_state=42)
    km.fit(patterns_scaled)
    inertias.append(km.inertia_)

K_opt = get_elbow(inertias, K_range)
print(f"Selected K = {K_opt}")
K = K_opt
kmeans = KMeans(n_clusters=K, n_init=50, random_state=42)
labels = kmeans.fit_predict(patterns_scaled)
fig, ax = plt.subplots(nrows=K, figsize=(8,12), tight_layout=True)
for k in range(K):
    ax[k].plot(time_grid, patterns_scaled[labels == k].T, alpha=0.3)
    ax[k].plot(time_grid, patterns_scaled[labels == k].mean(axis=0), color='black', linewidth=2)
    ax[k].set_ylabel(f"Cluster {k}")

fig.savefig(snakemake.params.outputfolder + "/clustered_patterns.png", dpi=300)
sig_results = results_df[results_df['FDR'] < 0.05]
sig_results['labels'] = labels
sig_results.to_csv(snakemake.output.sig_clustered, sep='\t')
results_df.to_csv(snakemake.output.table, sep='\t')