import numpy as np
from scipy.special import softmax
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, Kernel
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

class OrdinalKernel(Kernel):

    def __init__(self, base_kernel, unique_times, log_spacings=None, log_spacings_bounds=(-5.0, 5.0)):
        self.base_kernel = base_kernel
        self.unique_times = np.asarray(unique_times, dtype=float)
        self.log_spacings_bounds = log_spacings_bounds
        n = len(self.unique_times)
        if log_spacings is None:
            self.log_spacings = np.zeros(max(0, n - 1), dtype=float)
        else:
            self.log_spacings = np.atleast_1d(np.asarray(log_spacings, dtype=float))

    def remap_time(self, t_raw):
        spacings = np.atleast_1d(softmax(self.log_spacings))
        remapped_time = np.concatenate([[0.0], np.cumsum(spacings)])
        idx = np.searchsorted(self.unique_times, t_raw.ravel())
        idx = np.clip(idx, 0, len(self.unique_times) - 1)
        return remapped_time[idx]

    def remap_X(self, X):
        """
        remap time column, leave other covariates unchanged (if there are any.)
        """
        time_remapped = self.remap_time(X[:, 0])
        if X.shape[1] > 1:
            return np.column_stack([time_remapped, X[:, 1:]])
        return time_remapped[:, None]

    @property
    def theta(self):
        return np.concatenate([self.log_spacings, self.base_kernel.theta])

    @theta.setter
    def theta(self, value):
        n = len(self.log_spacings)
        self.log_spacings = np.atleast_1d(np.asarray(value[:n], dtype=float))
        self.base_kernel.theta = value[n:]

    @property
    def bounds(self):
        n = len(self.log_spacings)
        warp_bounds = np.full((n, 2), self.log_spacings_bounds)
        return np.vstack([warp_bounds, self.base_kernel.bounds])

    def __call__(self, X, Y=None, eval_gradient=False):
        X_w = self.remap_X(X)
        Y_w = self.remap_X(Y) if Y is not None else None

        if not eval_gradient:
            return self.base_kernel(X_w, Y_w)

        K, dK_base = self.base_kernel(X_w, Y_w, eval_gradient=True)

        n_spacings = len(self.log_spacings)
        n_rows = X_w.shape[0]
        n_cols = Y_w.shape[0] if Y_w is not None else n_rows
        dK_warp = np.zeros((n_rows, n_cols, n_spacings))
        eps = 1e-4
        orig = self.log_spacings.copy()
        for i in range(n_spacings):
            self.log_spacings = orig.copy(); self.log_spacings[i] += eps
            K_plus = self.base_kernel(self.remap_X(X), self.remap_X(Y) if Y is not None else None)
            self.log_spacings = orig.copy(); self.log_spacings[i] -= eps
            K_minus = self.base_kernel(self.remap_X(X), self.remap_X(Y) if Y is not None else None)
            dK_warp[:, :, i] = (K_plus - K_minus) / (2 * eps)
        self.log_spacings = orig

        return K, np.concatenate([dK_warp, dK_base], axis=2)

    def diag(self, X):
        return self.base_kernel.diag(self.remap_X(X))

    def is_stationary(self):
        return False

    def get_params(self, deep=True):
        params = dict(
            base_kernel=self.base_kernel,
            unique_times=self.unique_times,
            log_spacings=self.log_spacings,
            log_spacings_bounds=self.log_spacings_bounds,
        )
        if deep:
            params.update(
                {"base_kernel__" + k: v
                 for k, v in self.base_kernel.get_params(deep=True).items()}
            )
        return params



def fit_peak_gp(time, y, time_grid_pred, covariates=None, perms=100):
    '''
    Fit a gaussian process for a single peak.
    Performs a likelihood ratio against a flat kernel to get an idea on the effect size,
    and a permutation test to obtain a p-value (since distribution of LR statistic is ill-defined).

    time - array-like time points
    y - array-like counts.
    time_grid_pred - 
    '''
    time_col = np.asarray(time, dtype=float)[:, None]
    X, X_null = build_design_matrices(time_col, covariates)

    kernel = C(1.0, (1e-5, 1e5)) * RBF(
        length_scale=np.ones(X.shape[1]),
        length_scale_bounds=(1e-2, 1e2)
    )
    return fit_peak_gp_core(
        time_col,
        y,
        np.asarray(time_grid_pred, dtype=float),
        covariates,
        perms,
        kernel,
        X,
        X_null,
    )


def build_design_matrices(time_col, covariates):
    if covariates is not None:
        X = np.hstack([time_col, covariates.values])
        X_null = covariates.values
    else:
        X = time_col
        X_null = np.ones((len(time_col), 1))
    return X, X_null


def predict_on_grid(gp, time_grid_pred, covariates):
    if covariates is not None:
        cov_ref = covariates.mean(axis=0).values
        X_test = np.hstack([
            time_grid_pred,
            np.tile(cov_ref, (len(time_grid_pred), 1))
        ])
    else:
        X_test = time_grid_pred
    return gp.predict(X_test, return_std=True)


def permutation_lr(time_col, y, covariates, perms, kernel_fixed, cte_kernel):
    rng = np.random.default_rng()
    lr_perm = []
    for _ in range(perms):
        perm_idx = rng.permutation(len(time_col))
        time_perm = time_col[perm_idx]

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

        lr_perm.append(
            2 * (
                gp_perm.log_marginal_likelihood()
                - gp_cte_perm.log_marginal_likelihood()
            )
        )

    return np.asarray(lr_perm)


def fit_peak_gp_core(
    time_col,
    y,
    time_grid_pred,
    covariates,
    perms,
    kernel,
    X,
    X_null,
    return_fitted_kernel=False,
):
    cte_kernel = C(1.0, (1e-3, 1e3))
    gp_cte = GaussianProcessRegressor(kernel=cte_kernel, alpha=0.1, optimizer=None)
    gp_cte.fit(X_null, y)
    lml_cte = gp_cte.log_marginal_likelihood()

    gp = GaussianProcessRegressor(kernel=kernel, alpha=0.1, n_restarts_optimizer=1)
    gp.fit(X, y)
    lml_gp = gp.log_marginal_likelihood()

    lr_obs = 2 * (lml_gp - lml_cte)
    y_pred, y_std = predict_on_grid(gp, time_grid_pred, covariates)

    lr_perm = permutation_lr(time_col, y, covariates, perms, gp.kernel_, cte_kernel)
    p_value = (np.sum(lr_perm >= lr_obs) + 1) / (len(lr_perm) + 1)
    if return_fitted_kernel:
        return lr_obs, p_value, y_pred, y_std, gp.kernel_
    return lr_obs, p_value, y_pred, y_std



def fit_peak_gp_ordinal(time, y, time_grid_pred, covariates=None, perms=100):
    """
    Fit a Gaussian process for a single peak with learned time warping.

    The warping learns a monotonic mapping of the observed timepoints so
    that the GP kernel operates in a 'perceptual' time space rather than
    assuming equal spacing between measurements.
    """
    time = np.asarray(time, dtype=float)
    unique_times = np.sort(np.unique(time))
    time_col = time[:, None]
    X, X_null = build_design_matrices(time_col, covariates)

    base_kernel = C(1.0, (1e-5, 1e5)) * RBF(
        length_scale=np.ones(X.shape[1]),
        length_scale_bounds=(1e-2, 1e2)
    )
    kernel = OrdinalKernel(base_kernel, unique_times)

    lr_obs, p_value, y_pred, y_std, fitted_kernel = fit_peak_gp_core(
        time_col,
        y,
        np.asarray(time_grid_pred, dtype=float),
        covariates,
        perms,
        kernel,
        X,
        X_null,
        return_fitted_kernel=True,
    )

    distances = softmax(fitted_kernel.log_spacings)
    return lr_obs, p_value, y_pred, y_std, distances


def run_parallel_gp(
    count_matrix_norm,
    time,
    time_grid,
    cov_encoded,
    fit_func,
    n_jobs,
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
        delayed(gp_row_test)(row)
        for row in count_matrix_norm.values
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
    _, results_df['FDR'], _, _ = multipletests(results_df['pvalue'], method="fdr_bh")
    return results_df

count_matrix = pd.read_csv(snakemake.input.mat, sep='\t')
peaks = count_matrix[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
samplesheet = pd.read_csv(snakemake.params.samplesheet, sep='\t', index_col=0)
    
# Decide the type of analysis
_comparison = snakemake.params.comparison
match _comparison['time_type']:
    case 'continuous': 
        count_matrix = count_matrix[list(samplesheet.index)]
        count_matrix_norm = np.log1p( count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6 )
        time = samplesheet[_comparison['time']].values
        time_grid = np.linspace(time.min(), time.max(), 100)[:, None]
        if len(samplesheet.columns) > 1:
            covars = samplesheet.drop(columns=[_comparison['time']])
            cov_encoded = pd.get_dummies(covars, drop_first=False)
            assert cov_encoded.shape[0] == count_matrix_norm.shape[1]
        else:
            cov_encoded = None
        assert len(time) == count_matrix_norm.shape[1]

        lrs, pvals, ypreds, ystds, _ = run_parallel_gp(
            count_matrix_norm,
            time,
            time_grid,
            cov_encoded,
            fit_peak_gp,
            snakemake.threads,
        )
        results_df = build_results_df(peaks, lrs, pvals, ypreds, ystds)
    case 'ordinal':
        levels = _comparison['order']
        samples_of_interest = [sam for sam in samplesheet.index if samplesheet.loc[sam, _comparison['time']] in levels]
        samplesheet = samplesheet.loc[samples_of_interest]
        count_matrix = count_matrix[ samples_of_interest ]
        count_matrix_norm = np.log1p( count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6 )
        timepoints = levels
        time_str = samplesheet[_comparison['time']].values
        time_map = {t: i for i, t in enumerate(levels)}
        time = np.array([time_map[t] for t in time_str], dtype=float)

        unique_times = np.arange(len(levels), dtype=float)
        time_grid = unique_times[:, None]

        if len(samplesheet.columns) > 1:
            covars = samplesheet.drop(columns=[_comparison['time']])
            cov_encoded = pd.get_dummies(covars, drop_first=False)
            assert cov_encoded.shape[0] == count_matrix_norm.shape[1]
        else:
            cov_encoded = None
        assert len(time) == count_matrix_norm.shape[1]

        lrs, pvals, ypreds, ystds, extras = run_parallel_gp(
            count_matrix_norm,
            time,
            time_grid,
            cov_encoded,
            fit_peak_gp_ordinal,
            snakemake.threads,
            extra_metric_names=["distances"],
        )
        results_df = build_results_df(peaks, lrs, pvals, ypreds, ystds)
        transition_labels = [f"{levels[i]}->{levels[i + 1]}" for i in range(len(levels) - 1)]
        distances_df = pd.DataFrame(
            np.vstack(extras["distances"]),
            index=peaks,
            columns=transition_labels,
        )
        distances_df.to_csv(snakemake.params.outputfolder + '/transition_distances.tsv', sep='\t')




# Cluster significant patterns, plot, and save results
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
    if _comparison['time_type'] == 'ordinal':
        ax[k].set_xticks(unique_times)
        ax[k].set_xticklabels(levels, rotation=90)

fig.savefig(snakemake.params.outputfolder + "/clustered_patterns.png", dpi=300)
sig_results = results_df[results_df['FDR'] < 0.05]
sig_results['labels'] = labels
sig_results.to_csv(snakemake.output.sig_clustered, sep='\t')
results_df.to_csv(snakemake.output.table, sep='\t')