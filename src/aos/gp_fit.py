import numpy as np
from scipy.special import softmax
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C

from aos.gp_ordinal import OrdinalKernel


def build_design_matrices(time_col, covariates):
    if covariates is not None:
        X = np.hstack([time_col, covariates.values])
        X_null = covariates.values
    else:
        X = time_col
        X_null = np.ones((len(time_col), 1))
    return X, X_null


def predict_acc(gp, time_grid_pred, covariates):
    if covariates is not None:
        cov_ref = covariates.mean(axis=0).values
        X_test = np.hstack([time_grid_pred, np.tile(cov_ref, (len(time_grid_pred), 1))])
    else:
        X_test = time_grid_pred
    return gp.predict(X_test, return_std=True)


def fit_peak_gp(time, y, time_grid_pred, covariates=None, perms=1000):
    """
    Fit a gaussian process for a single peak.
    Calculates likelihood ratio against a flat kernel to get an idea on the effect size,
    and a permutation test to obtain a p-value (since distribution of LR statistic is ill-defined).

    time - array-like time points
    y - array-like counts.
    time_grid_pred - array-like time points for prediction
    """
    time_col = np.asarray(time, dtype=float)[:, None]
    X, X_null = build_design_matrices(time_col, covariates)

    kernel = C(1.0, (1e-5, 1e5)) * RBF(
        length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-2, 1e2)
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
    y_pred, y_std = predict_acc(gp, time_grid_pred, covariates)

    lr_perm = permutation_lr(time_col, y, covariates, perms, gp.kernel_, cte_kernel)
    p_value = (np.sum(lr_perm >= lr_obs) + 1) / (len(lr_perm) + 1)
    if return_fitted_kernel:
        return lr_obs, p_value, y_pred, y_std, gp.kernel_
    return lr_obs, p_value, y_pred, y_std


def fit_peak_gp_ordinal(time, y, time_grid_pred, covariates=None, perms=1000):
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
        length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-2, 1e2)
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
            2
            * (
                gp_perm.log_marginal_likelihood()
                - gp_cte_perm.log_marginal_likelihood()
            )
        )

    return np.asarray(lr_perm)
