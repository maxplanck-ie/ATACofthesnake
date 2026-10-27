import numpy as np
import pandas as pd
import pytest
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C

from aos import gp_fitter
from aos.gp_custom_kernels import OrdinalKernel

# 12 samples: 4 timepoints x 3 replicates.
TIME = np.repeat([0.0, 1.0, 2.0, 3.0], 3)[:, None]
UNIQUE_TIMES = np.array([0.0, 1.0, 2.0, 3.0])
COV = np.tile([0.0, 1.0], 6)[:, None]
NUIS = np.tile([0.0, 0.0, 1.0], 4)[:, None]
GROUP2 = np.tile([0.0, 1.0], 6)[:, None]
GROUP3 = np.tile([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], (4, 1))
Y = np.random.default_rng(0).normal(size=12)
ALPHA = 0.1
PERMS = 50


def sklearn_lml(kernel, X, y):
    gp = GaussianProcessRegressor(kernel=kernel, alpha=ALPHA, optimizer=None)
    return gp.fit(X, y).log_marginal_likelihood()


def old_perm_lr_gp(kernel, X, lml_null):
    """Reference: the original loop, permuting the time column."""
    lr = []
    for p in gp_fitter._permutations(len(X), PERMS):
        X_perm = np.hstack([X[p, :1], X[:, 1:]])
        lr.append(2 * (sklearn_lml(kernel, X_perm, Y) - lml_null))
    return np.array(lr)


def old_perm_lr_interaction(kernel_null, kernel_full, X):
    """Reference: Freedman-Lane loop, permuting the null model's residuals."""
    gp_null = GaussianProcessRegressor(kernel=kernel_null, alpha=ALPHA, optimizer=None)
    y_hat = gp_null.fit(X, Y).predict(X)
    resid = Y - y_hat
    lr = []
    for p in gp_fitter._permutations(len(X), PERMS):
        y_perm = y_hat + resid[p]
        lml_full = sklearn_lml(kernel_full, X, y_perm)
        lml_null = sklearn_lml(kernel_null, X, y_perm)
        lr.append(2 * (lml_full - lml_null))
    return np.array(lr)


def interaction_kernels(fit_type, int_cols, nuis_cols):
    k_t, k_i, k_n = gp_fitter._make_kernel_triple(
        fit_type, [0], int_cols, nuis_cols, UNIQUE_TIMES
    )
    kernel_null = k_t + k_i if k_n is None else k_t + k_i + k_n
    kernel_null = C(1.5) * kernel_null
    k_t2, k_i2, _ = gp_fitter._make_kernel_triple(
        fit_type, [0], int_cols, nuis_cols, UNIQUE_TIMES
    )
    kernel_full = kernel_null + C(0.5) * (k_t2 * k_i2)
    return kernel_null, kernel_full


def test_perm_gp_continuous():
    kernel = C(2.0) * RBF([1.5])
    new = gp_fitter._perm_lr_gp(kernel, TIME, Y, ALPHA, -10.0, PERMS)
    old = old_perm_lr_gp(kernel, TIME, -10.0)
    np.testing.assert_allclose(new, old, rtol=1e-7, atol=1e-8)


def test_perm_gp_continuous_covariates():
    X = np.hstack([TIME, COV])
    kernel = C(2.0) * RBF([1.5, 0.7])
    new = gp_fitter._perm_lr_gp(kernel, X, Y, ALPHA, -10.0, PERMS)
    old = old_perm_lr_gp(kernel, X, -10.0)
    np.testing.assert_allclose(new, old, rtol=1e-7, atol=1e-8)


def test_perm_gp_ordinal():
    kernel = OrdinalKernel(C(2.0) * RBF([1.5]), UNIQUE_TIMES, [0.3, -0.2, 0.1])
    new = gp_fitter._perm_lr_gp(kernel, TIME, Y, ALPHA, -10.0, PERMS)
    old = old_perm_lr_gp(kernel, TIME, -10.0)
    np.testing.assert_allclose(new, old, rtol=1e-7, atol=1e-8)


def test_perm_gp_ordinal_covariates():
    X = np.hstack([TIME, COV])
    kernel = OrdinalKernel(C(2.0) * RBF([1.5, 0.7]), UNIQUE_TIMES, [0.3, -0.2, 0.1])
    new = gp_fitter._perm_lr_gp(kernel, X, Y, ALPHA, -10.0, PERMS)
    old = old_perm_lr_gp(kernel, X, -10.0)
    np.testing.assert_allclose(new, old, rtol=1e-7, atol=1e-8)


def test_perm_interaction_continuous_nuisance():
    X = np.hstack([TIME, GROUP2, NUIS])
    kernel_null, kernel_full = interaction_kernels("continuous", [1], [2])
    new = gp_fitter._perm_lr_gp_interaction(
        kernel_null, kernel_full, X, Y, ALPHA, PERMS
    )
    old = old_perm_lr_interaction(kernel_null, kernel_full, X)
    np.testing.assert_allclose(new, old, rtol=1e-7, atol=1e-8)


def test_perm_interaction_ordinal_three_levels():
    X = np.hstack([TIME, GROUP3])
    kernel_null, kernel_full = interaction_kernels("ordinal", [1, 2], None)
    new = gp_fitter._perm_lr_gp_interaction(
        kernel_null, kernel_full, X, Y, ALPHA, PERMS
    )
    old = old_perm_lr_interaction(kernel_null, kernel_full, X)
    np.testing.assert_allclose(new, old, rtol=1e-7, atol=1e-8)


def test_non_separable_kernel_raises():
    X = np.hstack([TIME, COV])
    kernel = RBF([1.0, 0.5]) + RBF([0.3, 2.0])
    with pytest.raises(ValueError):
        gp_fitter._perm_lr_gp(kernel, X, Y, ALPHA, -10.0, PERMS)


def test_fit_gp_pvalue():
    paramdic = {
        "cov_encoded": pd.DataFrame({"cov": COV.ravel()}),
        "perms": 20,
        "time": TIME,
        "time_grid": UNIQUE_TIMES[:, None],
        "unique_times": UNIQUE_TIMES,
        "fit_type": "ordinal",
        "alpha": ALPHA,
    }
    p_value = gp_fitter.fit_gp(Y, paramdic)[-1]
    assert 0 < p_value <= 1
    paramdic["perms"] = 0
    assert gp_fitter.fit_gp(Y, paramdic)[-1] == 1


def test_fit_gp_interaction_pvalue():
    paramdic = {
        "time": TIME,
        "int_encoded": pd.DataFrame({"group": GROUP2.ravel()}),
        "nuis_encoded": None,
        "perms": 20,
        "alpha": ALPHA,
        "time_grid": UNIQUE_TIMES[:, None],
        "fit_type": "continuous",
        "unique_times": UNIQUE_TIMES,
    }
    p_value = gp_fitter.fit_gp_interaction(Y, paramdic)[-1]
    assert 0 < p_value <= 1
    paramdic["perms"] = 0
    assert gp_fitter.fit_gp_interaction(Y, paramdic)[-1] == 1
