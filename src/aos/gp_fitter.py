import sys
import time
import warnings
from datetime import timedelta

import numpy as np
from scipy.special import softmax
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C

from aos.gp_custom_kernels import OrdinalKernel, SlicedOrdinalKernel, SlicedRBF

_PERMUTATION_SEED = 1337
_PERM_CHUNK = 250


def _permutations(n, perms):
    rng = np.random.default_rng(_PERMUTATION_SEED)
    return np.array(
        [rng.permutation(n) for _ in range(perms)], dtype=np.intp
    ).reshape(perms, n)


def _lml_stack(K, y, alpha):
    """sklearn's GPR log marginal likelihood for a stack of kernel matrices."""
    n = K.shape[-1]
    K[..., np.arange(n), np.arange(n)] += alpha
    L = np.linalg.cholesky(K)
    z = np.linalg.solve(L, y[:, None])[..., 0]
    return (
        -0.5 * np.sum(z**2, axis=-1)
        - np.sum(np.log(np.diagonal(L, axis1=-2, axis2=-1)), axis=-1)
        - 0.5 * n * np.log(2 * np.pi)
    )


def _perm_lml(build_K, y, alpha, perm_idx):
    lml = np.empty(len(perm_idx))
    for s in range(0, len(perm_idx), _PERM_CHUNK):
        lml[s : s + _PERM_CHUNK] = _lml_stack(
            build_K(perm_idx[s : s + _PERM_CHUNK]), y, alpha
        )
    return lml


def _split_time_kernel(kernel, X):
    """Split a fitted kernel into K(X) = K_time * K_cov. Time is column 0."""
    X = np.asarray(X, dtype=float)
    X_time = X.copy()
    X_time[:, 1:] = 0.0
    X_cov = X.copy()
    X_cov[:, 0] = X[0, 0]
    K_time = kernel(X_time)
    K_cov = kernel(X_cov) / kernel.diag(X[:1])[0]
    if not np.allclose(K_time * K_cov, kernel(X)):
        raise ValueError("Fitted GP kernel is not separable in time and covariates.")
    return K_time, K_cov


def _gather_levels(K_lv, codes, p):
    c = codes[p]
    idx = np.arange(len(codes))
    return K_lv[c[:, :, None], c[:, None, :], idx[:, None], idx[None, :]]


def _level_pair_kernels(kernel, X, cols, levels, codes):
    """K_lv[g, h] = kernel(X with cols set to level g, X with cols set to level h)."""
    X_lv = []
    for level in levels:
        X_l = X.copy()
        X_l[:, cols] = level
        X_lv.append(X_l)
    K_lv = np.array([[kernel(X_g, X_h) for X_h in X_lv] for X_g in X_lv])
    identity = np.arange(X.shape[0])[None, :]
    if not np.allclose(_gather_levels(K_lv, codes, identity)[0], kernel(X)):
        raise ValueError("Level-pair kernels do not match the fitted GP kernel.")
    return K_lv


def _perm_lr_gp(kernel, X, y, alpha, lml_null, perms):
    K_time, K_cov = _split_time_kernel(kernel, X)
    lml = _perm_lml(
        lambda p: K_time[p[:, :, None], p[:, None, :]] * K_cov,
        y,
        alpha,
        _permutations(len(X), perms),
    )
    return 2 * (lml - lml_null)


def _perm_lr_gp_interaction(kernel_null, kernel_full, X, y, alpha, int_cols, perms):
    X = np.asarray(X, dtype=float)
    levels, codes = np.unique(X[:, int_cols], axis=0, return_inverse=True)
    codes = codes.reshape(-1)
    K_null = _level_pair_kernels(kernel_null, X, int_cols, levels, codes)
    K_full = _level_pair_kernels(kernel_full, X, int_cols, levels, codes)
    perm_idx = _permutations(len(X), perms)
    lml_null = _perm_lml(lambda p: _gather_levels(K_null, codes, p), y, alpha, perm_idx)
    lml_full = _perm_lml(lambda p: _gather_levels(K_full, codes, p), y, alpha, perm_idx)
    return np.maximum(0.0, 2 * (lml_full - lml_null))


def fit_gp(y, paramdic):
    cov_encoded = paramdic["cov_encoded"]
    perms = paramdic["perms"]
    time = paramdic["time"]
    time_grid = paramdic["time_grid"]
    unique_times = paramdic["unique_times"]
    fit_type = paramdic["fit_type"]
    alpha = paramdic["alpha"]

    if cov_encoded is not None:
        X = np.hstack([time, cov_encoded.values])
        X_null = cov_encoded.values
    else:
        X = time
        X_null = np.ones((len(time), 1))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)

        # Set kernels
        if cov_encoded is not None:
            cte_kernel = C(1.0, (1e-5, 1e5)) * RBF(
                length_scale=np.ones(cov_encoded.shape[1]),
                length_scale_bounds=(1e-3, 1e3),
            )
        else:
            cte_kernel = C(1.0, (1e-3, 1e3))
        if fit_type == "ordinal":
            base_kernel = C(1.0, (1e-5, 1e5)) * RBF(
                length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-3, 1e3)
            )
            kernel = OrdinalKernel(base_kernel, unique_times)
        else:
            assert fit_type == "continuous"
            kernel = C(1.0, (1e-5, 1e5)) * RBF(
                length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-3, 1e3)
            )
        gp_cte = GaussianProcessRegressor(
            kernel=cte_kernel, alpha=alpha, optimizer=None
        )
        gp_cte.fit(X_null, y)
        lml_cte = gp_cte.log_marginal_likelihood()

        gp = GaussianProcessRegressor(
            kernel=kernel, alpha=alpha, n_restarts_optimizer=5
        )
        gp.fit(X, y)
        lml_gp = gp.log_marginal_likelihood()

        lr_obs = 2 * (lml_gp - lml_cte)

        # PREDICTION
        if cov_encoded is not None:
            cov_ref = cov_encoded.mean(axis=0).values
            X_test = np.hstack([time_grid, np.tile(cov_ref, (len(time_grid), 1))])
        else:
            X_test = time_grid
        y_pred, y_std = gp.predict(X_test, return_std=True)

        # DISTANCES
        if fit_type == "ordinal":
            distances = softmax(gp.kernel_.log_spacings)
        else:
            distances = np.array([])

        # PERM - null model does not see time, so lml_cte holds for every permutation.
        lr_perm = _perm_lr_gp(gp.kernel_, X, y, alpha, lml_cte, perms)
        p_value = (np.sum(lr_perm >= lr_obs) + 1) / (len(lr_perm) + 1)
        return y_pred, y_std, distances, lr_obs, p_value


def _make_kernel_triple(fit_type, t_idx, i_idx, n_idx, unique_times):
    if fit_type == "ordinal":
        base = C(1.0, (1e-5, 1e5)) * RBF(
            length_scale=1.0, length_scale_bounds=(1e-3, 1e3)
        )
        k_t = SlicedOrdinalKernel(base, unique_times)
    else:
        k_t = SlicedRBF(
            t_idx, length_scale=np.ones(len(t_idx)), length_scale_bounds=(1e-3, 1e3)
        )
    k_i = SlicedRBF(
        i_idx, length_scale=np.full(len(i_idx), 0.5), length_scale_bounds="fixed"
    )
    k_n = (
        SlicedRBF(
            n_idx, length_scale=np.ones(len(n_idx)), length_scale_bounds=(1e-3, 1e3)
        )
        if n_idx is not None
        else None
    )
    return k_t, k_i, k_n


def fit_gp_interaction(y, paramdic):
    time = paramdic["time"]
    int_encoded = paramdic["int_encoded"]
    nuis_encoded = paramdic["nuis_encoded"]
    perms = paramdic["perms"]
    alpha = paramdic["alpha"]
    time_grid = paramdic["time_grid"]
    fit_type = paramdic["fit_type"]
    unique_times = paramdic["unique_times"]

    if nuis_encoded is not None:
        X = np.hstack([time, int_encoded, nuis_encoded])
    else:
        X = np.hstack([time, int_encoded])

    t_idx = list(range(time.shape[1]))
    i_start = time.shape[1]
    i_idx = list(range(i_start, i_start + int_encoded.shape[1]))
    n_idx = (
        list(
            range(
                i_start + int_encoded.shape[1],
                i_start + int_encoded.shape[1] + nuis_encoded.shape[1],
            )
        )
        if nuis_encoded is not None
        else None
    )

    k_t_null, k_i_null, k_n_null = _make_kernel_triple(
        fit_type, t_idx, i_idx, n_idx, unique_times
    )
    kernel_null = k_t_null + k_i_null
    if k_n_null is not None:
        kernel_null += k_n_null
    kernel_null = C(1.0, (1e-5, 1e5)) * kernel_null

    k_t_int, k_i_int, _ = _make_kernel_triple(
        fit_type, t_idx, i_idx, n_idx, unique_times
    )
    kernel_full = kernel_null + C(1.0, (1e-2, 1e2)) * (k_t_int * k_i_int)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        gp_null = GaussianProcessRegressor(
            kernel=kernel_null, alpha=alpha, n_restarts_optimizer=5
        )
        gp_null.fit(X, y)
        lml_null = gp_null.log_marginal_likelihood()

        gp_full = GaussianProcessRegressor(
            kernel=kernel_full, alpha=alpha, n_restarts_optimizer=5
        )
        gp_full.fit(X, y)
        lml_full = gp_full.log_marginal_likelihood()

    lr_obs = max(0.0, 2 * (lml_full - lml_null))

    # PREDICTION

    ref_level = np.zeros((1, int_encoded.shape[1]))
    int_levels = np.vstack([ref_level, np.eye(int_encoded.shape[1])])

    if nuis_encoded is not None:
        nuis_ref = np.mean(nuis_encoded, axis=0, keepdims=True)

    y_pred_list, y_std_list = [], []
    for int_ref in int_levels:
        int_ref = int_ref[np.newaxis, :]
        if nuis_encoded is not None:
            X_test = np.hstack(
                [
                    time_grid,
                    np.repeat(int_ref, len(time_grid), axis=0),
                    np.repeat(nuis_ref, len(time_grid), axis=0),
                ]
            )
        else:
            X_test = np.hstack(
                [
                    time_grid,
                    np.repeat(int_ref, len(time_grid), axis=0),
                ]
            )
        yp, ys = gp_full.predict(X_test, return_std=True)
        y_pred_list.append(yp)
        y_std_list.append(ys)

    y_pred_arr = np.vstack(y_pred_list)
    y_std_arr = np.vstack(y_std_list)

    # ORDINAL DISTANCES
    if fit_type == "ordinal":

        def _find_ordinal(k):
            if isinstance(k, SlicedOrdinalKernel):
                return k
            for attr in ("k1", "k2"):
                if hasattr(k, attr):
                    found = _find_ordinal(getattr(k, attr))
                    if found is not None:
                        return found
            return None

        ok = _find_ordinal(gp_full.kernel_)
        distances = softmax(ok.log_spacings) if ok is not None else np.array([])
    else:
        distances = np.array([])

    # PERMUTATION TEST
    lr_perm = _perm_lr_gp_interaction(
        gp_null.kernel_, gp_full.kernel_, X, y, alpha, i_idx, perms
    )
    p_value = (np.sum(lr_perm >= lr_obs) + 1) / (len(lr_perm) + 1)
    return y_pred_arr, y_std_arr, distances, lr_obs, p_value


def log_progress(results, total, label, logfile, every=60):
    """Yield results, writing progress and ETA to stderr and logfile."""
    start = last = time.monotonic()
    with open(logfile, "w") as log:
        _line = f"{label}: fitting {total} peaks"
        print(_line, file=sys.stderr, flush=True)
        print(_line, file=log, flush=True)
        for done, result in enumerate(results, 1):
            yield result
            now = time.monotonic()
            if now - last >= every or done == total:
                last = now
                elapsed = now - start
                eta = elapsed / done * (total - done)
                _line = (
                    f"{label}: {done}/{total} peaks ({done / total:.0%}) | "
                    f"{done / elapsed:.1f} peaks/s | "
                    f"elapsed {timedelta(seconds=round(elapsed))} | "
                    f"ETA {timedelta(seconds=round(eta))}"
                )
                print(_line, file=sys.stderr, flush=True)
                print(_line, file=log, flush=True)
