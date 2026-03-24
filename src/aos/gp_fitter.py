import numpy as np
import warnings
from sklearn.exceptions import ConvergenceWarning

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C
from scipy.special import softmax

from aos.gp_custom_kernels import OrdinalKernel, SlicedRBF, SlicedOrdinalKernel

def fit_gp(y, paramdic):
    cov_encoded = paramdic['cov_encoded']
    perms = paramdic['perms']
    time = paramdic['time']
    time_grid = paramdic['time_grid']
    unique_times = paramdic['unique_times']
    fit_type = paramdic['fit_type']
    
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
                length_scale_bounds=(1e-3, 1e3)
            )
        else:
            cte_kernel = C(1.0, (1e-3, 1e3))
        if fit_type == 'ordinal':
            base_kernel = C(1.0, (1e-5, 1e5)) * RBF(
                    length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-3, 1e3)
                )
            kernel = OrdinalKernel(base_kernel, unique_times)
        else:
            assert fit_type == 'continuous'
            kernel = C(1.0, (1e-5, 1e5)) * RBF(
                length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-3, 1e3)
            )
        gp_cte = GaussianProcessRegressor(kernel=cte_kernel, alpha=0.1, optimizer=None)
        gp_cte.fit(X_null, y)
        lml_cte = gp_cte.log_marginal_likelihood()
        
        gp = GaussianProcessRegressor(kernel=kernel, alpha=0.1, n_restarts_optimizer=5)
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
        if fit_type == 'ordinal':
            distances = softmax(gp.kernel_.log_spacings)
        else:
            distances = np.array([])
        
        # PERM
        rng = np.random.default_rng()
        lr_perm = []
    
        for _ in range(perms):
            perm_idx = rng.permutation(len(time))
            time_perm = time[perm_idx]
        
            if cov_encoded is not None:
                X_perm = np.hstack([time_perm, cov_encoded.values])
                X_null_perm = cov_encoded.values
            else:
                X_perm = time_perm
                X_null_perm = np.ones((len(time_perm), 1))
        
            gp_perm = GaussianProcessRegressor(
                kernel=gp.kernel_,
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
        lr_perm = np.asarray(lr_perm)
        p_value = (np.sum(lr_perm >= lr_obs) + 1) / (len(lr_perm) + 1)
        return y_pred, y_std, distances, lr_obs, p_value

def _make_kernel_triple(fit_type, t_idx, i_idx, n_idx, unique_times):
    if fit_type == 'ordinal':
        base = C(1.0, (1e-5, 1e5)) * RBF(
            length_scale=1.0, length_scale_bounds=(1e-3, 1e3)
        )
        k_t = SlicedOrdinalKernel(base, unique_times)
    else:
        k_t = SlicedRBF(t_idx, length_scale=np.ones(len(t_idx)),
                        length_scale_bounds=(1e-3, 1e3))
    k_i = SlicedRBF(i_idx,
                length_scale=np.full(len(i_idx), 0.5),
                length_scale_bounds="fixed")
    k_n = (SlicedRBF(n_idx, length_scale=np.ones(len(n_idx)),
                     length_scale_bounds=(1e-3, 1e3))
           if n_idx is not None else None)
    return k_t, k_i, k_n

def fit_gp_interaction(y, paramdic):
    time        = paramdic["time"]
    int_encoded = paramdic["int_encoded"]
    nuis_encoded = paramdic["nuis_encoded"]
    perms       = paramdic["perms"]
    time_grid   = paramdic["time_grid"]
    fit_type    = paramdic["fit_type"]
    unique_times = paramdic["unique_times"]

    if nuis_encoded is not None:
        X = np.hstack([time, int_encoded, nuis_encoded])
    else:
        X = np.hstack([time, int_encoded])

    t_idx  = list(range(time.shape[1]))
    i_start = time.shape[1]
    i_idx  = list(range(i_start, i_start + int_encoded.shape[1]))
    n_idx  = (list(range(i_start + int_encoded.shape[1],
                         i_start + int_encoded.shape[1] + nuis_encoded.shape[1]))
              if nuis_encoded is not None else None)

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
            kernel=kernel_null, alpha=0.5, n_restarts_optimizer=3
        )
        gp_null.fit(X, y)
        lml_null = gp_null.log_marginal_likelihood()

        gp_full = GaussianProcessRegressor(
            kernel=kernel_full, alpha=0.5, n_restarts_optimizer=3
        )
        gp_full.fit(X, y)
        lml_full = gp_full.log_marginal_likelihood()

    lr_obs = max(0.0, 2 * (lml_full - lml_null))

    # PREDICTION

    ref_level  = np.zeros((1, int_encoded.shape[1]))
    int_levels = np.vstack([ref_level, np.eye(int_encoded.shape[1])])

    if nuis_encoded is not None:
        nuis_ref = np.mean(nuis_encoded, axis=0, keepdims=True)

    y_pred_list, y_std_list = [], []
    for int_ref in int_levels:
        int_ref = int_ref[np.newaxis, :]
        if nuis_encoded is not None:
            X_test = np.hstack([
                time_grid,
                np.repeat(int_ref,  len(time_grid), axis=0),
                np.repeat(nuis_ref, len(time_grid), axis=0),
            ])
        else:
            X_test = np.hstack([
                time_grid,
                np.repeat(int_ref, len(time_grid), axis=0),
            ])
        yp, ys = gp_full.predict(X_test, return_std=True)
        y_pred_list.append(yp)
        y_std_list.append(ys)

    y_pred_arr = np.vstack(y_pred_list)
    y_std_arr  = np.vstack(y_std_list)

    # ORDINAL DISTANCES
    if fit_type == 'ordinal':
        def _find_ordinal(k):
            if isinstance(k, SlicedOrdinalKernel):
                return k
            for attr in ('k1', 'k2'):
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
    rng = np.random.default_rng()
    lr_perm = []
    for _ in range(perms):
        perm_idx = rng.permutation(len(time))
        int_perm = int_encoded.iloc[perm_idx]
        if nuis_encoded is not None:
            X_perm = np.hstack([time, int_perm, nuis_encoded])
        else:
            X_perm = np.hstack([time, int_perm])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            gp_null_perm = GaussianProcessRegressor(
                kernel=gp_null.kernel_, alpha=0.5, optimizer=None
            )
            gp_full_perm = GaussianProcessRegressor(
                kernel=gp_full.kernel_, alpha=0.5, optimizer=None
            )
            gp_null_perm.fit(X_perm, y)
            gp_full_perm.fit(X_perm, y)
            lr_perm.append(max(0.0, 2 * (
                gp_full_perm.log_marginal_likelihood()
                - gp_null_perm.log_marginal_likelihood()
            )))
    lr_perm = np.asarray(lr_perm)
    p_value = (np.sum(lr_perm >= lr_obs) + 1) / (len(lr_perm) + 1)
    return y_pred_arr, y_std_arr, distances, lr_obs, p_value