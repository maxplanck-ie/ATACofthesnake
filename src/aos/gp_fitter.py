import numpy as np
import warnings
from sklearn.exceptions import ConvergenceWarning

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C
from scipy.special import softmax

from aos.gp_ordinal import OrdinalKernel

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
        
        gp = GaussianProcessRegressor(kernel=kernel, alpha=0.1, n_restarts_optimizer=1)
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
                cov_perm = cov_encoded.values[perm_idx]
                X_perm = np.hstack([time_perm, cov_perm])
                X_null_perm = cov_perm
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