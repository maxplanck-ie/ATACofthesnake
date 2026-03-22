import numpy as np
from scipy.special import softmax
from sklearn.gaussian_process.kernels import Kernel
from sklearn.gaussian_process.kernels import RBF

class OrdinalKernel(Kernel):
    def __init__(
        self,
        base_kernel,
        unique_times,
        log_spacings=None,
        log_spacings_bounds=(-2.0, 2.0),
    ):
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
            self.log_spacings = orig.copy()
            self.log_spacings[i] += eps
            K_plus = self.base_kernel(
                self.remap_X(X), self.remap_X(Y) if Y is not None else None
            )
            self.log_spacings = orig.copy()
            self.log_spacings[i] -= eps
            K_minus = self.base_kernel(
                self.remap_X(X), self.remap_X(Y) if Y is not None else None
            )
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
                {
                    "base_kernel__" + k: v
                    for k, v in self.base_kernel.get_params(deep=True).items()
                }
            )
        return params


class SlicedRBF(RBF):
    """RBF that only operates on a slice of X."""
    def __init__(self, indices, length_scale=1.0, length_scale_bounds=(1e-3, 1e3)):
        self.indices = indices
        super().__init__(length_scale=length_scale, length_scale_bounds=length_scale_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        Xs = X[:, self.indices]
        Ys = Y[:, self.indices] if Y is not None else None
        return super().__call__(Xs, Ys, eval_gradient=eval_gradient)

    def diag(self, X):
        return super().diag(X[:, self.indices])

    def get_params(self, deep=True):
        params = super().get_params(deep=deep)
        params['indices'] = self.indices
        return params


class SlicedOrdinalKernel(OrdinalKernel):
    """OrdinalKernel that extracts only column 0 (time) from the full X.
    All other columns are discarded"""
    def remap_X(self, X):
        t = self.remap_time(X[:, 0])
        return t[:, None]