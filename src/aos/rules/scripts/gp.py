import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import os
from aos.gp_fitter import fit_gp
from statsmodels.stats.multitest import multipletests
from rich.progress import track

# Avoid overthreading.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# Snakemake variables
_comparison = snakemake.params.comparison
mat = snakemake.input.mat
ss = snakemake.params.ss
THREADS = snakemake.threads
PERMS = snakemake.params.permutations
GP_TIMESTEPS = snakemake.params.gp_timesteps
table_output = snakemake.output.table
norm_matrix_output = snakemake.output.norm_matrix

if _comparison['time_type'] == 'ordinal':
    # Set up variables
    count_matrix = pd.read_csv(mat, sep="\t")
    peaks = count_matrix[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
    samplesheet = pd.read_csv(ss, sep="\t", index_col=0)
    levels = _comparison["order"]

    # In case of ordinal from cont. data, match case in levels for sake of subsampling
    _dtype = samplesheet[_comparison['time']].dtype
    _levels_case = [_dtype.type(i) for i in _comparison['order']]
    
    samples_of_interest = [
        sam
        for sam in samplesheet.index
        if samplesheet.loc[sam, _comparison["time"]] in _levels_case
    ]
    samplesheet = samplesheet.loc[samples_of_interest]
    count_matrix = count_matrix[samples_of_interest]

    _time_str = samplesheet[_comparison["time"]].values.astype(str)
    _time_map = {t: i for i, t in enumerate(levels)}
    _time = np.array([_time_map[t] for t in _time_str], dtype=float)
    time = _time[:, None]
    unique_times = np.arange(len(levels), dtype=float)
    time_grid = unique_times[:, None]

else:
    assert _comparison['time_type'] == 'continuous'
    # Set up variables
    count_matrix = pd.read_csv(mat, sep="\t")
    peaks = count_matrix[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
    samplesheet = pd.read_csv(ss, sep="\t", index_col=0)

    count_matrix = count_matrix[list(samplesheet.index)]
    
    time = samplesheet[_comparison["time"]].values[:, None]
    time_grid = np.linspace(time.min(), time.max(), GP_TIMESTEPS)[:, None]
    # Just needed in ordinal case, not here.
    unique_times = None

# Common prep - norm counts
count_matrix_norm = np.log1p(
    count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
)
count_matrix_norm.index = peaks
count_matrix_norm.to_csv(norm_matrix_output, sep='\t', index=True, header=True)
# Get dummies for other covariates, if they exist
if len(samplesheet.columns) > 1:
    covars = samplesheet.drop(columns=[_comparison["time"]])
    cov_encoded = pd.get_dummies(covars, drop_first=True)
    assert cov_encoded.shape[0] == count_matrix_norm.shape[1]
else:
    cov_encoded = None
assert len(time) == count_matrix_norm.shape[1]

#Paramdic to fit GPs
fit_parameters = {
    'time': time,
    'time_grid': time_grid,
    'unique_times': unique_times,
    'cov_encoded': cov_encoded,
    'perms': PERMS,
    'fit_type': _comparison['time_type']
}

results = Parallel(n_jobs=THREADS)(
    delayed(fit_gp)(row.values, fit_parameters)
    for _, row in track(
        count_matrix_norm.iterrows(),
        total=count_matrix_norm.shape[0],
        description="Fitting GP"
    )
)

y_preds, y_stds, distances_list, lr_obs_list, pvals = zip(*results)

_index = count_matrix_norm.index
if fit_parameters['fit_type'] == 'ordinal':
    colnames = levels
else:
    colnames = time_grid.ravel()
y_pred_df = pd.DataFrame(np.vstack(y_preds), index=_index, columns=colnames)
y_std_df = pd.DataFrame(np.vstack(y_stds), index=_index, columns=colnames)
distances_df = pd.DataFrame(np.vstack(distances_list), index=_index)
if not distances_df.empty:
    distances_df.columns = [f"{a}-{b}" for a, b in zip(levels, levels[1:])]
results_df = pd.DataFrame({
    "lr_obs": lr_obs_list,
    "p_value": pvals
}, index=_index)
_, results_df["FDR"], _, _ = multipletests(results_df["p_value"], method="fdr_bh")

# Save results
results_df.to_csv(table_output, sep='\t', index=True, header=True)
y_pred_df.to_csv(table_output.replace('_gp_results.tsv', '_acc_pred.tsv'), sep='\t', index=True, header=True)
y_std_df.to_csv(table_output.replace('_gp_results.tsv', '_acc_pred_std.tsv'), sep='\t', index=True, header=True)
if not distances_df.empty:
    distances_df.to_csv(table_output.replace('_gp_results.tsv', '_distances.tsv'), sep='\t', index=True, header=True)
