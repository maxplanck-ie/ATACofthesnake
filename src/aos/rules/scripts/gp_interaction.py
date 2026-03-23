# lib imports
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import os
from aos.gp_fitter import fit_gp_interaction
from statsmodels.stats.multitest import multipletests


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

INT_STRING = 'treatment'

print(f"Interaction GP for comparison: {_comparison}")

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
    unique_times = sorted(samplesheet[_comparison["time"]].unique())

# Common prep

count_matrix_norm = np.log1p(
    count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
)
count_matrix_norm.index = peaks

# Check for nuisance variables to control for
if len(samplesheet.columns) > 2:
    covars = samplesheet.drop(columns=[_comparison["time"], INT_STRING])
    nuis_encoded = pd.get_dummies(covars, drop_first=True)
    assert nuis_encoded.shape[0] == count_matrix_norm.shape[1]
else:
    nuis_encoded = None
# Encode interaction
int_encoded = pd.get_dummies(samplesheet[INT_STRING], drop_first=True)
assert len(time) == count_matrix_norm.shape[1]

# Get levels for prediction.
_all_levels = sorted(samplesheet[INT_STRING].unique())
_ref = [level for level in _all_levels if level not in int_encoded.columns]

fit_parameters = {
    'time': time,
    'time_grid': time_grid,
    'unique_times': unique_times,
    'nuis_encoded': nuis_encoded,
    'int_encoded': int_encoded,
    'int_labels':  _ref + list(int_encoded.columns),
    'perms': PERMS,
    'fit_type': _comparison['time_type']
}

results = Parallel(n_jobs=THREADS)(
    delayed(fit_gp_interaction)(row.values, fit_parameters)
    for _, row in tqdm(count_matrix_norm.iterrows(), total=count_matrix_norm.shape[0])
)

y_preds, y_stds, distances_list, lr_obs_list, pvals = zip(*results)

_index      = count_matrix_norm.index
int_labels  = fit_parameters["int_labels"]
n_levels    = len(int_labels)

if fit_parameters['fit_type'] == 'ordinal':
    colnames = levels
else:
    colnames = time_grid.ravel()

midx = pd.MultiIndex.from_arrays(
    [np.repeat(_index, n_levels), np.tile(int_labels, len(_index))],
    names=["gene", "interaction"]
)

y_pred_df = pd.DataFrame(np.vstack(y_preds), index=midx, columns=colnames)
y_std_df  = pd.DataFrame(np.vstack(y_stds),  index=midx, columns=colnames)

distances_df = pd.DataFrame(np.vstack(distances_list), index=_index)
if not distances_df.empty:
    distances_df.columns = [f"{a}-{b}" for a, b in zip(levels, levels[1:])]

results_df = pd.DataFrame({
    "lr_obs":   lr_obs_list,
    "p_value":  pvals,
}, index=_index)
_, results_df["FDR"], _, _ = multipletests(results_df["pvalue"], method="fdr_bh")

# Save results
results_df.to_csv(table_output, sep='\t', index=True, header=True)
y_pred_df.to_csv(table_output.replace('_gp_results.tsv', '_acc_pred.tsv'), sep='\t', index=True, header=True)
y_std_df.to_csv(table_output.replace('_gp_results.tsv', '_acc_pred_std.tsv'), sep='\t', index=True, header=True)
if not distances_df.empty:
    distances_df.to_csv(table_output.replace('_gp_results.tsv', '_distances.tsv'), sep='\t', index=True, header=True)
