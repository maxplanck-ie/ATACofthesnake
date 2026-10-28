# lib imports
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests

from aos.gp_fitter import fit_gp_interaction, log_progress

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
ALPHA = snakemake.params.gp_alpha
INT_STRING = snakemake.params.interaction
# First timepoint is shared by all interaction levels (e.g. untreated baseline).
COMMON_BASELINE = _comparison.get("common_baseline", False)
table_output = snakemake.output.table

if _comparison["time_type"] == "ordinal":
    # Set up variables
    count_matrix = pd.read_csv(mat, sep="\t", dtype={"#chr": str})
    peaks = count_matrix[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
    samplesheet = pd.read_csv(ss, sep="\t", index_col=0)
    levels = _comparison["order"]

    # In case of ordinal from cont. data, match case in levels for sake of subsampling
    _dtype = samplesheet[_comparison["time"]].dtype
    _levels_case = [_dtype.type(i) for i in _comparison["order"]]

    samples_of_interest = [
        sam
        for sam in samplesheet.index
        if samplesheet.loc[sam, _comparison["time"]] in _levels_case
    ]
    samplesheet = samplesheet.loc[samples_of_interest]

    _n_int = samplesheet[INT_STRING].nunique()
    _keep = [
        samplesheet.loc[samplesheet[_comparison["time"]] == lvl, INT_STRING].nunique()
        == _n_int
        for lvl in _levels_case
    ]
    # A common baseline doesn't need every level.
    if COMMON_BASELINE:
        _keep[0] = True
    if not all(_keep):
        _msg = (
            f"Interaction {INT_STRING}: dropping timepoints "
            f"{[lvl for lvl, k in zip(levels, _keep) if not k]}, not observed "
            f"for every {INT_STRING} level."
        )
        if not _keep[0]:
            _msg += (
                f" If {levels[0]} is a baseline shared by all {INT_STRING} levels "
                "(e.g. untreated), set 'common_baseline: true' in the comparison."
            )
        print(_msg)
        levels = [lvl for lvl, k in zip(levels, _keep) if k]
        _levels_case = [lvl for lvl, k in zip(_levels_case, _keep) if k]
        samplesheet = samplesheet[samplesheet[_comparison["time"]].isin(_levels_case)]
    if len(levels) < 2:
        raise ValueError(
            f"Interaction {INT_STRING}: fewer than 2 timepoints are observed for "
            f"every {INT_STRING} level ({levels}); cannot test the interaction."
        )
    count_matrix = count_matrix[list(samplesheet.index)]

    _time_str = samplesheet[_comparison["time"]].values.astype(str)
    _time_map = {t: i for i, t in enumerate(levels)}
    _time = np.array([_time_map[t] for t in _time_str], dtype=float)
    time = _time[:, None]
    unique_times = np.arange(len(levels), dtype=float)
    time_grid = unique_times[:, None]
    baseline = 0.0 if COMMON_BASELINE else None

else:
    assert _comparison["time_type"] == "continuous"
    # Set up variables
    count_matrix = pd.read_csv(mat, sep="\t", dtype={"#chr": str})
    peaks = count_matrix[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
    samplesheet = pd.read_csv(ss, sep="\t", index_col=0)

    count_matrix = count_matrix[list(samplesheet.index)]

    time = samplesheet[_comparison["time"]].values[:, None]
    time_grid = np.linspace(time.min(), time.max(), GP_TIMESTEPS)[:, None]
    # Just needed in ordinal case, not here.
    unique_times = sorted(samplesheet[_comparison["time"]].unique())
    # time_grid starts at time.min(), so predictions there are masked too.
    baseline = float(time.min()) if COMMON_BASELINE else None

# Common prep

count_matrix_norm = np.log1p(count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6)
count_matrix_norm.index = peaks

# Check for nuisance variables to control for
if len(samplesheet.columns) > 2:
    covars = samplesheet.drop(columns=[_comparison["time"], INT_STRING])
    nuis_encoded = pd.get_dummies(covars, drop_first=True)
    if nuis_encoded.shape[0] != count_matrix_norm.shape[1]:
        raise ValueError(
            f"Nuisance covariate table has {nuis_encoded.shape[0]} samples but "
            f"the count matrix has {count_matrix_norm.shape[1]}; samplesheet "
            "and count matrix are misaligned."
        )
    # Covariates can become constant after dropping timepoints.
    if nuis_encoded.shape[1] == 0:
        nuis_encoded = None
else:
    nuis_encoded = None
# Encode interaction
int_encoded = pd.get_dummies(samplesheet[INT_STRING], drop_first=True)
if len(time) != count_matrix_norm.shape[1]:
    raise ValueError(
        f"Time vector has {len(time)} samples but the count matrix has "
        f"{count_matrix_norm.shape[1]}; samplesheet and count matrix are "
        "misaligned."
    )

# Get levels for prediction.
_all_levels = sorted(samplesheet[INT_STRING].unique())
_ref = [level for level in _all_levels if level not in int_encoded.columns]

fit_parameters = {
    "time": time,
    "time_grid": time_grid,
    "unique_times": unique_times,
    "nuis_encoded": nuis_encoded,
    "int_encoded": int_encoded,
    "int_labels": _ref + list(int_encoded.columns),
    "perms": PERMS,
    "alpha": ALPHA,
    "fit_type": _comparison["time_type"],
    "baseline": baseline,
}

results = list(
    log_progress(
        Parallel(n_jobs=THREADS, return_as="generator")(
            delayed(fit_gp_interaction)(row.values, fit_parameters)
            for _, row in count_matrix_norm.iterrows()
        ),
        total=count_matrix_norm.shape[0],
        label=f"GP interaction {snakemake.wildcards.comparison} - {INT_STRING}",
        logfile=snakemake.log[0],
    )
)

y_preds, y_stds, distances_list, lr_obs_list, pvals, lr_shape_list, p_shape_list = zip(
    *results
)

_index = count_matrix_norm.index
int_labels = fit_parameters["int_labels"]
n_levels = len(int_labels)

if fit_parameters["fit_type"] == "ordinal":
    colnames = levels
else:
    colnames = time_grid.ravel()

midx = pd.MultiIndex.from_arrays(
    [np.repeat(_index, n_levels), np.tile(int_labels, len(_index))],
    names=["gene", "interaction"],
)

y_pred_df = pd.DataFrame(np.vstack(y_preds), index=midx, columns=colnames)
y_std_df = pd.DataFrame(np.vstack(y_stds), index=midx, columns=colnames)

distances_df = pd.DataFrame(np.vstack(distances_list), index=_index)
if not distances_df.empty:
    distances_df.columns = [f"{a}-{b}" for a, b in zip(levels, levels[1:])]

# lr_obs/p_value/FDR: any difference between levels (offset and/or time-varying).
# *_shape: time-varying differences only, adjusted stage-wise in gp_postprocess.
results_df = pd.DataFrame(
    {
        "lr_obs": lr_obs_list,
        "p_value": pvals,
        "lr_shape": lr_shape_list,
        "p_shape": p_shape_list,
    },
    index=_index,
)
_, results_df["FDR"], _, _ = multipletests(results_df["p_value"], method="fdr_bh")

# Save results
results_df.to_csv(table_output, sep="\t", index=True, header=True)
y_pred_df.to_csv(
    table_output.replace("_gp_results.tsv", "_acc_pred.tsv"),
    sep="\t",
    index=True,
    header=True,
)
y_std_df.to_csv(
    table_output.replace("_gp_results.tsv", "_acc_pred_std.tsv"),
    sep="\t",
    index=True,
    header=True,
)
if not distances_df.empty:
    distances_df.to_csv(
        table_output.replace("_gp_results.tsv", "_distances.tsv"),
        sep="\t",
        index=True,
        header=True,
    )
