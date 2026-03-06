import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

from aos.gp_fit import fit_peak_gp, fit_peak_gp_ordinal
from aos.gp_utils import build_results_df, run_parallel_gp
from aos.helper import get_elbow

# Avoid overthreading.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# Snakemake defined variables
mat = snakemake.input.mat
ss = snakemake.params.samplesheet
_comparison = snakemake.params.comparison
gp_timesteps = snakemake.params.gp_timesteps
threads = snakemake.threads
perms = snakemake.params.permutations
perm_cutoff = snakemake.params.permutation_cutoff
outputfolder = snakemake.params.outputfolder
sig_clustered = snakemake.params.sig_clustered
otable = snakemake.output.table

count_matrix = pd.read_csv(mat, sep="\t")
peaks = count_matrix[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
samplesheet = pd.read_csv(ss, sep="\t", index_col=0)

# Simple continuous fitting or ordinal fitting
match _comparison["time_type"]:
    case "continuous":
        count_matrix = count_matrix[list(samplesheet.index)]
        count_matrix_norm = np.log1p(
            count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
        )
        time = samplesheet[_comparison["time"]].values
        time_grid = np.linspace(time.min(), time.max(), gp_timesteps)[:, None]
        if len(samplesheet.columns) > 1:
            covars = samplesheet.drop(columns=[_comparison["time"]])
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
            threads,
            perms,
        )
        results_df = build_results_df(peaks, lrs, pvals, ypreds, ystds)
    case "ordinal":
        levels = _comparison["order"]
        samples_of_interest = [
            sam
            for sam in samplesheet.index
            if samplesheet.loc[sam, _comparison["time"]] in levels
        ]
        samplesheet = samplesheet.loc[samples_of_interest]
        count_matrix = count_matrix[samples_of_interest]
        count_matrix_norm = np.log1p(
            count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
        )
        timepoints = levels
        time_str = samplesheet[_comparison["time"]].values
        time_map = {t: i for i, t in enumerate(levels)}
        time = np.array([time_map[t] for t in time_str], dtype=float)

        unique_times = np.arange(len(levels), dtype=float)
        time_grid = unique_times[:, None]

        if len(samplesheet.columns) > 1:
            covars = samplesheet.drop(columns=[_comparison["time"]])
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
            threads,
            permutations=perms,
            extra_metric_names=["distances"],
        )
        results_df = build_results_df(peaks, lrs, pvals, ypreds, ystds)
        transition_labels = [
            f"{levels[i]}->{levels[i + 1]}" for i in range(len(levels) - 1)
        ]
        distances_df = pd.DataFrame(
            np.vstack(extras["distances"]),
            index=peaks,
            columns=transition_labels,
        )
        distances_df.to_csv(
            snakemake.params.outputfolder + "/transition_distances.tsv", sep="\t"
        )


# Cluster significant patterns, plot, and save results
if not results_df[results_df["FDR"] < perm_cutoff].empty:
    patterns = np.vstack(results_df.loc[results_df["FDR"] < perm_cutoff, "y_pred"])
    patterns_scaled = patterns - patterns.mean(axis=1, keepdims=True)
    patterns_scaled /= patterns.std(axis=1, keepdims=True) + 1e-8

    K_range = range(1, 20)
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
    fig, ax = plt.subplots(nrows=K, figsize=(8, 12), tight_layout=True)
    for k in range(K):
        ax[k].plot(time_grid, patterns_scaled[labels == k].T, alpha=0.3)
        ax[k].plot(
            time_grid,
            patterns_scaled[labels == k].mean(axis=0),
            color="black",
            linewidth=2,
        )
        ax[k].set_ylabel(f"Cluster {k}")
        if _comparison["time_type"] == "ordinal":
            ax[k].set_xticks(unique_times)
            ax[k].set_xticklabels(levels, rotation=90)

    fig.savefig(outputfolder + "/clustered_patterns.png", dpi=300)
    sig_results = results_df[results_df["FDR"] < perm_cutoff]
    sig_results["labels"] = labels
    sig_results.to_csv(sig_clustered, sep="\t")
else:
    print(f"No peaks meet cutoff threshold of {perm_cutoff}")
results_df.to_csv(otable, sep="\t")
