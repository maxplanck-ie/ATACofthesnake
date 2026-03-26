import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import os
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# parameters
#config['cutoffs']['permutation_cutoff']
perm_cutoff = snakemake.params.permutation_cutoff
comp_name = snakemake.params.comp_name
min_sigpeaks = snakemake.params.min_sigpeaks

# I
results = pd.read_table(snakemake.input.results, sep='\t', index_col=0)
y_pred = pd.read_table(snakemake.params.y_pred, sep='\t', index_col=0)

odir = snakemake.params.odir
# O
if getattr(snakemake.params, "int", None):
    k_table = Path(odir) / f"inttest_{comp_name}_{snakemake.params.int}_k_table.tsv"
    k_plot = Path(odir) / f"inttest_{comp_name}_{snakemake.params.int}_k_plot.png"
else:
    k_table = Path(odir) / f"{comp_name}_k_table.tsv"
    k_plot = Path(odir) / f"{comp_name}_k_plot.png"

# In case output exists already (params, not I/O in snakemake), delete them.
for path in [k_table, k_plot]:
    if path.exists():
        path.unlink()

sig = results[results['FDR'] < perm_cutoff]
if len(sig) < min_sigpeaks:
    print(f"Only {len(sig)} significant peaks found for {comp_name} with permutation cutoff {perm_cutoff}. Need at least {min_sigpeaks} to continue.")
    sys.exit(0)

y_pred = y_pred.loc[sig.index]

# Avoid overthreading.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

def get_elbow(inertias, K_range):
    x1, y1 = K_range[0], inertias[0]
    x2, y2 = K_range[-1], inertias[-1]

    distances = []
    for x0, y0 in zip(K_range, inertias):
        num = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1)
        den = np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)
        distances.append(num / den)

    distances = np.array(distances)
    elbow_idx = np.argmax(distances)
    return K_range[elbow_idx]


def compute_inertia(k):
    km = KMeans(n_clusters=k, n_init=20, random_state=42)
    km.fit(patterns_scaled)
    return km.inertia_

if 'interaction' in y_pred.columns:
    _valcols = y_pred.columns.drop("interaction")
    patterns = np.vstack(
        y_pred.groupby(y_pred.index)[_valcols]
          .apply(lambda x: np.hstack(x.to_numpy()))
          .values
    )
else:
    patterns = np.vstack(y_pred.values)
patterns_scaled = patterns - patterns.mean(axis=1, keepdims=True)
patterns_scaled /= patterns.std(axis=1, keepdims=True) + 1e-8
K_range = range(2, min(20, len(patterns)))
inertias = Parallel(n_jobs=snakemake.threads)(
    delayed(compute_inertia)(k) for k in K_range
)

K_opt = get_elbow(inertias, K_range)
print(f"Selected K = {K_opt}")

K = K_opt
kmeans = KMeans(n_clusters=K, n_init=50, random_state=1337)
labels = kmeans.fit_predict(patterns_scaled)
sig['k'] = labels
sig.to_csv(k_table, sep='\t', index=True, header=True)

if 'interaction' in y_pred.columns:
    _valcols
    n_timepoints = len(_valcols)
    n_interactions = y_pred.groupby(y_pred.index).size().iloc[0]
    colors = plt.cm.tab10(range(n_interactions))
    def unstack_pattern(row):
        return row.reshape(n_interactions, n_timepoints)
    interaction_labels = (
        y_pred.groupby(y_pred.index)["interaction"]
              .apply(list)
              .iloc[0]
    )
    
    fig, ax = plt.subplots(nrows=K, figsize=(8, 12), tight_layout=True)
    
    for k in range(K):
        cluster_patterns = patterns_scaled[labels == k]
        for pattern in cluster_patterns:
            reshaped = pattern.reshape(n_interactions, n_timepoints)
            for i in range(n_interactions):
                ax[k].plot(
                    _valcols,
                    reshaped[i],
                    color="gray",
                    alpha=0.3
                )
        mean_pattern = cluster_patterns.mean(axis=0).reshape(n_interactions, n_timepoints)
    
        for i, name in enumerate(interaction_labels):
            ax[k].plot(
                _valcols,
                mean_pattern[i],
                color=colors[i],
                linewidth=2,
                label=name if k == 0 else None
            )
    
        ax[k].set_ylabel(f"Cluster {k}")
    ax[0].set_title(comp_name)
    ax[0].legend()
    fig.savefig(k_plot, dpi=300)
else:
    fig, ax = plt.subplots(nrows=K, figsize=(8, 12), tight_layout=True)
    for k in range(K):
        ax[k].plot(y_pred.columns, patterns_scaled[labels == k].T, alpha=0.3)
        ax[k].plot(
            y_pred.columns,
            patterns_scaled[labels == k].mean(axis=0),
            color="black",
            linewidth=2,
        )
        ax[k].set_ylabel(f"Cluster {k}")
    ax[0].set_title(comp_name)
    fig.savefig(k_plot, dpi=300)