import matplotlib.pyplot as plt
from gimmemotifs.motif import read_motifs
from pathlib import Path
import pandas as pd
import numpy as np
import sys

# params
comp = snakemake.params.motifcomp
motdir = Path(f'motifs/{comp}/')
# I/O
motiffile = snakemake.input.motiffile
gimmemotiffile = motdir / 'clusteredmotifs_consensus_motifs.pwm'
opng = motdir / 'motif_enrichment.png'

amefiles = [p / 'ame.tsv' for p in motdir.iterdir() if p.is_dir() and "_ame" in p.name]
fdr_cutoff = snakemake.params.fdr_cutoff

with open(motiffile,'r') as f:
    with open(gimmemotiffile, 'w') as o:
        ret = False
        for line in f:
            if line.startswith('MOTIF'):
                name = '>' + line.strip().split(' ')[1]
                o.write(name + '\n')
                ret = True
            elif line.strip() and ret and not line.startswith('letter'):
                o.write(line.strip().replace('  ', '\t') + '\n')

def parse_ame(amefiles, fdr_cutoff):
    dfs = []
    groups = []
    for ame in amefiles:
        try:
            a = pd.read_table(ame, sep='\t', comment='#')
        except pd.errors.EmptyDataError:
            continue
        _grp = ame.parts[-2].replace('_ame', '')
        a['group'] = _grp
        groups.append(_grp)
        a['score'] = np.log2( (a['%TP'] + 1e-6) / (a['%FP'] + 1e-6) )
        a = a[a['adj_p-value'] < fdr_cutoff][['motif_ID', 'motif_alt_ID', 'adj_p-value', '%TP', '%FP', 'group', 'score']]
        dfs.append(a)
    df = pd.concat(dfs)
    # Fill 'missing ones' with 0
    df = df.pivot(index=["motif_ID", "motif_alt_ID"], columns="group", values="score").reset_index().fillna(0)
    # Get motif order based on scores
    df = df.sort_values(by=[c for c in df.columns if c not in ['motif_ID', 'motif_alt_ID']], ascending=False)
    motif_order = df['motif_ID'].tolist()
    df = df.melt(id_vars=['motif_ID', 'motif_alt_ID'])
    df.columns = ['motif_ID', 'motif_alt_ID', 'group', 'score']
    df['group'] = pd.Categorical(df["group"], categories=groups, ordered=True)
    return df, motif_order

ame_results, plotmotifs = parse_ame(amefiles, fdr_cutoff)
if ame_results.empty:
    sys.exit(0)

pfms = dict([(m.id, m) for m in read_motifs(str(gimmemotiffile)) if m.id in plotmotifs])
motifs = [m for m in pfms.keys()]
assert len(pfms) == len(plotmotifs), "not all motifs in enrichment results found back in original meme file."

group_order = ame_results['group'].cat.categories.tolist()

fig, axs = plt.subplots(figsize=(8,10), constrained_layout=True)
axs.set_axis_off()
gs = fig.add_gridspec(len(motifs),4)

# Plot motifs
for i,mot in enumerate(plotmotifs):
    # Visualize motif
    _x = fig.add_subplot(gs[i, 0:2])
    pfms[mot].plot_logo(ax=_x, title=False, ylabel=True)
    _x.set_axis_off()

# Visualize scores
ax_hm = fig.add_subplot(gs[0:, 2:4], frameon=False)

ame_pivot = ame_results.pivot(index='motif_ID', columns='group', values='score').loc[plotmotifs]
cax = ax_hm.imshow(ame_pivot.values, aspect='auto', cmap='Reds', origin='upper')
name_map = dict(zip(ame_results["motif_ID"], ame_results["motif_alt_ID"]))
ax_hm.set_xticks(np.arange(len(group_order)))
ax_hm.set_xticklabels(group_order, rotation=90)
ax_hm.set_yticks(np.arange(len(plotmotifs)))
ax_hm.set_yticklabels([name_map[i] for i in plotmotifs])
cbar = fig.colorbar(cax, ax=ax_hm)
cbar.set_label('log2(%TP/%FP)')
fig.savefig(opng, dpi=300)
