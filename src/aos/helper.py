import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import repeat
import shutil
import numpy as np

def idx_to_mit(_f):
    count = 0
    totalcount = 0
    with open(_f) as f:
        for line in f:
            _l = line.strip().split()
            totalcount += int(_l[1])
            if 'mt' in _l[0].lower():
                count = int(_l[1])
    return(round(count/totalcount, 2))

def plotfragsize(fs, of):
    df = pd.read_csv(
        fs,
        sep='\t',
        index_col=None,
        comment='#'
    )
    df['Sample'] = [
        i.replace('.bam', '').replace('input/', '') for i in list(df['Sample'])
    ]
    reps = []
    for i,r in df.iterrows():
        size = int(r['Size'])
        occ = int(r['Occurrences'])
        if size < 1000:
            for k in repeat(size, occ):
                reps.append([r['Sample'], size])
    df = pd.DataFrame(reps)
    df.columns = ['sample', 'size']
    df.sort_values(by=['sample'], inplace=True)
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    d = sns.FacetGrid(df, row="sample", hue='sample',height=1.1, aspect=8, palette=pal)
    d.set(xlim=(0, 1000))
    # Densities.
    d.map(sns.kdeplot, "size",
        bw_adjust=.5, clip_on=False,
        fill=True, alpha=0.8, linewidth=1.5)
    # Nucleosome lines.
    d.map(plt.axvline, x=147, ls='--', c='red', alpha=0.3)
    d.map(plt.axvline, x=147*2, ls='--', c='red', alpha=0.3)
    d.map(plt.axvline, x=147*3, ls='--', c='red', alpha=0.3)
    # Clean up plots.
    d.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="right", va="center", transform=ax.transAxes)
    d.map(label, "sample")
    d.set_titles("")
    d.set(yticks=[], ylabel="", xlabel='fragment size (bps)')
    d.despine(bottom=True, left=True)
    # Overlap plots.
    d.figure.subplots_adjust(hspace=-0.4)
    # save fig.
    d.figure.savefig(of, dpi=300, bbox_inches='tight')

def plotfrip(fs, of):
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    df = pd.read_csv(
        fs,
        sep='\t',
        index_col=None,
        header=None,
        comment='#'
    )
    df.sort_values(by=[1], inplace=True, ascending=False)
    g = sns.barplot(
        data=df,
        x=0,
        y=1,
        color=pal[5]
    )
    g.tick_params(axis='x', labelrotation=90)
    g.set(xlabel='', ylabel='frip score')
    g.figure.savefig(of, dpi=300, bbox_inches='tight')

def plotixs(ixs, mitostring):
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    # ixs
    df = pd.read_csv(
        ixs,
        sep='\t',
        header=0,
        index_col=0,
        comment='#'
    )
    df = pd.DataFrame(round(df.loc[mitostring]/df.sum(), 2))
    df['sample'] = df.index
    df.sort_values(by=[0], inplace=True, ascending=False)
    g = sns.barplot(
        data=df,
        x='sample',
        y=0,
        color=pal[5]
    )
    g.tick_params(axis='x', labelrotation=90)
    g.set(xlabel='', ylabel='fraction mitochondrial reads')
    g.figure.savefig('figures/mitofraction.png', dpi=300, bbox_inches='tight')

def plotsieve(sieve):
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    df = pd.read_csv(
        sieve,
        sep='\t',
        header=0,
        index_col=0,
        comment='#'
    ).T
    df.sort_values(by=['fraction'], inplace=True, ascending=False)
    g = sns.barplot(
        data=df,
        x=df.index,
        y='fraction',
        color=pal[5]
    )
    g.tick_params(axis='x', labelrotation=90)
    g.set(xlabel='', ylabel='Fraction of surviving reads')
    g.figure.savefig('figures/alignmentsieve.png', dpi=300, bbox_inches='tight')

def maplot(tsv: str, of: str, compentry: dict):
    compentry.pop("type", None)
    compentry.pop("design", None)
    gr1 = list(compentry.keys())[0]
    gr2 = list(compentry.keys())[1]
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    df = pd.read_csv(tsv, sep='\t', header=0)
    sig_mask = df["FDR"] < 0.05
    df["sig"] = sig_mask.map({True: "S", False: "NS"})
    upc = int((sig_mask & (df["logFC"] > 0)).sum())
    downc = int((sig_mask & (df["logFC"] <= 0)).sum())
    g = sns.scatterplot(
        data=df,
        x='logCPM',
        y='logFC',
        edgecolor = None,
        s=5,
        hue='sig',
        alpha=0.5,
        palette={'NS':'#7f7f7f', 'S': pal[5]}
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.text(.99, .99, '{} - {} diff peaks.'.format(gr2, upc), ha='right', va='top', transform=g.transAxes)
    plt.text(.99, .01, '{} - {} diff peaks.'.format(gr1, downc), ha='right', va='bottom', transform=g.transAxes)
    g.figure.savefig(
        of,
        dpi=300,
        bbox_inches='tight'
    )

def merge_idx(i, o):
    _i = list(i)
    _o = str(o)
    print(list(_i))
    print(str(_o))
    rnames = ['samples']
    with open(_i[0]) as f:
      for line in f:
        contig = line.strip().split('\t')[0]
        if contig != '*':
          rnames.append(contig)
    dflists = [rnames]
    for sample in _i:
      slis = [sample.split('/')[1].replace('_ix.tsv', '')]
      with open(sample) as f:
        for line in f:
          if not line.startswith('*'):
            slis.append(int(line.strip().split()[1]))
      dflists.append(slis)
    df = pd.DataFrame(dflists)
    df.T.to_csv(_o, sep='\t', header=False, index=False)

def merge_sieve(i, o):
    _i = list(i)
    _o = str(o)
    data = [['sample', 'surviving', 'initial', 'fraction']]
    for sample in _i:
        with open(sample) as f:
            for line in f:
                if not line.startswith('#'):
                    lis = line.strip().split()
                    data.append(
                        [
                            lis[0].replace('input/', '').replace('.bam', ''),
                            int(lis[1]),
                            int(lis[2]),
                            int(lis[1])/int(lis[2])
                        ]
                    )
    df = pd.DataFrame(data)
    df.T.to_csv(_o , sep='\t', header=False, index=False)


def peak_boundaries(peaks, genomefa, peakset, of):
    if not peakset:
        chromdic = {}
        with open(genomefa) as f:
            header = None
            for line in f:
                if line.startswith('>'):
                    header = str(line.strip().replace('>', '').split(' ')[0])
                    chromdic[header] = 0
                else:
                    chromdic[header] += len(line.strip())
        bedlis = []
        peakchange = 0 
        with open(peaks) as f:
            for line in f:
                chrom = str(line.strip().split()[0])
                start = int(line.strip().split()[1])
                end = int(line.strip().split()[2])
                if end > chromdic[chrom]:
                    bedlis.append(
                        [chrom, start, chromdic[chrom]]
                    )
                    peakchange += 1
                else:
                    bedlis.append(
                        [chrom, start, end]
                    )
        print("Changed {} peaks.".format(peakchange))
        beddf = pd.DataFrame(bedlis)
        beddf.to_csv(
            of,
            sep='\t',
            header=False,
            index=False
        )
    else:
        shutil.copyfile(peakset, of)

def PCA_colors(samplesheet, samples):
    colors = [
        '#1f77b4',
        '#ff7f0e',
        '#2ca02c',
        '#d62728',
        '#9467bd',
        '#8c564b',
        '#e377c2',
        '#7f7f7f',
        '#bcbd22',
        '#17becf'
    ]
    if samplesheet:
        sdf = pd.read_csv(
            samplesheet,
            sep='\t',
            header=0
        )
        sdf = sdf.set_index('sample')
        sdf = sdf.loc[samples]
        colDic = {}
        colIx = 0
        for s in sdf.iloc[:,[0]].values:
            if s[0] not in colDic:
                colDic[s[0]] = colors[colIx]
                colIx += 1
        PCAstr = "--colors"
        for s in sdf.iloc[:,[0]].values:
            PCAstr += f" \"{colDic[s[0]]}\""
        return (PCAstr)
    return ("")

def get_elbow(inertias, K_range):
    x1, y1 = K_range[0], inertias[0]
    x2, y2 = K_range[-1], inertias[-1]

    distances = []
    for x0, y0 in zip(K_range, inertias):
        num = abs((y2 - y1)*x0 - (x2 - x1)*y0 + x2*y1 - y2*x1)
        den = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        distances.append(num / den)

    distances = np.array(distances)
    elbow_idx = np.argmax(distances)
    return K_range[elbow_idx]