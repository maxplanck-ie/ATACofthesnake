import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import repeat

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

def plotfragsize(frags):
    df = pd.read_csv(
        frags,
        sep='\t',
        index_col=None,
        comment='#'
    )
    df['Sample'] = [
        i.replace('.bam', '').replace('input/', '') for i in list(df['Sample'])
    ]
    reps = []
    for i,r in df.iterrows():
        size = int(r[0])
        occ = int(r[1])
        for k in repeat(size, occ):
            reps.append([r['Sample'], size])
    df = pd.DataFrame(reps)
    df.columns = ['sample', 'size']
    df.sort_values(by=['sample'], inplace=True)
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    d = sns.FacetGrid(df, row="sample", hue='sample',height=1.1, aspect=8, palette=pal)
    d.map(sns.kdeplot, "size",
        bw_adjust=.5, clip_on=False,
        fill=True, alpha=0.8, linewidth=1.5)
    d.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="right", va="center", transform=ax.transAxes)
    d.map(label, "sample")
    d.set_titles("")
    d.set(yticks=[], ylabel="", xlabel='fragment size (bps)')
    d.despine(bottom=True, left=True)
    d.figure.subplots_adjust(hspace=-0.4)
    d.figure.savefig('figures/fragmentsizes.png', dpi=300, bbox_inches='tight')

def plotfrip(frips):
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    df = pd.read_csv(
        frips,
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
    g.figure.savefig('figures/fripscores.png', dpi=300, bbox_inches='tight')

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

def maplot(tsv, of, gr1, gr2):
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    df = pd.read_csv(
        tsv,
        sep='\t',
        header=0
    )
    diflis = []
    upc = 0
    downc = 0
    for i,r in df.iterrows():
        if r['FDR'] < 0.05:
            diflis.append('S')
            if r['logFC'] > 0:
                upc += 1
            else:
                downc += 1
        else:
            diflis.append('NS')
    df['sig'] = diflis
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
    gr2df = pd.DataFrame(
        [i.split('_') for i in list(df[(df['logFC'] > 0) & (df['FDR'] < 0.05)]['peak_id'])]
    )
    gr1df = pd.DataFrame(
        [i.split('_') for i in list(df[(df['logFC'] < 0) & (df['FDR'] < 0.05)]['peak_id'])]
    )
    basefolder = os.path.abspath(tsv).split('/')[-2]

    gr2df.to_csv(
        os.path.join(basefolder, gr2 + '.bed'),
        sep='\t',
        header=False,
        index=False
    )
    gr1df.to_csv(
        os.path.join(basefolder, gr1 + '.bed'),
        sep='\t',
        header=False,
        index=False
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
