import os
import pandas as pd
import sys
import rich
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings


def checkNumDiff(paramDic):
    for Comp in paramDic['Comp']:
        compFolder = "AOS/diffAcc_" + str(Comp)
        countUp = 0
        countDown = 0
        with open(os.path.join(compFolder,
                  str(Comp) + '_edgeR_annotated_UP.tsv')) as f:
            for line in f:
                countUp += 1
        with open(os.path.join(compFolder,
                  str(Comp) + '_edgeR_annotated_DOWN.tsv')) as f:
            for line in f:
                countDown += 1
        if countUp > 11 and countUp > 11:
            if 'diffComp' not in paramDic:
                paramDic['diffComp'] = {}
            paramDic['diffComp'][Comp] = paramDic['Comp'][Comp]
    return paramDic


def setCeil(pdSer):
    maxSer = max(pdSer)
    maxSerScale = float((maxSer*10/10))
    return np.ceil(maxSerScale)


def plotter(what, inFiles, outFile, conds=None):
    colors = ["windows blue",
              "amber", "greyish",
              "faded green", "dusty purple"]
    sns.set_palette(sns.xkcd_palette(colors))
    if what == 'frip':
        res = []
        for sample in inFiles:
            for i in sample:  # snakemake returns a nested list.
                filestr = 'AOS/QC/' + str(i) + '.FRiP.txt'
                with open(filestr) as f:
                    for line in f:
                        if line.strip().split()[0] != 'sample':
                            sample = str(line.strip().split()[0])
                            frip = float(line.strip().split()[2])
                            gcov = float(line.strip().split()[3])
                            res.append([sample, frip, gcov])
        df = pd.DataFrame(res)
        df.columns = ['sample', 'frip', 'peak_genome_coverage']
        df = pd.melt(df, id_vars='sample')
        fig, ax1 = plt.subplots()
        g = sns.barplot(x=df['sample'], y=df['value'],
                        hue=df['variable'], data=df, ax=ax1)
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
        g = ax1.set_ylabel('Frip score')
        frip = df[df.variable == 'frip']['value']
        g = ax1.set_ylim((0,
                          0.2 + setCeil(frip)))
        g = ax2 = ax1.twinx()
        peakGCov = df[df.variable == 'peak_genome_coverage']['value']
        g = ax2.set_ylim((0,
                          setCeil(peakGCov)))
        g = ax2.set_ylabel('Peak genome coverage')
        plt.tight_layout()
        g.figure.savefig(outFile, dpi=300)
    if what == 'idxstat':
        res = []
        for sample in inFiles:
            for i in sample:
                filestr = 'AOS/QC/' + str(i) + '.idxstat.txt'
                with open(filestr) as f:
                    chromcount = 0
                    for line in f:
                        chrom = str(line.strip().split()[0])
                        count = int(line.strip().split()[1])
                        if chrom.lower().startswith('m') or \
                           'mito' in chrom.lower():
                            mtcount = count
                        else:
                            chromcount += count
                    res.append([i,
                               mtcount/(chromcount+mtcount),
                               chromcount/(chromcount+mtcount)])
        df = pd.DataFrame(res)
        df.columns = ['sample', 'MTfrac', 'Chromfrac']
        df = pd.melt(df, id_vars='sample')
        g = sns.barplot(x=df['sample'],
                        y=df['value'],
                        hue=df['variable'], data=df)
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
        g.set_ylabel('% of Alignments')
        plt.tight_layout()
        g.figure.savefig(outFile, dpi=300)
    if what == 'maPlot':
        deDF = pd.read_csv(inFiles, sep='\t', index_col=False)
        # Fetch a string containing the condition, and how many regions
        upStr = conds[0][1] + \
            '_Open n=' + \
            str(len(deDF.loc[(deDF['FDR'] < 0.05) &
                (deDF['logFC'] > 0), ]))
        downStr = conds[0][0] + \
            '_Open n=' + \
            str(len(deDF.loc[(deDF['FDR'] < 0.05) &
                (deDF['logFC'] < 0), ]))
        # Define status column and fill conditionaly
        deDF['Status'] = 'nonSign n=' + \
                         str(len(deDF.loc[deDF['FDR'] > 0.05, ]))
        deDF.loc[(deDF['FDR'] < 0.05) &
                 (deDF['logFC'] > 0), 'Status'] = upStr
        deDF.loc[(deDF['FDR'] < 0.05) &
                 (deDF['logFC'] < 0), 'Status'] = downStr
        # Plot and save
        g = sns.scatterplot(data=deDF,
                            x='logCPM',
                            y='logFC',
                            alpha=0.4,
                            hue='Status')
        g.figure.savefig(outFile, dpi=300)


def diffCount(Complist):
    retainComp = []
    for Comp in Complist:
        counter = 0
        diffFile = "diffAcc_" + Comp + "/" + Comp + "_DESeq2_annotated.tsv"
        if os.path.exists(diffFile):
            with open(diffFile) as f:
                for line in f:
                    counter += 1
            if counter > 50:
                retainComp.append(Comp)
        else:
            retainComp.append(Comp)
    return retainComp


def returnCompfromSample(sample, paramDic):
    compList = []
    for comp in paramDic['Comp']:
        if sample in paramDic['Comp'][comp]['Samples']:
            compList.append(
                "AOS/MACS2/{}_Merged_peaks.narrowPeak".format(comp))
    if len(compList) == 1:
        return compList[0]


def summitIncorp(finalhits, summits, output):
    hits = pd.read_csv(finalhits, sep='\t', header=0)
    summits = pd.read_csv(summits, sep='\t', header=None)
    summits.columns = ['peak_chr',
                       'summit_start',
                       'summit_stop',
                       'peak_id',
                       'abs_peaksummit']
    summits['summit_start'] = summits['summit_start'] - 100
    summits['summit_stop'] = summits['summit_stop'] + 100
    del summits['peak_chr']
    mergeDF = pd.merge(hits, summits, on='peak_id')
    mergeDF.to_csv(output, header=True, index=True, sep='\t')


def mergeDiff_Ann(annotation, diff, outName):
    diffDF = pd.read_csv(diff,
                         sep='\t',
                         header=0,
                         index_col=False,
                         dtype={0: 'str'})
    diffDF = diffDF.set_index('peak_id')
    annDF = pd.read_csv(annotation,
                        sep='\t',
                        header=0,
                        dtype={0: 'str'},
                        index_col=None)
    annDF.index = annDF[['peak_chr',
                         'peak_start',
                         'peak_end']].apply(
                         lambda row: '_'.join(row.values.astype(str)), axis=1)
    res = pd.merge(diffDF,
                   annDF,
                   left_index=True,
                   right_index=True)
    if len(diffDF.index) != len(res.index):
        print("n(merged) != n(input). Double check peak annotations..")
    res.to_csv(outName, header=True, index=True, sep='\t')


def sortGTF(GTF):
    GTF = pd.read_csv(GTF,
                      sep='\t',
                      comment='#',
                      header=None,
                      dtype={0: 'str'})
    GTF.columns = ['chr', 'source', 'feature',
                   'start', 'end', 'score',
                   'strand', 'frame', 'attribute']
    GTF = GTF[GTF['feature'] == 'transcript']
    GTF = GTF.sort_values(["chr", "start"],
                          ascending=(True, True))
    GTF.to_csv('genes.sort.gtf',
               header=False,
               index=False,
               sep='\t')


def conditionsfromCount(countmat, paramDic):
    # exception to return if countmat doesn't exist -> dryrun fails
    if not os.path.exists(countmat):
        return -1
    else:
        with open(countmat) as f:
            header = f.readline().strip().split()
            header = header[3:]
        conditionOrder = []
        flipDic = {}
        for cond in paramDic:
            for sample in paramDic[cond]:
                flipDic[sample] = cond
        for sample in header:
            print(sample)
            conditionOrder.append(flipDic[sample])
        return ','.join(conditionOrder)

def batchesfromCount(countmat, paramDic):
    # exception to return if countmat doesn't exist -> dryrun fails
    if not os.path.exists(countmat):
        return -1
    else:
        with open(countmat) as f:
            header = f.readline().strip().split()
            header = header[3:]
        if paramDic['batchStatus'] == 1:
            flipDic = {}
            for i in range(len(paramDic['Samples'])):
                flipDic[paramDic['Samples'][i]] = paramDic['Batch'][i]
            batchOrder = []
            for sample in header:
                print(sample)
                batchOrder.append(flipDic[sample])
            return ','.join(batchOrder)


def GTFtoTSS(GTF):
    TSS = []
    linecount = 0
    with open(GTF) as f:
        for line in f:
            linecount += 1
    with open(GTF) as f:
        for line in f:
            if not line.startswith('#'):
                liLis = line.strip().split()
                if liLis[2] == 'transcript':
                    if liLis[6] == '+':
                        TSS.append([liLis[0],
                                    int(liLis[3]),
                                    int(liLis[3]) + 1])
                    elif liLis[6] == '-':
                        TSS.append([liLis[0],
                                    int(liLis[4]) - 1,
                                    int(liLis[4])])
    TSSdf = pd.DataFrame(TSS)
    TSSdf = TSSdf.drop_duplicates()
    TSSdf.to_csv('TSS.bed', sep='\t', index=False, header=False)


def readBamDir(bamDir):
    bams = []
    for i in os.listdir(bamDir):
        if i.endswith('bam'):
            bams.append(i.replace(".bam", ""))
    for bamFile in bams:
        if '-' in i:
            return "Rename bamFiles, no - allowed."
            sys.exit()
    return bams


def setdefault_readss(ss, bams):
    warnings.simplefilter(action='ignore', category=FutureWarning)
    ss = pd.read_csv(ss, sep='\t', header=0)
    if 'Sample' in list(ss.columns) and \
        'Cond' in list(ss.columns) and \
        'Comp' in list(ss.columns):
        batchStatus = 0
        if len(ss.columns) > 3:
            if 'Batch' in list(ss.columns):
                print('Batch column found.')
                batchStatus = 1
            else:
                return "More then three columns, one of them is not Batch."
                sys.exit()
        elif len(ss.Cond.unique()) != 2:
            return "Only two conditions allowed in sampleSheet.."
            sys.exit()
        for sample in ss.Sample:
            if sample.replace(".bam", "") not in bams:
                return "Can't find {}. Please correct.".format(sample)
                sys.exit()
        diffDic = {}
        diffDic["Cond"] = list(ss.Cond.unique())
        diffDic["Comp"] = {}
        diffDic["Samples"] = list(ss.Sample.str.replace(".bam", ""))
        if batchStatus == 1:
            diffDic["Batch"] = list(ss.Batch)
        for Comp in ss.Comp.unique():
            tempss = ss[ss['Comp'] == Comp]
            diffDic["Comp"][Comp] = {}
            key1 = tempss.Cond.unique()[0]
            key2 = tempss.Cond.unique()[1]
            diffDic["Comp"][Comp]['Cond'] = {}
            diffDic["Comp"][Comp]['Cond'][key1] = \
                list(
                tempss[tempss['Cond'] == key1]['Sample'].str.replace(
                    ".bam", ""))
            diffDic["Comp"][Comp]['Cond'][key2] = \
                list(
                tempss[tempss['Cond'] == key2]['Sample'].str.replace(
                    ".bam", ""))
            samplesList = \
                diffDic["Comp"][Comp]['Cond'][key1] + \
                diffDic["Comp"][Comp]['Cond'][key2]
            diffDic["Comp"][Comp]['Samples'] = samplesList
        for Comp in ss.Comp.unique():
            table = \
                rich.table.Table(
                    title="Samples - Conditions: {}".format(Comp))
            table.add_column(key1, justify="center", style="cyan")
            table.add_column(key2, justify="center", style="green")
            for i in range(max(
                        [len(diffDic["Comp"][Comp]['Cond'][key1]),
                        len(diffDic["Comp"][Comp]['Cond'][key2])])):
                try:
                    table.add_row(diffDic["Comp"][Comp]['Cond'][key1][i],
                                    diffDic["Comp"][Comp]['Cond'][key2][i])
                except Exception:
                    try:
                        table.add_row(
                            diffDic["Comp"][Comp]['Cond'][key1][i], "")
                    except Exception:
                        table.add_row(
                            "", diffDic["Comp"][Comp]['Cond'][key2][i])
            console = rich.console.Console()
            console.print(table)
        diffDic['baseDir'] = os.path.dirname(__file__)
        return diffDic, batchStatus
    else:
        return "Column headers not ok, (expected [Sample, Cond, Comp])"
        sys.exit()
