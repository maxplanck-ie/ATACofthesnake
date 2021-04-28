import os
import pandas as pd
import sys
import glob
import rich
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings

def checkExist(filesList):
    notExistent = []
    for fileName in filesList:
        if not os.path.exists(fileName):
            notExistent.append(fileName)
    if len(notExistent) > 0:
        rich.print("[red]File(s) not found![/red]")
        rich.print("Check your path for: [red]{}[/red]".format(notExistent))
    else:
        rich.print("Required input found. Moving on...")
        return True

def returnScriptPath():
    return os.path.dirname(__file__)

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


def plotter(what, inFiles, outFile, conds=None, outDir=None):
    colors = ["windows blue",
              "amber", "greyish",
              "faded green", "dusty purple"]
    sns.set_palette(sns.xkcd_palette(colors))
    if what == 'frip':
        res = []
        for sample in inFiles:
            with open(sample) as f:
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
        print("idxStat plot invoked.")
        res = []
        print("input files: {}".format(inFiles))
        for sample in inFiles:
            with open(sample) as f:
                chromcount = 0
                for line in f:
                    chrom = str(line.strip().split()[0])
                    count = int(line.strip().split()[1])
                    if chrom.lower().startswith('m') or \
                        'mito' in chrom.lower():
                        mtcount = count
                    else:
                        chromcount += count
                baseName = sample.split('/')[-1].strip('.idxstat.txt')
                res.append([baseName,
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


def returnCompfromSample(sample, ss, outDir):
    for comparison in ss['Comp']:
        if sample in ss[comparison]:
            return os.path.join(outDir, 'Peaks', comparison + "_peaks.bed")


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


def sortGTF(GTF, outDir):
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
    GTF.to_csv(outDir + '/genes.sort.gtf',
               header=False,
               index=False,
               sep='\t')

def GTFtoTSS(GTF, outDir):
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
    TSSdf.to_csv(outDir + '/TSS.bed', sep='\t', index=False, header=False)


def readBamDir(bamDir):
    bams = []
    if not os.path.exists(bamDir):
        print('Directory {} not found.'.format(bamDir))
        sys.exit()
    for i in os.listdir(bamDir):
        if i.endswith('bam'):
            bams.append(i.replace(".bam", ""))
    for bamFile in bams:
        if '-' in i:
            return "Rename bamFiles, no '-' allowed."
            sys.exit()
    if len(bams) == 0:
        print('No bamfiles found. Check if {} contains files with a bam extension.')
        sys.exit()
    return sorted(bams)

def returnPeaks(outDir, Comp, mergeStatus):
    print(Comp)
    peakDir = 'MACS2_mergeBAM' if mergeStatus else 'MACS2'
    return os.path.join(outDir, peakDir, Comp + "_peaks.bed")

def readss(ss, bams):
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
        for sample in ss.Sample:
            if sample.replace(".bam", "") not in bams:
                return "Can't find {}. Please correct.".format(sample)
                sys.exit()
        diffDic = {}
        diffDic['Comp'] = list(ss.Comp.unique())
        diffDic['CompCond'] = []
        diffDic['CompCondDic'] = {}
        # Populate samples in comparison.
        for comparison in ss.Comp.unique():
            diffDic[comparison] = list(sample.replace('.bam','') for sample in ss[ss['Comp'] == comparison].Sample)
            conds = list(ss[ss['Comp'] == comparison].Cond.unique())
            diffDic['CompCondDic'][comparison] = conds
            if len(conds) != 2:
                print("More then 1 condition within a comparison. Exiting..")
                sys.exit()
            else:
                tempdf = ss[ss['Comp'] == comparison]
                CondComp1 = comparison + '_' + conds[0]
                diffDic['CompCond'].append(CondComp1)
                diffDic[CondComp1] = list(sample.replace('.bam','') for sample in tempdf[tempdf['Cond'] == conds[0]].Sample)
                CondComp2 = comparison + '_' + conds[1]
                diffDic['CompCond'].append(CondComp2)
                diffDic[CondComp2] = list(sample.replace('.bam','') for sample in tempdf[tempdf['Cond'] == conds[1]].Sample)              

                # Include batchStatus:
                if batchStatus:
                    diffDic['Batch'] = {}
                    for row in range(len(ss)):
                        diffDic['Batch'][ss.Sample[row].replace('.bam','')] = ss.Batch[row]
        return diffDic
    else:
        print("Column headers not ok, (expected [Sample, Cond, Comp])")
        sys.exit()

def conditionsfromCount(countmat, ss, comp):
    # exception to return if countmat doesn't exist -> dryrun fails
    if not os.path.exists(countmat):
        return -1
    else:
        with open(countmat) as f:
            header = f.readline().strip().split()
            header = header[3:]
        # Get the two condition names for a comparison.
        conditions = ss['CompCondDic'][comp]
        CompCond1 = comp + '_' + conditions[0]
        CompCond2 = comp + '_' + conditions[1]
        flipDic = {}
        for sample in ss[CompCond1]:
            flipDic[sample] = CompCond1
        for sample in ss[CompCond2]:
            flipDic[sample] = CompCond2
        conditionOrder = []
        for sample in header:
            conditionOrder.append(flipDic[sample])
        return ','.join(conditionOrder)

def batchesfromCount(countmat, ss, comp):
    # exception to return if countmat doesn't exist -> dryrun fails
    if not os.path.exists(countmat):
        return -1
    else:
        with open(countmat) as f:
            header = f.readline().strip().split()
            header = header[3:]
        if 'Batch' in ss:
            batchOrder = []
            for sample in header:
                batchOrder.append(ss['Batch'][sample])
            return ','.join(batchOrder)
        else:
            return None