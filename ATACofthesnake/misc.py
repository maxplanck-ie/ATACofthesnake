import os
import pandas as pd
import sys
import rich
from rich.progress import Progress
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings
import collections

def checkNumDiff(paramDic):
    #diffDownstream.yaml
    diffDownstream = {}
    samples = []
    for Comp in paramDic['Comp']:
        compFolder = "diffAcc_" + str(Comp)
        countUp = 0
        with open(os.path.join(compFolder, str(Comp) + '_edgeR_annotated_UP.tsv')) as f:
            for line in f:
                countUp += 1
        with open(os.path.join(compFolder, str(Comp) + '_edgeR_annotated_DOWN.tsv')) as f:
            for line in f:
                countDown += 1
        if countUp > 11 and countUp > 11:
            paramDic['diffComp'][Comp] = paramDic['Comp']
    return paramDic


def plotter(what, inFiles, outFile, conds=None):
    colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
    sns.set_palette(sns.xkcd_palette(colors))
    if what == 'frip':
        res = []
        for sample in inFiles:
            for i in sample: #snakemake returns a nested list.
                filestr = 'QC/' + str(i) + '.FRiP.txt'
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
        g = sns.barplot(x=df['sample'], y=df['value'], hue=df['variable'],data=df, ax=ax1)
        g.set_xticklabels(g.get_xticklabels(),rotation=90)
        g = ax1.set_ylabel('Frip score')
        g = ax1.set_ylim((0,0.2+np.ceil(float(max(df[df.variable=='frip']['value']))*10)/10))
        g = ax2 = ax1.twinx()
        g = ax2.set_ylim((0,np.ceil(float(max(df[df.variable=='peak_genome_coverage']['value']))*10)/10))
        g = ax2.set_ylabel('Peak genome coverage')
        plt.tight_layout()
        g.figure.savefig(outFile, dpi=300)
    if what == 'idxstat':
        res = []
        for sample in inFiles:
            for i in sample:
                filestr = 'QC/' + str(i) + '.idxstat.txt'
                with open(filestr) as f:
                    chromcount = 0
                    for line in f:
                        chrom = str(line.strip().split()[0])
                        count = int(line.strip().split()[1])
                        if chrom.lower() == 'mt':
                            mtcount = count
                        else:
                            chromcount += count
                    res.append([i, mtcount/(chromcount+mtcount), chromcount/(chromcount+mtcount)])
        df = pd.DataFrame(res)
        df.columns = ['sample', 'MTfrac', 'Chromfrac']
        df = pd.melt(df, id_vars='sample')
        g = sns.barplot(x=df['sample'], y=df['value'], hue=df['variable'], data=df)
        g.set_xticklabels(g.get_xticklabels(),rotation=90)
        g.set_ylabel('% of Alignments')
        plt.tight_layout()
        g.figure.savefig(outFile, dpi=300)
    if what == 'maPlot':
        deDF = pd.read_csv(inFiles, sep='\t', index_col=False)
        #Fetch a string containing the condition, and how many regions
        upStr = conds[0][1] + '_Open n=' + str(len(deDF.loc[(deDF['FDR'] < 0.05) & (deDF['logFC'] > 0),])) # snakeMake returns nested list
        downStr = conds[0][0] + '_Open n=' + str(len(deDF.loc[(deDF['FDR'] < 0.05) & (deDF['logFC'] < 0),])) # ditto
        #Define status column and fill conditionaly
        deDF['Status'] = 'nonSign n=' + str(len(deDF.loc[deDF['FDR'] > 0.05,]))
        deDF.loc[(deDF['FDR'] < 0.05) & (deDF['logFC'] > 0), 'Status'] = upStr
        deDF.loc[(deDF['FDR'] < 0.05) & (deDF['logFC'] < 0), 'Status'] = downStr
        #Plot and save
        g = sns.scatterplot(data=deDF, x='logCPM', y='logFC', alpha=0.4, hue='Status')
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
            compList.append("MACS2/{}_Merged_peaks.narrowPeak".format(comp))
    if len(compList) == 1:
        return compList[0]

def summitIncorp(finalhits, summits, output):
    hits = pd.read_csv(finalhits, sep='\t', header=0)
    summits = pd.read_csv(summits, sep='\t', header=None)
    summits.columns = ['peak_chr', 'summit_start', 'summit_stop', 'peak_id', 'abs_peaksummit']
    summits['summit_start'] = summits['summit_start'] - 100
    summits['summit_stop'] = summits['summit_stop'] + 100
    del summits['peak_chr']
    mergeDF = pd.merge(hits, summits, on='peak_id')
    mergeDF.to_csv(output, header=True, index=True, sep='\t')

def mergeDiff_Ann(annotation, diff, outName):
    diffDF = pd.read_csv(diff, sep='\t', header=0, index_col=False, dtype={0:'str'})
    diffDF = diffDF.set_index('peak_id')
    annDF = pd.read_csv(annotation, sep='\t', header=0, dtype={0:'str'}, index_col=None)
    annDF.index = annDF[['peak_chr','peak_start','peak_end']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    res = pd.merge(diffDF, annDF, left_index=True, right_index=True)
    if len(diffDF.index) != len(res.index):
        print("Merged set contains a different number of peaks than the input set. Double check peak annotations..")
    res.to_csv(outName, header=True, index=True, sep='\t')

def sortGTF(GTF):
    GTF = pd.read_csv(GTF, sep='\t', comment='#', header=None, dtype={0:'str'})
    GTF.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand','frame','attribute']
    GTF = GTF[GTF['feature'] == 'gene']
    GTF = GTF.sort_values(["chr", "start"], ascending = (True, True))
    GTF.to_csv('genes.sort.gtf', header=False, index=False, sep='\t')


def conditionsfromCount(countmat, paramDic):
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


def GTFtoTSS(GTF):
    TSS = []
    linecount = 0
    with open(GTF) as f:
        for line in f:
            linecount += 1
    with rich.progress.Progress() as progress:
        task = progress.add_task("Extracting TSS from GTF", total=linecount)
        with open(GTF) as f:
            for line in f:
                if not line.startswith('#'):
                    if line.strip().split()[2] == 'gene':
                        geneid = line.strip().split()[8].split(';')[0].replace('gene_id ',"").replace('"', '')
                        if line.strip().split()[6] == '+':
                            TSS.append([line.strip().split()[0], int(line.strip().split()[3]), int(line.strip().split()[3]) + 1])
                        elif line.strip().split()[6] == '-':
                            TSS.append([line.strip().split()[0], int(line.strip().split()[4]) - 1, int(line.strip().split()[4])])
                progress.advance(task)
    TSSdf = pd.DataFrame(TSS)
    TSSdf.to_csv('TSS.bed', sep='\t', index=False, header=False)

def readBamDir(bamDir):
    bams = []
    for i in os.listdir(bamDir):
        if i.endswith('bam'):
            bams.append(i.replace(".bam",""))
    for bamFile in bams:
        if '-' in i:
            return "Illegal character '-' found in bamfiles. Rename them and try again."
            sys.exit()
    return bams

def setdefault_readss(ss, bams):
    warnings.simplefilter(action='ignore', category=FutureWarning)
    ss = pd.read_csv(ss, sep='\t', header=0)
    if ss.columns[0] == 'Sample' and ss.columns[1] == 'Cond' and ss.columns[2] == 'Comp':
        if len(ss.columns) > 3:
            return "Error parsing sampleSheet, more than 3 columns found."
            sys.exit()
        elif len(ss.Cond.unique()) != 2:
            return "Error parsing sampleSheet, I need exactly 2 conditions."
            sys.exit()
        else:
            for sample in ss.Sample:
                if sample.replace(".bam","") not in bams:
                    return "I didn't find {} in your bam directory. Please correct.".format(sample)
                    sys.exit()
            diffDic = {}
            diffDic["Cond"] = list(ss.Cond.unique())
            diffDic["Comp"] = {}
            diffDic["Samples"] = list(ss.Sample.str.replace(".bam",""))
            for Comp in ss.Comp.unique():
                tempss = ss[ss['Comp'] == Comp]
                diffDic["Comp"][Comp] = {}
                key1 = tempss.Cond.unique()[0]
                key2 = tempss.Cond.unique()[1]
                diffDic["Comp"][Comp]['Cond'] = {}
                diffDic["Comp"][Comp]['Cond'][key1] = list(tempss[tempss['Cond'] == key1]['Sample'].str.replace(".bam",""))
                diffDic["Comp"][Comp]['Cond'][key2] = list(tempss[tempss['Cond'] == key2]['Sample'].str.replace(".bam",""))
                samplesList = diffDic["Comp"][Comp]['Cond'][key1] + diffDic["Comp"][Comp]['Cond'][key2]
                diffDic["Comp"][Comp]['Samples'] = samplesList
            for Comp in ss.Comp.unique():
                table = rich.table.Table(title="Samples - Conditions: {}".format(Comp))
                table.add_column(key1, justify="center", style="cyan")
                table.add_column(key2, justify="center", style="green")
                for i in range(max([len(diffDic["Comp"][Comp]['Cond'][key1]), len(diffDic["Comp"][Comp]['Cond'][key2])])):
                    try:
                        table.add_row(diffDic["Comp"][Comp]['Cond'][key1][i], diffDic["Comp"][Comp]['Cond'][key2][i])
                    except:
                        try:
                            table.add_row(diffDic["Comp"][Comp]['Cond'][key1][i], "")
                        except:
                            table.add_row("", diffDic["Comp"][Comp]['Cond'][key2][i])
                console = rich.console.Console()
                console.print(table)
            diffDic['baseDir'] = os.path.dirname(__file__)
            return diffDic
    else:
        return "Column headers not ok, (expected [Sample, Cond, Comp])"
        sys.exit()

def createTexfromTemplate(texfile, paramDic):
    if os.path.exists("Report.pdf"):
        print("Report exists already, cleaning and recompiling.")
        os.remove("Report.pdf")
    texOut = []
    with open(texfile) as f:
        for line in f:
            texOut.append(line.strip())
    #now start appending lines based on number of comparisons.
    for i in paramDic['Comp']:
        texOut.append("\n")
        texOut.append("\section{" + "{}".format(i) + "}")
        texOut.append("This section contains the comparison {}.".format(i))
        texOut.append("\n")
        texOut.append("\\begin{figure}[!htb]")
        texOut.append("\\begin{center}")
        texOut.append("\includegraphics[width=1\\textwidth]{Figures/" + "{}".format(i) + "_PCA.png}")
        texOut.append("\caption{" + "{}".format(i) + " PCA.}")
        texOut.append("\end{center}")
        texOut.append("\end{figure}")
        texOut.append("\n")
        texOut.append("\\begin{figure}[!htb]")
        texOut.append("\\begin{center}")
        texOut.append("\includegraphics[width=1\\textwidth]{Figures/" + "{}".format(i) + "_plotCorr.png}")
        texOut.append("\caption{" + "{}".format(i) + " correlation plot.}")
        texOut.append("\end{center}")
        texOut.append("\end{figure}")
        texOut.append("\n")
        texOut.append("\\begin{figure}[!htb]")
        texOut.append("\\begin{center}")
        texOut.append("\includegraphics[width=1\\textwidth]{Figures/" + "{}".format(i) + "_Heatmap.png}")
        texOut.append("\caption{" + "{}".format(i) + " heatmap.}")
        texOut.append("\end{center}")
        texOut.append("\end{figure}")
        texOut.append("\pagebreak")
    texOut.append("\end{document}")
    with open("Report.tex", "w") as f:
        for texLine in texOut:
            f.write("%s\n" % texLine)
    subprocess.call(['tectonic', 'Report.tex'])
    if os.path.exists("Report.pdf"):
        print("Report compilation done. Removing temporary file.")
        os.remove("Report.tex")
    print("Done")    