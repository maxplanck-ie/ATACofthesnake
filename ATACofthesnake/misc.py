import os
import pandas as pd
import sys
import rich
from rich.progress import Progress

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


def GTFtoTSS(bed):
    TSS = []
    linecount = 0
    with open(bed) as f:
        for line in f:
            linecount += 1
    with rich.progress.Progress() as progress:
        task = progress.add_task("Extracting TSS from GTF", total=linecount)
        with open(bed) as f:
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