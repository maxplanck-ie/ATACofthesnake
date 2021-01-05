import os
import pandas as pd
import sys
import rich

from rich.progress import Progress
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

def setdefault_readss(ss):
    ss = pd.read_csv(ss, sep='\t', header=0)
    if ss.columns[0] == 'Sample' and ss.columns[1] == 'Cond':
        if len(ss.columns) > 2:
            return "Error parsing sampleSheet, more than 2 columns found."
            sys.exit()
        elif len(ss.Cond.unique()) != 2:
            return "Error parsing sampleSheet, I need exactly 2 conditions."
            sys.exit()
        else:
            diffDic = {}
            key1 = ss.Cond.unique()[0]
            key2 = ss.Cond.unique()[1]
            diffDic['Cond'] = {}
            diffDic['Cond'][key1] = list(ss[ss['Cond'] == key1]['Sample'].str.replace(".bam",""))
            diffDic['Cond'][key2] = list(ss[ss['Cond'] == key2]['Sample'].str.replace(".bam",""))
            samplesList = diffDic['Cond'][key1] + diffDic['Cond'][key2]
            table = rich.table.Table(title="Samples - Conditions.")
            table.add_column(key1, justify="center", style="cyan")
            table.add_column(key2, justify="center", style="green")
            for i in range(max([len(diffDic['Cond'][key1]), len(diffDic['Cond'][key2])])):
                try:
                    table.add_row(diffDic['Cond'][key1][i], diffDic['Cond'][key2][i])
                except:
                    try:
                        table.add_row(diffDic['Cond'][key1][i], "")
                    except:
                        table.add_row("", diffDic['Cond'][key2][i])
            console = rich.console.Console()
            console.print(table)
            diffDic['Samples'] = samplesList
            diffDic['baseDir'] = os.path.dirname(__file__)
            return diffDic