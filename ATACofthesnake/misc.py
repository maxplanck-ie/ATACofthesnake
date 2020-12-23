import os
import pandas as pd
import sys
import rich

def readBamDir(bamDir):
    bams = []
    for i in os.listdir(bamDir):
        if i.endswith('bam'):
            bams.append(i)
    return bams

def readSS(ss):
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
            diffDic[key1] = list(ss[ss['Cond'] == key1]['Sample'])
            diffDic[key2] = list(ss[ss['Cond'] == key2]['Sample'])
            samplesList = diffDic[key1] + diffDic[key2]
            table = rich.table.Table(title="Samples - Conditions.")
            table.add_column(key1, justify="center", style="cyan")
            table.add_column(key2, justify="center", style="green")
            for i in range(max([len(diffDic[key1]), len(diffDic[key2])])):
                try:
                    table.add_row(diffDic[key1][i], diffDic[key2][i])
                except:
                    try:
                        table.add_row(diffDic[key1][i], "")
                    except:
                        table.add_row("", diffDic[key2][i])
            console = rich.console.Console()
            console.print(table)
            return diffDic, samplesList
