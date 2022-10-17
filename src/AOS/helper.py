import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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


def qcplotter(qcdir, _of):
    plotdic = {}
    for sample in glob.glob(
        os.path.join(qcdir, '*ix.tsv')
    ):
        plotdic[
            os.path.basename(sample).replace('_ix.tsv', '')  
        ] = idx_to_mit(sample)
    pldf = pd.DataFrame(plotdic.items())                                                                                                                                                                     
    pldf.columns = ['Sample', 'mitofraction']                                                                                                                                                                    
    pldf.sort_values(by='mitofraction', inplace=True)                                                                                                                                                            
    g = sns.barplot(data=pldf, x='Sample', y='mitofraction')                                                                                                                                                     
    g.set_xticklabels(g.get_xticklabels(), rotation=90)                                                                                                                                                          
    g.get_figure().savefig(_of, dpi=300, bbox_inches='tight') 
