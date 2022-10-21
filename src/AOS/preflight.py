import glob
import os
import shutil
from rich import print
import yaml
import pandas as pd
import sys

class Preflight():
    '''
    Bundles the config, preflight checks & samplesheet parsing
    '''
    def __init__(
        self,
        bamdir,
        outputdir,
        gtf,
        genomefasta,
        genomesize,
        readattractingregions,
        motifs,
        fragsize,
        snakemakeprofile,
        samplesheet,
        comparison,
        interaction,
        mitostring
    ):
        def retabspath(_p):
            if _p:
                return(os.path.abspath(_p))
            else:
                return('')

        print("[red][bold]---------- preflight ----------[/red][/bold]")
        if comparison and not samplesheet:
            print("You need to provide a samplesheet when providing asking for a comparison.")
            sys.exit()
        self.dirs = {
            'bamdir': os.path.abspath(bamdir),
            'outputdir': os.path.abspath(outputdir),
            'scriptsdir': os.path.dirname(__file__)
        }
        self.files = {
            'readattractingregions': retabspath(readattractingregions),
            'gtf': retabspath(gtf),
            'fna': retabspath(genomefasta),
            'motif': retabspath(motifs),
            'samplesheet': retabspath(samplesheet),
            'comparison': retabspath(comparison)
        }
        self.vars = {
            'genomesize': genomesize,
            'fragsize': fragsize,
            'snakemakeprofile': snakemakeprofile,
            'samplesheet': samplesheet,
            'mitostring': mitostring
        }
        self.samples = [os.path.basename(x).replace('.bam', '') for x in glob.glob(
            os.path.join(os.path.abspath(bamdir), '*.bam')
        )]
        self.envs = {
            'seqtools': os.path.join(
                self.dirs['scriptsdir'], 'envs', 'seqtools.yml'
            ),
            'tobias': os.path.join(
                self.dirs['scriptsdir'], 'envs', 'tobias.yml'
            )
        }
        self.rules = {
            'wf': os.path.join(
                self.dirs['scriptsdir'], 'workflow.smk'
            ),
            'peaks': os.path.join(
                self.dirs['scriptsdir'], 'rules', 'peaks.smk'
            ),
            'qc': os.path.join(
                self.dirs['scriptsdir'], 'rules', 'qc.smk'
            ),
            'de': os.path.join(
                self.dirs['scriptsdir'], 'rules', 'DE.smk'
            )
        }
        self.rscripts = {
            'scalefactors': os.path.join(
                self.dirs['scriptsdir'], 'rscripts', 'scalefactors.R'
            ),
            'edger': os.path.join(
                self.dirs['scriptsdir'], 'rscripts', 'edger.R'
            )
        }
        if interaction:
            self.interaction = '*'
        else:
            self.interaction = '+'

    def dumpconf(self):
        '''
        makes output directory if not existing, dumps config yaml in there.
        '''
        _odir = self.dirs['outputdir']
        _conf = os.path.join(
            _odir,
            'config.yaml'
        )
        if not os.path.exists(_odir):
            print("Creating {}".format(_odir))
            os.mkdir(_odir)
        with open(_conf, 'w') as f:
            yaml.dump(
                self.retconf(),
                f,
                default_flow_style=False
            )
        self.files['configfile'] = _conf

    def retconf(self):
        '''
        Returns itself as a dictionary.
        '''
        return(
            {
                'dirs': self.dirs,
                'files': self.files,
                'vars': self.vars,
                'samples': self.samples,
                'rules': self.rules,
                'envs': self.envs,
                'rscripts': self.rscripts,
                'comparison': self.comparison,
                'factors': self.factors,
                'interaction': self.interaction
            }
        )
    
    def checkcomps(self):
        if self.files['samplesheet']:
            if not self.files['comparison']:
                sys.exit("Providing a samplesheet requires a comparisons files as well.")
            else:
                ssdf = pd.read_csv(self.files['samplesheet'], sep='\t', header=0)
                if list(ssdf.columns)[0].lower() != 'sample':
                    sys.exit("First column in samplesheet needs to be sample. Exiting.")
                factors = list(ssdf.columns)[1::]
                with open(self.files['comparison']) as f:
                    cdat = yaml.safe_load(f)
                # Double Check.
                for comp in cdat:
                    for group in cdat[comp]:
                        for factor in cdat[comp][group]:
                            if factor not in ssdf.columns:
                                sys.exit("{} not in samplesheet. exiting.".format(factor))
                self.comparison = cdat
                self.factors = factors
        else:
            self.comparison = ''
            self.factors = ''

    def genTSS(self):
        _sortgtfo = os.path.join(
            self.dirs['outputdir'],
            'genes.sorted.gtf'
        )
        _tsso = os.path.join(
            self.dirs['outputdir'],
            'tss.bed'
        )
        if not os.path.exists(self.dirs['outputdir']):
            os.mkdir(
                self.dirs['outputdir']
            )
        if not os.path.exists(_sortgtfo) or not os.path.exists(_tsso):
            GTF = pd.read_csv(
                self.files['gtf'],
                sep='\t',
                comment='#',
                header=None,
                dtype={0: 'str'}
            )
            GTF.columns = [
                'chr',
                'source',
                'feature',
                'start',
                'end',
                'score',
                'strand',
                'frame',
                'attribute'
            ]
            GTF = GTF[GTF['feature'] == 'transcript']
            GTF = GTF.sort_values(
                ["chr", "start"],
                ascending=(True, True)
            )
            GTF.to_csv(
                _sortgtfo,
                header=False,
                index=False,
                sep='\t'
            )
            TSS = []
            for ix, row in GTF.iterrows():
                l = list(row)
                if l[6] == '+':
                    TSS.append(
                        [
                            l[0],
                            int(l[3]) - 1,
                            int(l[3])
                        ]
                    )
                elif l[6] == '-':
                    TSS.append(
                        [
                            int(l[4]) - 1,
                            int(l[4]),
                        ]
                    )
            TSSdf = pd.DataFrame(TSS).drop_duplicates()
            TSSdf.to_csv(
                _tsso,
                sep='\t',
                index=False,
                header=False
            )
            self.files['gtf'] = _sortgtfo
            self.files['tss'] = _tsso
        self.files['gtf'] = _sortgtfo