import glob
import os
import shutil
from rich import print
import yaml

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
        design
    ):
        print("[red][bold]---------- preflight ----------[/red][/bold]")
        self.dirs = {
            'bamdir': os.path.abspath(bamdir),
            'outputdir': os.path.abspath(outputdir)
        }
        self.files = {
            'readattractingregions': readattractingregions,
            'gtf': gtf,
            'fna': genomefasta,
            'motif': motifs,
            'samplesheet': samplesheet,
        }
        self.vars = {
            'genomesize': genomesize,
            'fragsize': fragsize,
            'snakemakeprofile': snakemakeprofile,
            'samplesheet': samplesheet,
            'design': design
        }
        self.samples = glob.glob(
            os.path.join(os.path.abspath(bamdir), '*.bam')
        )
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

    def retconf(self):
        return(
            {
                'dirs': self.dirs,
                'files': self.files,
                'variables': self.vars,
                'samples': self.samples
            }
        )
