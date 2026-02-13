from pathlib import Path
import yaml
import pandas as pd
import jsonschema
from aos.logger import setup_logger
from importlib.metadata import version

class Preflight():
    '''
    Bundles the config, preflight checks & samplesheet parsing
    '''
    def __init__(self, clickdict: dict) -> None:
        # Store paths - required
        self.bamdir = Path(clickdict['bamdir'])
        self.outputdir = Path(clickdict['outputdir'])

        # Create outputdir if not existing
        self.outputdir.mkdir(exist_ok=True, parents=True)
        # Set up log
        self.logfile = self.outputdir / 'aos.log'
        self.logger = setup_logger(self.logfile)
        self.logger.info(f"AOS - version {version('ATACofthesnake')}")
        self.logger.info("--" * 10)
        
        self.logger.info("Checking mandatory paths...")
        self.rar = Path(clickdict['readattractingregions'])
        self.gtf = Path(clickdict['gtf'])
        self.fna = Path(clickdict['genomefasta'])

        # Figure out workflow file path
        self.logger.info("Finding snakemake workflow file...")
        self.workflowfile = Path(__file__).parents[0] / 'workflow.smk'
        assert self.workflowfile.exists(), f"{self.workflowfile} does not exist next to {__file__}. Please check your installation."


        # Store paths - optional
        self.logger.info("Parsing optional paths...")
        self.motif = self.optional_paths(clickdict['motifs'])
        self.samplesheet = self.optional_paths(clickdict['samplesheet'])
        self.comparison = self.optional_paths(clickdict['comparison'])
        self.peakset = self.optional_paths(clickdict['peakset'])

        
        if self.samplesheet and self.comparison:
            self.logger.info("Validating samplesheet and comparison file...")
            # Validate samplesheet and comparison file.
            #self.validate_comparison(self.samplesheet, self.comparison)
            self.validate_samplesheet(self.samplesheet, self.bamdir)
        else:
            self.logger.info("No differential testing requested...")
            self.interaction = None

        # Store variables
        self.logger.info("Parsing other variables...")
        self.fragsize = clickdict['fragsize']
        self.mitostring = clickdict['mitostring']
        self.uro = [
            clickdict['featureuro'],
            clickdict['upstreamuro'],
            clickdict['downstreamuro'],
        ]
        self.pseudocount = clickdict['pseudocount']


        # Run settings
        self.logger.info("Parsing run settings...")
        self.smk_profile = clickdict['snakemakeprofile']
        self.threads = clickdict['threads']

        # Parse fasta file
        self.logger.info("Validating fasta file and calculating effective genome size...")
        self.parse_fasta()
        self.logger.info(f"Effective genome size calculated to be {self.ESS}.")

        # Parse gtf file
        self.logger.info("Parsing GTF file and creating TSS file...")
        self.parse_gtf()

        self.logger.info("Preflight complete. Configuration:")
        self.logger.debug(self.__dict__)

        # Validate mitostring
        self.logger.info("Making sure mitostring is present in fasta file...")
        self.validate_mitostring(self.fna, self.mitostring, self.rar)

        # Write out config file
        self.logger.info("Writing config file to output directory...")
        self.configfile = self.outputdir / 'aos_config.yaml'

        with open(self.configfile, 'w') as f:
            yaml.dump(self.to_dict(),f,default_flow_style=False,sort_keys=False)
    
    @staticmethod
    def optional_paths(path: str|None) -> Path|None:
        if path:
            return Path(path)
        else:
            return None
    
    # @staticmethod
    # def validate_comparison(ss: Path, comp: Path) -> None:
    #     ssdf = pd.read_csv(ss, sep='\t', header=0)
    #     assert any([c == 'sample' for c in ssdf.columns]), "No 'sample' column found in samplesheet. Exiting."
    #     factors = [ c for c in ssdf.columns if c != 'sample' ]
    #     with open(comp) as f:
    #         cdat = yaml.safe_load(f)
    #     # Double check that all factors in the comparison file are actually present.
    #     _factors = []
    #     for comp in cdat:
    #         for group in cdat[comp]:
    #             for factor in cdat[comp][group]:
    #                 _factors.append(factor)
    #     # Check that all factors in the columns
    #     assert all([factor in factors for factor in _factors]), f"Not all factors in comparison file {set(_factors)} are present in samplesheet {ssdf.columns}. Exiting."
    
    @staticmethod
    def validate_mitostring(fna, mitostring, rarfile):
        with open(fna, 'r') as f:
            headers = [ line.strip() for line in f if line.startswith('>') ]
        assert any([ mitostring in header for header in headers ]), f"Provided mitostring {mitostring} not found in fasta file. Please check your fasta file and the --mitostring parameter. Exiting."
        with open(rarfile, 'r') as f:
            rarheaders = [ line.strip() for line in f if not line.startswith('#') ]
        assert any([ mitostring in header for header in rarheaders ]), f"Provided mitostring {mitostring} not found in read attracting regions file. Please check your rar file and the --mitostring parameter. Exiting."
    
    @staticmethod
    def validate_samplesheet(ss: Path, bamdir: Path) -> None:
        df = pd.read_csv(ss, sep='\t', header=0)
        for sample in df['sample']:
            bamfile = bamdir / f"{sample}.bam"
            cramfile = bamdir / f"{sample}.cram"
            assert bamfile.exists() or cramfile.exists(), f"Bamfile or cramfile not found for sample {sample} in samplesheet. Exiting."
        # All unique factor combinations should have at least two replicates.
        counts = df.groupby(df.columns.drop("sample").tolist()).size()
        assert (counts >= 2).all(), \
            f"Some factor combinations have < 2 replicates:\n{counts[counts < 2]}"

    @staticmethod
    def validate_comparisonentry(compentry: dict, samplesheet: Path) -> None:
        # Validate type is present
        # Validate design -> ~, +, :, *, factors in samplesheet
        # if type == twogroup, only two keys if type and design are dropped.
        # if type == twogroup, validate factors and values in groups are present in samplesheet
        None

    def parse_fasta(self):
        ESS = 0
        with open(self.fna, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    assert '|' not in line, f" '|' are not allowed in fasta headers. Please reformat. Offending header: {line.strip()}"
                else:
                    ESS += len(line.strip()) - line.strip().lower().count('n')
        self.ESS = ESS

    def parse_gtf(self):
        gtf = pd.read_csv(self.gtf, sep='\t', comment='#', header=None, dtype={0: 'str'})
    
        gtf.columns = [
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
        gtf = gtf[gtf['feature'] == self.uro[0]]
        assert not gtf.empty, f"Feature {self.uro[0]} not found in GTF file. Please check your GTF file and the --featureuro parameter. Exiting."

        gtf = gtf.sort_values(["chr", "start"],ascending=(True, True))
        gtf.to_csv(self.outputdir / 'genes.sorted.gtf', header=False, index=False, sep='\t')

        tss = []
        for _, row in gtf.iterrows():
            gtfentry = list(row)
            if gtfentry[6] == '+':
                tss.append( [ gtfentry[0], int(gtfentry[3]) - 1, int(gtfentry[3]) ] )
            elif gtfentry[6] == '-':
                tss.append( [ gtfentry[0], int(gtfentry[4]), int(gtfentry[4]) + 1 ] )
        tss = pd.DataFrame(tss).drop_duplicates()
        tss.to_csv(self.outputdir / 'tss.bed', sep='\t', index=False, header=False)
        self.gtf_sorted = self.outputdir / 'genes.sorted.gtf'
        self.tss = self.outputdir / 'tss.bed'
    
    def snakemake_arguments(self) -> list:
        if self.smk_profile:
            return ['--profile', self.smk_profile]
        else:
            return ['--cores', f"{self.threads}", "--use-conda"]
    
    def to_dict(self):
        d = {}
        for k, v in self.__dict__.items():
            if k == 'logger':
                continue
            if isinstance(v, Path):
                d[k] = str(v.absolute())
            else:
                d[k] = v
        return d