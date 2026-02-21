from pathlib import Path
import yaml
import pandas as pd
from aos.logger import setup_logger
from importlib.metadata import version
import re

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
            self.validate_samplesheet(self.samplesheet, self.bamdir)
            self.validate_comparison(self.comparison, self.samplesheet)
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
    def _extract_design_vars(formula: str) -> list[str]:
        rhs = formula.replace('~', '')
        parts = re.split(r'[+*:]', rhs)
        vars_ = []
        for part in parts:
            var = part.strip()
            if var and var != '1':
                vars_.append(var)
        return vars_

    @staticmethod
    def validate_comparison(comp_path: Path, ss: Path) -> None:
        samplesheet = pd.read_csv(ss, sep='\t', header=0)
        covariates = samplesheet.columns.drop('sample').tolist()
        with open(comp_path, 'r') as f:
            comps = yaml.safe_load(f)
        for _ in comps:
            comp = comps[_]
            # Type is present
            assert 'type' in comp, f"Comparison {comp} does not have `type` key. Exiting."
            assert comp['type'] in ['twogroup', 'lrt', 'timecourse'], f"Comparison {comp} has invalid type {comp['type']}. Type should be one of 'twogroup', 'lrt' or 'timecourse'. Exiting."
            # Design is valid if present
            if 'design' in comp:
                assert '~' in comp['design'], f"Comparison {comp} has a design that does not contain '~'. Exiting."
                _designvars = Preflight._extract_design_vars(comp['design'])
                for var in _designvars:
                    assert var in covariates, f"Comparison {comp} has design variable {var} that is not present in samplesheet. Exiting."
                # Assert principle of marginality is not violated
            if comp['type'] == 'twogroup':
                _ = comp.copy()
                _.pop('design', None)
                _.pop('type', None)
                assert len(_.keys()) == 2, f"Comparison {comp} of type 'twogroup' should only have two groups specified in the comparison file. Exiting."
                for group in _:
                    for factor in _[group]:
                        assert factor in covariates, f"Comparison {comp} has factor {factor} that is not present in samplesheet. Exiting."
                        if isinstance(_[group][factor], str):
                            assert _[group][factor] in samplesheet[factor].unique(), f"Comparison {comp} has level {_[group][factor]} for factor {factor} that is not present in samplesheet. Exiting."
                        else:
                            for level in _[group][factor]:
                                assert level in samplesheet[factor].unique(), f"Comparison {comp} has level {level} for factor {factor} that is not present in samplesheet. Exiting."
            if comp['type'] == 'lrt':
                assert 'reduced' in comp, f"Comparison {comp} of type 'lrt' should have a 'reduced' key specifying the reduced model. Exiting."
                assert '~' in comp['reduced'], f"Comparison {comp} has a reduced design that does not contain '~'. Exiting."
                _designvars = Preflight._extract_design_vars(comp['reduced'])
                for var in _designvars:
                    assert var in covariates, f"Comparison {comp} has reduced design variable {var} that is not present in samplesheet. Exiting."
                # Check that reduced is in fact, reduced.
            if comp['type'] == 'timecourse':
                assert 'time' in comp, f"Comparison {comp} of type 'timecourse' should have a 'time' key specifying the time variable. Exiting."
                assert 'time_type' in comp, f"Comparison {comp} of type 'timecourse' should have a 'time_type' key specifying the time variable type (continuous or ordinal). Exiting."
                assert comp['time_type'] in ['continuous', 'ordinal'], f"Comparison {comp} has invalid time_type {comp['time_type']}. Should be 'continuous' or 'ordinal'. Exiting."
                assert comp['time'] in covariates, f"Comparison {comp} has time variable {comp['time']} that is not present in samplesheet. Exiting."
                if comp['time_type'] == 'ordinal':
                    assert 'order' in comp, f"Comparison {comp} of type 'timecourse' with time_type 'ordinal' should have an 'order' key specifying the order of the time points. Exiting."
                    for level in comp['order']:
                        assert level in samplesheet[comp['time']].unique(), f"Comparison {comp} has level {level} for time variable {comp['time']} that is not present in samplesheet. Exiting."
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