import rich_click as click
from rich.console import Console
from rich import inspect
import os
import snakemake
import shutil
import yaml
from AOS.preflight import Preflight


@click.command(
        context_settings=dict(
            help_option_names=["-h", "--help"]
            )
    )
@click.option(
    '-i',
    '--bamdir',
    required=True,
    type=click.Path(exists=True),
    help='Specify directory that contains your bamfiles.'
)
@click.option(
    '-o',
    '--outputdir',
    required=True,
    help='Specify output directory.'
)
@click.option(
    '-g',
    '--gtf',
    type=click.Path(exists=True),
    required=True,
    help='Specify a gtf file containing gene annotation. Will be used to extract TSS.'
)
@click.option(
    '-r',
    '--genomefasta',
    type=click.Path(exists=True),
    required=True,
    help='Specify a fasta file that contains the reference genome.'
)
@click.option(
    '-s',
    '--genomesize',
    required=True,
    type=float,
    help='Specify the effective genome size, either as int (1870000000) or float (1.87e9)'
)
@click.option(
    '-p',
    '--snakemakeprofile',
    required=True,
    default='',
    help='specify the name of your snakemake profile.'
)
@click.option(
    '-b',
    '--readattractingregions',
    type=click.Path(exists=True),
    required=True,
    help='Specify a bed file containing read attracting regions. Should contain the mitochondrial genome at least.'
)
@click.option(
    '-m',
    '--motifs',
    type=click.Path(exists=True),
    help='Specify a file containing motifs. Needs to be in meme format.'
)
@click.option(
    '-f',
    '--fragsize',
    default=150,
    type=int,
    help='Specify the maximum fragment size (bps). Sits at 150 bps by default to capture only NFR (nucleosome-free region).'
)
@click.option(
    '--samplesheet',
    default='',
    help='specify a samplesheet (as a tsv file). See Readme for formatting.'
)
@click.option(
    '--comparison',
    default='',
    help='specify yaml file with comparisons. Required if a samplesheet is given.'
)
@click.option(
    '--interaction',
    default=False,
    is_flag=True,
    help='Wether or not to add interactions in the differential calculations. (e.g. ~factor1*factor2 is set, ~factor1+factor2 if not set).'
)
@click.option(
    '--mitostring',
    required=False,
    default='MT',
    help='Name of the mitochondrial contig. Defaults to MT.'
)
def main(bamdir,
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
    # Init
    pf = Preflight(**locals())
    # GTF
    print("Sorting GTF & creating TSS.bed..")
    pf.genTSS()
    # comparisons.
    pf.checkcomps()
    # Write conf
    print("Writing config file in {}..".format(
        os.path.basename(pf.dirs['outputdir'])
    ))
    #
    pf.dumpconf()
    #inspect(pf)
    console = Console()
    with console.status("[bold green] Running snakemake..."): 
        snakemake.main(
            [
                '-s', pf.rules['wf'],
                '--profile', pf.vars['snakemakeprofile'],
                '--max-jobs-per-second', '1',
                '-p',
                '--configfile', pf.files['configfile'],
                '-d', pf.dirs['outputdir']
            ]
        )