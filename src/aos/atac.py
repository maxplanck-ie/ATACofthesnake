import subprocess

import click
from rich.console import Console

from aos.preflight import Preflight


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-i",
    "--bamdir",
    required=True,
    type=click.Path(exists=True),
    help="Specify directory that contains your bamfiles.",
)
@click.option("-o", "--outputdir", required=True, help="Specify output directory.")
@click.option(
    "-g",
    "--gtf",
    type=click.Path(exists=True),
    required=True,
    help="Specify a gtf file containing gene annotation. Will be used to extract TSS.",
)
@click.option(
    "-r",
    "--genomefasta",
    type=click.Path(exists=True),
    required=True,
    help="Specify a fasta file that contains the reference genome.",
)
@click.option(
    "-b",
    "--readattractingregions",
    type=click.Path(exists=True),
    required=True,
    help="Specify a bed file containing read attracting regions. Should contain the mitochondrial genome at least.",
)
@click.option(
    "-p",
    "--snakemakeprofile",
    required=False,
    help="specify the name of your snakemake profile.",
)
@click.option(
    "-@",
    "--threads",
    default=1,
    type=int,
    show_default=True,
    help="specify the number of threads to use. Only relevant if no snakemake profile is given.",
)
@click.option(
    "-m",
    "--motifs",
    type=click.Path(exists=True),
    help="Specify a file containing motifs. Needs to be in meme format. If not provided, no motif analyses will be ran.",
)
@click.option(
    "-f",
    "--fragsize",
    default=150,
    type=int,
    show_default=True,
    help="Specify the maximum fragment size (bps) to be considered for peak calling. Sits at 150 bps by default to only use reads from nucleosome-free regions.",
)
@click.option(
    "--samplesheet",
    type=click.Path(exists=True),
    help="specify a samplesheet (as a tsv file). Refer to the documentation on how to format this file.",
)
@click.option(
    "--comparison",
    type=click.Path(exists=True),
    help="specify yaml file with comparisons. Required if a samplesheet is given. Refer to the documentation on how to format this file.",
)
@click.option(
    "--mitostring",
    required=False,
    default="MT",
    show_default=True,
    help="Name of the mitochondrial contig (as in the reference genome / BAM file). Defaults to MT.",
)
@click.option(
    "--upstreamuro",
    required=False,
    default=50000,
    show_default=True,
    help="Maximum permitted distance upstream of a feature (peak annotation).",
)
@click.option(
    "--downstreamuro",
    required=False,
    default=50000,
    show_default=True,
    help="Maximum permitted distance downstream of a feature (peak annotation).",
)
@click.option(
    "--featureuro",
    required=False,
    default="gene",
    show_default=True,
    help="the feature in the GTF file (column 3) to use for peak annotation.",
)
@click.option(
    "--pseudocount",
    required=False,
    default=8,
    show_default=True,
    help="Pseudocount to add to the count matrix prior to differential calling.",
)
@click.option(
    "--peakset",
    required=False,
    type=click.Path(exists=True),
    help="Include an external peak file (bed format). If not provided, the union of all peaks across all samples will be used to generate a count matrix.",
)
@click.option(
    "--permutation_cutoff",
    required=False,
    default=1e-2,
    show_default=True,
    help="A p-value cutoff to determine significance after permutation testing (relevant for timecourse mode).",
)
@click.option(
    "--permutation_iterations",
    required=False,
    default=1000,
    show_default=True,
    help="Number of permutations to perform for significance testing (relevant for timecourse mode).",
)
@click.option(
    "--fdr_cutoff",
    required=False,
    default=1e-3,
    show_default=True,
    help="The FDR cutoff to determine significance after wald/LRT test (relevant for two-group and LRT modes).",
)
@click.option(
    "--lfc_cutoff",
    required=False,
    default=0,
    show_default=True,
    help="The log2FoldChange cutoff used to determine significance (combined with fdr_cutoff, only relevant for two-group mode).",
)
@click.option(
    "--lrt_peaks",
    required=False,
    default=1000,
    show_default=True,
    help="Minimum number of significant peaks in LRT mode to continue with downstream analysis.",
)
@click.option(
    "--gp_timesteps",
    required=False,
    default=10,
    show_default=True,
    help="Number of timesteps to use for gaussian process regression.",
)
def main(**kwargs) -> None:
    pf = Preflight(kwargs)
    console = Console()
    with console.status("[bold green] Running snakemake..."):
        subprocess.run(
            [
                "snakemake",
                "-s",
                pf.workflowfile,
                "--configfile",
                pf.configfile,
                "-d",
                pf.outputdir,
                "-p",
            ]
            + pf.snakemake_arguments()
        )
