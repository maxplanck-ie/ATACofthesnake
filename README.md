# ATACofthesnake

Downstream processing of ATAC data, including QC's and differential accessibility. Starting point are deduplicated bam files, could be obtained from snakePipes (https://github.com/maxplanck-ie/snakepipes).

## Important

All samples in a 'run' have to belong to the same group. That means that out of all peaks/sample, a union will be made.

Fasta headers un field 0 (space delimited) are not allowed to contain a pipe character '|'.

## Installation

  From github:
>  git clone git@github.com:maxplanck-ie/ATACofthesnake.git  
>  pixi run ATAC -h

  From pypi:
>  pip install ATACofthesnake

## Quickstart

standard analysis:

 > ATAC --bamdir folder/with/bamfiles/ --outputdir outputfolder \
   --gtf genes.gtf --genomefasta genome.fa \
   --snakemakeprofile profile -b read_attracting_regions.bed

Note that this pipeline depends on snakemake. Additionally, the snakemake environments are managed using conda, which means you need to have conda installed and configured. Make sure this is set in your snakemake profile (if you use one), when not using a profile, snakemake will be ran with the --use-conda flag by default.

The default analysis will generate:

 - sieved bamfiles.
(cfr. --fragsize/-f & --read_attracting_regions/-b). Note that read attracting regions are obliged. At minimum this should contain the mitochondrial genome. Note that the mitochondrial contig is assumed to be named 'MT'. You can change this using --mitostring.
 - peaks called per bamfile (under peaks/) and a union of all peaks (peakset/).
 - bigwigs normalized using scalefactors and RPKM.
 - a number of QC plots (figures/) and metrics (qc/)

## Differential analysis

Differential analysis requires the additional yaml file specifying the comparison (comparison.yaml) and the samplesheet.

 > ATAC --bamdir folder/with/bamfiles/ --outputdir outputfolder \
   --gtf genes.gtf --genomefasta genome.fa --genomesize 1.87e9 \
   --snakemakeprofile profile -b read_attracting_regions.bed \
   --comparison comparison.yaml --samplesheet samplesheet.tsv

The samplesheet has to be a tab-separated file with the first column containing the sample names (as the bam files are named), and the other columns containing the factors that will make up the design. For example:

| sample | genotype | treatment |
| -- | -- | -- |
| sample1 | WT | DMSO |
| sample2 | WT | DMSO |
| sample3 | WT | drug |
| sample4 | WT | drug |
| sample5 | KO | DMSO |
| sample6 | KO | DMSO |
| sample7 | KO | drug |
| sample8 | KO | drug |

will result in the following design:

 > ~ genoype + treatment

or if --interaction is set:

 > ~ genotype * treatment  

Note the comparisons made need to be specified in the comparison.yaml file, a samplesheet alone is not enough. An example together with the above samplesheet would be:

```bash
comparison1:
  group1:
    genotype: 'WT'
  group2:
    genotype: 'KO'
wtdmso_vs_wtdrug:
  wtdmso:
    genotype: 'WT'
    treatment: 'DMSO'
  wtdrug:
    genotype: 'WT'
    treatment: 'drug'
```

In this case, 2 analyses would be performed (under folders 'comparison1' & 'wtdmso_vs_wtdrug'). In both cases all of the factor columns in the samplesheet will be used in the design.

A motif file in [meme](https://meme-suite.org/meme/doc/meme-format.html) format can be supplied as well. If this is the case, the motifs will first be clustered, and motif enrichments will be calculated for all differential regions (using [ame](https://meme-suite.org/meme/doc/ame.html)). Note that the 'reciprocal' differential peaks will be used as background. E.g. group2 differential peaks are background for group1 enrichments and vice versa. Note that this approach could be biased if there is a strong disbalance between those groups.
