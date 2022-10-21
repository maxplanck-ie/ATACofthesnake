# ATACofthesnake

Downstream processing of ATAC data, including QC's and differential accessibility. Starting point are deduplicated bam files, could be obtained from snakePipes (https://github.com/maxplanck-ie/snakepipes).

## Important

All samples in a 'run' have to belong to the same group, e.g. out of all individual peaks, a union will be made.
If you have multiple groups, you need multiple runs of the program to analyse them appropriately.

## Installation

  set up the environment:  
>  git clone git@github.com:maxplanck-ie/ATACofthesnake.git  
>  cd ATACofthesnake  
>  conda create -n ATACofthesnake python=3  
>  conda activate ATACofthesnake  
>  pip install ./  

## Quickstart

 > ATAC -h  

standard analysis:

 > ATAC --bamdir folder/with/bamfiles/ --outputdir outputfolder \
   --gtf genes.gtf --genomefasta genome.fa --genomesize 1.87e9 \
   --snakemakeprofile profile -b read_attracting_regions.bed

This will generate:
 - sieved bamfiles.
(cfr. --fragsize & --read_attracting_regions). Note that -b is obliged. At minimum this should contain the mitochondrial genome. Note that the mitochondrial contig is assumed to be named 'MT'. You can change this using --mitostring.
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

In this case, 2 analysis would be performed (under folders 'comparison1' & 'wtdmso_vs_wtdrug'). Note that the full design will be used anyhow (which is good for incorporating e.g. batch factors).



## Running:  

> ATAC -h  

> ATAC --bamDir [bamDir] --outDir [outputDir] --blackList [blacklist.bed] --genes [genes.gtf] --genomeFasta [genome.fa] --motifs [motifs.meme] --genomeSize [genomesize] --condaPrefix [path/to/conda]

  - Flags (required):
    - bamDir: Directory containing the bam files.
    - outDir: Directory for output
    - blackList: bed file containing regions to blacklist. Should at least contain the mitochondrial genome.
    - genes: gtf file containing gene annotations.
    - genomeFasta: fasta file containing genome of interest (uncompressed).
    - motifs: transcription factor motifs (meme format).
    - genomeSize: float or integer specifying effective genome size.
    - condaPrefix: path to your conda installation.

  - Flags (optional):
    - sampleSheet: To invoke differential accessibility. Provide tsv file containing Sample, Cond, Comp columns (header is required.). An additional columns 'Batch' can be used to include batch effects.
    - fragSizeMax: integer to set a maximum allowed fragment size. Defaults to 150.
    - peakSet: Provide a bed file containing peaks. If set, pipeline will use those instead of calling peaks with MACS2.
    - mergeBam: Flag to force merging of bamFiles (all of them if no sampleSheet is given, per condition if a sampleSheet is given) prior to peak calling.
    - clusterCMD: grid submission. defaults to SlurmEasy.
    
  sampleSheet example:

  | Sample | Cond | Comp |
  | -- | -- | -- |
  | WT1 | WT | WTvsKOa |
  | WT2 | WT | WTvsKOa |
  | KOa1 | KO | WTvsKOa |
  | KOa2 | KO | WTvsKOa |
  | WT1 | WT | WTvsKOb |
  | WT2 | WT | WTvsKOb |
  | KOb1 | KO | WTvsKOb |
  | KOb2 | KO | WTvsKOb |

  or if there is a batch effect:
  
  | Sample | Cond | Comp | Batch |
  | -- | -- | -- | -- |
  | WT1 | WT | WTvsKOa | Batch1 |
  | WT2 | WT | WTvsKOa | Batch2 |
  | KOa1 | KO | WTvsKOa | Batch1 |
  | KOa2 | KO | WTvsKOa | Batch2 |
