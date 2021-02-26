# ATACofthesnake

Downstream processing of ATAC data, including QC's and differential accessibility. Starting point are deduplicated bam files, could be obtained from snakePipes (https://github.com/maxplanck-ie/snakepipes).


## Installation

  set up the environment:  
>  git clone git@github.com:WardDeb/ATACofthesnake.git  
>  cd ATACofthesnake  
>  conda create -n ATACofthesnake python=3  
>  conda activate ATACofthesnake  
>  pip install ./  

## Running:  

> ATAC --bamDir ./bam/ --outDir ./ --sampleSheet ss. --blackList blacklist.bed --Genes genes.gtf --genomeSize 1.87e9 --genomeFasta genome.fa --diffPeak

Currently creates output in working directory!

  - Flags:
    - bamDir: Directory containing the bam files.  
    - (outDir: Ouput directory (future).)  
    - sampleSheet: tsv file containing Sample, Cond, Comp columns (header is required.).  
    - blackList: bed file containing regions to remove from bam file (prior to peak calling).  
    - Genes: gtf file containing genome annotations.  
    - genomeSize: Effective genome size for organism of interest.  
    - genomeFasta: fasta file for genome.  
    - diffPeak: Wether or not to try and plot a heatmap of diff. Acc regions.
  
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


## Todo:

 - [ ] call-summits mode or not -> motif or no.
 - [ ] (TSS enrichment cutoff)
 - [x] in diffPlots, add Condition labels
 - [ ] TOBIAS
 - [ ] motif search: Implement MEME
 - [ ] specify norm options (scalefactors): background or signal (background need sparsity tests.)
 - [ ] diffheatmap function test rather than hardflag.
 - [ ] chromHMM - marks ?
 - [ ] chromVAR ?
 - [ ] clean rule all declinations / incorporate outDir.
 - [ ] ILP solver
 - [ ] Figures per comparison in a subfolder
 - [x] decryptify error message (e.g. stop retries).
 - [x] paramLogs put in yaml
 - [ ] lateX build PDF
 - [ ] Motif search vs background --> GC%test and #seq test --> warnings raised.
 - [ ] allow choice for DESeq2 and edgeR or not.
 - [ ] shift condition definition specific to comparison.
 - [ ] linter
 - [ ] clean conda envs

 
- Done:
 - [x] FrIPs
 - [x] MAplot
 - [x] index: generate bai files if they are missing.
 - [x] mitoBleed
 - [x] deeptools: fragmentsize
 - [x] for all arguments, make misc function that checks existence, otherwise fault before running.
 - [x] MACS2 - summits
 - [x] deal with illegal characters ?
 - [x] slurm submission
 - [x] produce report
 - [x] incorporate runID if more than 1 comparison per batch.
 - [x] Multiple comparisons
 - [x] annotate peaks
 - [x] Deal with multiple envs.
 - [x] Standardize logs
 - [x] scaleFactors
 - [x] bamCoverage
 - [x] computeMatrix
 - [x] deeptools: correlations
 - [x] shift all QC in a QC folder.
