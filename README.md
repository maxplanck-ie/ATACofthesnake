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

> ATAC --bamDir ./bam/ --sampleSheet ss.tsv --blackList blacklist.bed --genes genes.gtf --genomeSize 1.87e9 --genomeFasta genome.fa --condaPrefix /path/to/conda/ [--fragSizeMax] [--downStream] --motifs /path/to/motif.meme

Output is stored in AOS folder created in the invocation directory.

  - Flags:
    - bamDir: Directory containing the bam files.  
    - sampleSheet: tsv file containing Sample, Cond, Comp columns (header is required.). An additional columns 'Batch' can be used to include batch in the Diff. acc. analysis.
    - blackList: bed file containing regions to remove from bam file (prior to peak calling). This is required. If you don't want to mask any sequences, it should contain a mitochondrial genome alone.
    - genes: gtf file containing genome annotations (TSS enrichments are based on 'transcript' features.) 
    - genomeSize: Effective genome size for organism of interest.  
    - genomeFasta: fasta file for organism of interest.  
    - downStream: Try downstream plotting of differential sites, and AME enrichment of up and down-regulated peaks respectively.
    - condaPrefix: point to your conda installation.
    - dryRun: dry run snakemake.
    - motifs: path to (meme) motif files to scan for in AME. required.

  - Optional:  
    - clusterCMD: submission command, default to SlurmEasy.  
    - snakeOpts: additional snakemake options, given as a csv (API style)  
    - downStream: Flag to try downstream analysis, e.g. motif calling.  
    - motifs: point to a motif database of interest (meme format). Required if you invoke --downStream.  
    - fragSizeMax: defaults to 150. Filters fragment lengths to this integer.  
    
  
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
## Todo:

 - [ ] call-summits mode or not -> motif or no.
 - [ ] (TSS enrichment cutoff)
 - [ ] TOBIAS
 - [ ] specify norm options (scalefactors): background or signal (background need sparsity tests.)
 - [ ] chromHMM - marks ?
 - [ ] chromVAR ?
 - [ ] incorporate outDir
 - [ ] Clean up output declinations.
 - [ ] ILP solver
 - [ ] Figures per comparison in a subfolder
 - [ ] lateX build PDF - tectonic vs something else.
 - [ ] Motif search vs background --> GC%test and #seq test --> warnings raised.
 - [0] allow choice for DESeq2 and edgeR or not (For now I ignore DESeq2 alltogether).
 - [ ] linter fix - snakemake
 - [ ] outdir variable rather than fixed AOS str
 - [ ] fragSize distribution plots.
 - [ ] testData
 - [ ] pytests
 - [ ] move plotter into a class.
 - [ ] keep eye on meme bioconda installation.

 
- Done:
 - [x] motif search: Implement MEME
 - [x] diffheatmap function test rather than hardflag.
 - [x] decryptify error message (e.g. stop retries).
 - [x] paramLogs put in yaml
 - [x] shift condition definition specific to comparison.
 - [x] PEP8
 - [x] clean conda envs
 - [x] in diffPlots, add Condition labels
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
