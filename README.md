# ATACofthesnake

Downstream processing of ATAC data, including QC's and differential accessibility. Starting point are deduplicated bam files, could be obtained from snakePipes (https://github.com/maxplanck-ie/snakepipes).


## Installation

  set up the environment:  
>  git clone git@github.com:maxplanck-ie/ATACofthesnake.git  
>  cd ATACofthesnake  
>  conda create -n ATACofthesnake python=3  
>  conda activate ATACofthesnake  
>  pip install ./  

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
