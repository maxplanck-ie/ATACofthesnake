# ATACofthesnake

Downstream processing of ATAC data, including QC's and differential accessibility. Starting point are deduplicated bam files, could be obtained from snakePipes (https://github.com/maxplanck-ie/snakepipes).


  - Installation

  set up the environment:  
>  git clone git@github.com:WardDeb/ATACofthesnake.git  
>  cd ATACofthesnake  
>  conda create -n ATACofthesnake python=3  
>  conda activate ATACofthesnake  
>  pip install ./  

  - Running:  
> ATAC --bamDir ./bams/ --outDir ./ --sampleSheet ss.tsv --blackList blackList.bed --genomeSize 2652783500 --Genes genes.gtf

  - todo:

 - [x] index: generate bai files if they are missing.
 - [x] deeptools: fragmentsize
 - [x] mergeBams
 - [x] MACS2
 - [x] multibamsummary (counts)
 - [x] scaleFactors
 - [x] bamCoverage
 - [x] computeMatrix
 - [x] deeptools: TSS enrichment
 - [x] deeptools: correlations
 - [x] DESeq2
 - [x] Multiple comparisons
 - [x] annotate peaks
 - [x] Deal with multiple envs.
 - [x] Standardize logs
 - [ ] TOBIAS
 - [ ] motif search?
 - [ ] specify norm options (scalefactors): background or signal
 - [ ] diffheatmap
 - [ ] IDR values for replicates (https://github.com/kundajelab/idr)
 - [x] FrIPs
 - [ ] MAplot
 - [ ] chromHMM - marks ?

 - longer term to do:

 - [ ] EdgeR / csaw ?
 - [x] deal with illegal characters ?
 - [x] slurm submission
 - [x] produce report
 - [x] incorporate runID if more than 1 comparison per batch.
 - [ ] Support multiple genomes
      - [ ] blacklist DL + merge
      - [ ] GTF DL
      - [ ] genomeSizes

