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
 - [ ] for all arguments, make misc function that checks existence, otherwise fault before running.
 - [ ] MACS2 - summits
 - [ ] extend summits
 - [ ] TSS enrichment cutoff
 - [ ] peak merging on TSS cut passing samples
 - [ ] peak normalization / depth ?
 - [ ] multibamsummary (counts)
 - [x] scaleFactors
 - [x] bamCoverage
 - [x] computeMatrix
 - [x] deeptools: correlations
 - [x] DESeq2
 - [x] Multiple comparisons
 - [x] annotate peaks
 - [x] Deal with multiple envs.
 - [x] Standardize logs
 - [ ] TOBIAS
 - [ ] motif search - homer / MEME
 - [ ] specify norm options (scalefactors): background or signal
 - [ ] diffheatmap
 - [x] FrIPs
 - [ ] MAplot
 - [ ] chromHMM - marks ?
 - [ ] chromVAR
 - [ ] shift all QC in a QC folder.
 - [x] mitoBleed
 - [ ] clean rule all declinations.
 - [ ] Figures per comparison in a subfolder
 - [ ] diffHeat DAG problem ?
 - [ ] specify all input per 'module'

 - longer term to do:

 - [x] deal with illegal characters ?
 - [x] slurm submission
 - [x] produce report
 - [x] incorporate runID if more than 1 comparison per batch.
 - [ ] Support multiple genomes
      - [ ] blacklist DL + merge
      - [ ] GTF DL
      - [ ] genomeSizes

