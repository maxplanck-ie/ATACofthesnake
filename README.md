# ATACofthesnake

Downstream processing of ATAC data, including QC's and differential accessibility. Starting point are deduplicated bam files, could be obtained from snakePipes (https://github.com/maxplanck-ie/snakepipes).


  - Installation

  set up the environment:  
>  conda env create -f CondaEnv.yaml  
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
 - [ ] deal with illegal characters ?
 - [ ] Standardize logging
 - [x] scaleFactors
 - [x] bamCoverage
 - [x] computeMatrix
 - [x] deeptools: TSS enrichment
 - [x] deeptools: correlations
 - [ ] DESeq2
 - [ ] annotate
 - [ ] slurm submission
 - [ ] produce report
 - [ ] incorporate runID if more than 1 comparison per batch.
 - [ ] Support multiple genomes
      - [ ] blacklist DL + merge
      - [ ] GTF DL
      - [ ] genomeSizes

