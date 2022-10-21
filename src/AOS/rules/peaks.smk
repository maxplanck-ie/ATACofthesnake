rule lnBams:
  input: 
    os.path.join(config['dirs']['bamdir'], "{sample}.bam")
  output: 
    'input/{sample}.bam'
  threads: 1
  shell:'''
  ln -t input -s {input}
  '''

rule bamIx:
  input:
    'input/{sample}.bam'
  output:
    'input/{sample}.bam.bai'
  conda: config['envs']['seqtools']
  threads: 5
  shell:'''
  samtools index -@ {threads} {input}
  '''

rule alSieve:
  input:
    bam = 'input/{sample}.bam',
    bai = 'input/{sample}.bam.bai'
  output:
    shortb = 'sieve/{sample}.bam',
    qc = temp('qc/{sample}_sieve.txt')
  conda: config['envs']['seqtools']
  threads: 10
  params:
    rar = '--blackListFileName {}'.format(config['files']['readattractingregions']),
    size = '--maxFragmentLength {} --minFragmentLength 0'.format(config['vars']['fragsize'])
  shell:'''
  alignmentSieve --bam {input.bam} --outFile {output.shortb} --filterMetrics {output.qc} -p {threads} {params.rar} {params.size}
  '''

rule ixSieve:
  input:
    'sieve/{sample}.bam'
  output:
    'sieve/{sample}.bam.bai'
  conda: config['envs']['seqtools']
  threads: 5
  shell:'''
  samtools index -@ {threads} {input}
  '''


rule bamtobed:
  input:
    b = 'sieve/{sample}.bam',
    i = 'sieve/{sample}.bam.bai'
  output:
    'sieve/{sample}.bed'
  conda: config['envs']['seqtools']
  threads: 1
  shell:'''
  bedtools bamtobed -i {input.b} > {output}
  '''

rule peaks:
  input:
    'sieve/{sample}.bed'
  output:
    'peaks/{sample}_peaks.narrowPeak'
  params:
    gsize = config['vars']['genomesize'],
    outname = lambda wildcards: wildcards.sample,
    rar = config['files']['readattractingregions'],
    outdir = 'peaks'
  conda: config['envs']['seqtools']
  threads: 1
  shell:'''
  macs2 callpeak -t {input} \
    -f BED \
    --nomodel --shift -75 \
    --extsize 150 \
    -g {params.gsize} \
    -n {params.outname} \
    -q 0.01 \
    --outdir {params.outdir} \
    --keep-dup all
  '''
 
rule peakmerg:
  input:
    expand('peaks/{sample}_peaks.narrowPeak', sample=config['samples'])
  output:
    'peakset/peaks.bed'
  threads:1
  conda: config['envs']['seqtools']
  shell:'''
  cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
  '''

rule uropa:
  input:
    'peakset/peaks.bed'
  output:
    'peakset/peaks_uropa_finalhits.txt'
  params:
    gtf = config['files']['gtf'],
    prefix = "peaks_uropa",
    outdir = 'peakset'
  conda: config['envs']['seqtools']
  threads: 5
  shell:'''
  uropa -b {input} \
    -g {params.gtf} \
    --summary --feature transcript \
    --distance 20000 10000 \
    --internals 1 -p {params.prefix} \
    -o {params.outdir} -t {threads} --show-attributes gene_id transcript_id gene_name gene_type transcript_type
  '''

rule countmatrix:
  input:
    peaks = 'peakset/peaks.bed',
    bams = expand('sieve/{sample}.bam', sample=config['samples'])
  output:
    tsv = 'peakset/counts.tsv',
    npz = 'peakset/counts.npz'
  params:
    rar = '-bl {}'.format(
        config['files']['readattractingregions']
    )
  threads: 10
  conda: config['envs']['seqtools']
  shell:'''
  multiBamSummary BED-file --BED {input.peaks} {params.rar} \
    -p {threads} --outRawCounts {output.tsv} -o {output.npz} \
    -b {input.bams}
  sed -i "s/'//g" {output.tsv}
  sed -i 's/\.bam//g' {output.tsv}
  '''