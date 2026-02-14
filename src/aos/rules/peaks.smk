from aos.helper import peak_boundaries
from aos.helper import PCA_colors

rule lnBams:
  output: 
    'input/{sample}.bam'
  threads: 1
  run:
    if wildcards.sample in BAMSAMPLES:
      shell('ln -s {config[bamdir]}/{wildcards.sample}.bam {output}')
    elif wildcards.sample in CRAMSAMPLES:
      shell('samtools view -T {config[fna]} -b -h {config[bamdir]}/{wildcards.sample}.cram > {output}')

rule bamIx:
  input:
    'input/{sample}.bam'
  output:
    'input/{sample}.bam.bai'
  conda: "envs/seqtools.yml"
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
  conda: "envs/deeptools.yml"
  threads: 10
  params:
    rar = f"--blackListFileName {config['rar']}",
    size = f"--maxFragmentLength {config['fragsize']} --minFragmentLength 0"
  shell:'''
  alignmentSieve --bam {input.bam} --outFile {output.shortb} --filterMetrics {output.qc} -p {threads} {params.rar} {params.size}
  '''

rule ixSieve:
  input:
    'sieve/{sample}.bam'
  output:
    'sieve/{sample}.bam.bai'
  conda: "envs/seqtools.yml"
  threads: 5
  shell:'''
  samtools index -@ {threads} {input}
  '''


rule bamtobed:
  input:
    b = 'sieve/{sample}.bam',
    i = 'sieve/{sample}.bam.bai'
  output:
    temp('sieve/{sample}.bed')
  conda: "envs/seqtools.yml"
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
    gsize = config['ESS'],
    outname = lambda wildcards: wildcards.sample,
    rar = config['rar'],
    outdir = 'peaks'
  conda: "envs/seqtools.yml"
  threads: 1
  shell:'''
  macs3 callpeak -t {input} \
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
    expand('peaks/{sample}_peaks.narrowPeak', sample=SAMPLES)
  output:
    temp('peakset/peaks.cat.bed')
  threads:1
  conda: "envs/seqtools.yml"
  shell:'''
  cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
  '''

rule peakbounds:
  input:
    'peakset/peaks.cat.bed'
  output:
    'peakset/peaks.bed'
  params:
    fna = config['fna'],
    peakset = config['peakset']
  threads: 1
  run:
    peak_boundaries(input[0], params.fna, params.peakset, output[0])

rule uropa:
  input:
    'peakset/peaks.bed'
  output:
    'peakset/peaks_uropa_finalhits.txt'
  params:
    gtf = config['gtf_sorted'],
    prefix = "peaks_uropa",
    outdir = 'peakset',
    upstream = config['uro'][1],
    downstream = config['uro'][2],
    feature = config['uro'][0]
  conda: "envs/seqtools.yml"
  threads: 5
  shell:'''
  uropa -b {input} \
    -g {params.gtf} \
    --summary --feature {params.feature} \
    --distance {params.upstream} {params.downstream} \
    --internals 1 -p {params.prefix} \
    -o {params.outdir} -t {threads} --show-attributes gene_id transcript_id gene_name
  '''

rule countmatrix:
  input:
    peaks = 'peakset/peaks.bed',
    bams = expand('sieve/{sample}.bam', sample=SAMPLES)
  output:
    tsv = 'peakset/counts.tsv',
    npz = 'peakset/counts.npz'
  params:
    rar = '-bl {}'.format(
        config['rar']
    )
  threads: 40
  conda: "envs/deeptools.yml"
  shell:'''
  multiBamSummary BED-file --BED {input.peaks} {params.rar} \
    -p {threads} --outRawCounts {output.tsv} -o {output.npz} \
    -b {input.bams}
  sed -i "s/'//g" {output.tsv}
  sed -i 's/\\.bam//g' {output.tsv}
  '''

rule multibigwigsum:
  input:
    samples = expand('bw/{sample}.scalefac.bw', sample=SAMPLES),
    peaks = 'peakset/peaks.bed'
  output:
    'peakset/counts.bw.npz'
  threads: 1
  conda: "envs/deeptools.yml"
  shell:'''
  multiBigwigSummary BED-file --BED {input.peaks} -o {output} -b {input.samples}
  '''

rule plotPCA:
  input:
    peakset = 'peakset/counts.bw.npz'
  output:
    'figures/PCA.png'
  threads: 1
  params:
    colstr = PCA_colors(config['samplesheet'], SAMPLES)
  conda: "envs/deeptools.yml"
  shell:'''
  plotPCA --corData {input.peakset} -o {output} --transpose --ntop 5000 {params.colstr}
  '''