from aos.helper import merge_idx, merge_sieve
from aos.helper import plotfragsize, plotfrip, plotixs, plotsieve

rule ixStat:
  input:
    bam = 'input/{sample}.bam',
    bai = 'input/{sample}.bam.bai'
  output:
    temp('qc/{sample}_ix.tsv')
  conda: "envs/seqtools.yml"
  threads: 1
  shell:'''
  set +o pipefail;
  samtools idxstats {input.bam} | cut -f1,3 > {output}
  '''

rule mergeixStat:
  input:
    expand('qc/{sample}_ix.tsv', sample=SAMPLES)
  output:
    'qc/ixstat.tsv'
  run:
    merge_idx(input, output) 

rule mergesieve:
  input:
    expand('qc/{sample}_sieve.txt', sample=SAMPLES)
  output:
    'qc/sieve.tsv'
  run:
    merge_sieve(input, output)

rule fragsize:
  input:
    bams = expand('input/{sample}.bam', sample=SAMPLES),
    bais = expand('input/{sample}.bam.bai', sample=SAMPLES)
  output:
    'qc/fragsize.tsv'
  threads: 10
  conda: "envs/deeptools.yml"
  shell:'''
  bamPEFragmentSize -b {input.bams} \
    -p {threads} \
    --outRawFragmentLengths {output}
  '''

rule scalefactors:
  input:
    'peakset/counts.tsv'
  output:
    'peakset/scalefactors.txt'
  threads: 1
  conda: "envs/seqtools.yml"
  script:
    "scripts/scalefactors.R"

rule bigwigs:
  input:
    scalefactors = 'peakset/scalefactors.txt',
    bam = 'sieve/{sample}.bam',
    bai = 'sieve/{sample}.bam.bai'
  output:
    bigwigsf = 'bw/{sample}.scalefac.bw',
    bigwigrpkm = 'bw/{sample}.RPKM.bw'
  params:
    rar = config['rar'],
    sample = "{sample}"
  threads: 10
  conda: "envs/deeptools.yml"
  shell:'''
  SCALEFAC=$(grep "^{params.sample} " {input.scalefactors} | cut -f2 -d ' ')
  bamCoverage --scaleFactor $SCALEFAC \
    -b {input.bam} -o {output.bigwigsf} \
    -p {threads} -bs 1 -bl {params.rar}
  bamCoverage --normalizeUsing RPKM \
    -b {input.bam} -o {output.bigwigrpkm} \
    -p {threads} -bs 1 -bl {params.rar}
  '''

rule frips:
  input:
    bam = 'sieve/{sample}.bam',
    bai = 'sieve/{sample}.bam.bai',
    peaks = 'peakset/peaks.bed'
  output:
    temp('qc/{sample}.frip.txt')
  params:
    sample = '{sample}'
  conda: "envs/seqtools.yml"
  shell:'''
  mapped=$(samtools view -c -F 4 {input.bam})
  peakreads=$(samtools view -c -F 4 -L {input.peaks} {input.bam})
  frip=$(bc -l <<< $peakreads/$mapped)
  printf "%s\t%5.3f\n" {params.sample} $frip > {output}
  '''

rule fripcombine:
  input:
    expand('qc/{sample}.frip.txt', sample=SAMPLES)
  output:
    'qc/fripscores.txt'
  shell:'''
  cat {input} > {output}
  '''

rule plotfragsize:
  input:
    fs = 'qc/fragsize.tsv'
  output:
    of = 'figures/fragmentsizes.png'  
  run:
    plotfragsize(input.fs, output.of)

rule plotfrip:
  input:
    fs = 'qc/fripscores.txt'
  output:
    of = 'figures/fripscores.png'
  run:
    plotfrip(input.fs, output.of)

rule plotixs:
  input:
    'qc/ixstat.tsv'
  output:
    'figures/mitofraction.png'
  params:
    mito = config['mitostring']
  run:
    plotixs(input[0], params[0])

rule plotsieve:
  input:
   'qc/sieve.tsv'
  output:
    'figures/alignmentsieve.png'
  run:
    plotsieve(input[0])
