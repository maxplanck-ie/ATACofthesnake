rule subsample_bams:
  input:
    bams = lambda wildcards: expand('input/{sample}.bam', sample=FPDIC[wildcards.fp_group])
  output:
    finbam = 'footprints/bams/{fp_group}.bam',
    finbai = 'footprints/bams/{fp_group}.bam.bai',
  conda:
    'envs/seqtools.yml'
  threads: 10
  run:
    import pandas as pd
    from pathlib import Path

    Path(output.finbam).parent.mkdir(parents=True, exist_ok=True)

    countdic = {}
    for bam in input.bams:
      count = int(
        shell(
          f"samtools view -@ {threads} -c -F 2308 {bam}",
          iterable=False,
          read=True
        ).strip()
      )
      countdic[Path(bam).stem] = count
    target = min(countdic.values())
    print(f"Subsampling group {wildcards.fp_group} to {target} reads per sample.")
    temp_bams = []
    for bam in input.bams:
      sample = Path(bam).stem
      obam = f"footprints/bams/{sample}.bam"
      fraction = target / countdic[sample]
      fracstr = f"{fraction:.6f}".split(".")[1]
      fracstr = f"1337.{fracstr}"
      shell(f"samtools view -@ {threads} -b -s {fracstr} {bam} -o {obam}")
      shell(f"samtools index -@ {threads} {obam}")
      temp_bams.append(obam)
        
    # Merge all subsamples bams together.
    shell(f"samtools merge -@ {threads} {output.finbam} {' '.join(temp_bams)}")
    shell(f"samtools index -@ {threads} {output.finbam}")
    
    # Cleanup temp files.
    for bam in input.bams:
      sample = Path(bam).stem
      obam = f"footprints/bams/{sample}.bam"
      obai = f"footprints/bams/{sample}.bam.bai"
      Path(obam).unlink()
      Path(obai).unlink()

rule tobias_ataccorrect:
  input:
    bam = 'footprints/bams/{fp_group}.bam',
    peaks = 'peakset/peaks.bed'
  output:
    'footprints/ataccorrect/{fp_group}_corrected.bw'
  params:
    rar = config['rar'],
    fna = config['fna'],
    mitostr = config['mitostring']
  threads: 10
  conda:
    'envs/tobias.yml'
  shell:'''
  TOBIAS ATACorrect --bam {input.bam} \
    --genome {params.fna} \
    --peaks {input.peaks} \
    --blacklist {params.rar} \
    --outdir footprints/ataccorrect \
    --cores {threads} \
    --drop-chroms {params.mitostr}
  '''

rule tobias_scorebigwig:
  input:
    signal = 'footprints/ataccorrect/{fp_group}_corrected.bw',
    peaks = 'peakset/peaks.bed'
  output:
    bw = 'footprints/scores/{fp_group}_scores.bw'
  params:
    rar = config['rar'],
    fna = config['fna'],
    mitostr = config['mitostring']
  threads: 10
  conda:
    'envs/tobias.yml'
  shell:'''
  TOBIAS FootprintScores --signal {input.signal} \
    --regions {input.peaks} \
    --output {output.bw} \
    --cores {threads}
  '''