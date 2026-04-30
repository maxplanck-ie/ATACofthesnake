rule subsample_bams:
  input:
    bams = lambda wildcards: expand('input/{sample}.bam', sample=FPDIC[wildcards.fp_group])
  output:
    finbam = 'footprints/bams/{fp_group}.bam',
    finbai = 'footprints/bams/{fp_group}.bam.bai',
  conda:
    'envs/seqtools.yml'
  benchmark: "benchmarks/3_subsample-bams_{fp_group}.txt"
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
      if fraction == 1:
        shell(f"cp {bam} {obam}")
        shell(f"cp {bam}.bai {obam}.bai")
      else:
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
  benchmark: "benchmarks/3_tobias-ataccorrect_{fp_group}.txt"
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
  benchmark: "benchmarks/3_tobias-scorebigwig_{fp_group}.txt"
  shell:'''
  TOBIAS FootprintScores --signal {input.signal} \
    --regions {input.peaks} \
    --output {output.bw} \
    --cores {threads}
  '''

rule fimo_footprints:
  input:
    ame_done = '.snakemake_flags/all_plotame_done',
    motifs = 'motifs/{fpmotif}/{fpmotif}_consensus_motifs.meme',
    fna = 'motifs/{fpmotif}/{group_fna}.fna'
  output:
    'footprints/fimo/{fpmotif}/{group_fna}/fimo.tsv'
  params:
    of = lambda wildcards: f"footprints/fimo/{wildcards.fpmotif}/{wildcards.group_fna}"
  conda:
    'envs/meme.yml'
  benchmark: "benchmarks/3_fimo-footprints_{fpmotif}_{group_fna}.txt"
  shell:'''
  fimo --oc {params.of} {input.motifs} {input.fna}
  '''

checkpoint parse_fimo:
  input:
    tsvs = get_motifs_for_fp
  output:
    touch('.snakemake_flags/parse_fimo_done')
  script:
    'scripts/parse_fimo.py'

rule plot_aggregate:
    input:
        bws = expand('footprints/ataccorrect/{fp_group}_corrected.bw', fp_group=FPDIC.keys()),
    output:
        pdf = "footprints/plotaggregate/{comp}/{motif}.pdf",
        txt = "footprints/plotaggregate/{comp}/{motif}.txt"
    params:
        bedfiles = lambda wc: ' '.join(
            str(p) for p in 
            Path(f"footprints/plotaggregate/{wc.comp}/bedfiles").glob(f"{wc.motif}-*.bed")
        ),
        motname = lambda wc: pd.read_table(
            f'footprints/plotaggregate/{wc.comp}/{wc.comp}_motif_names.tsv', index_col=0
        ).loc[wc.motif, 'motif_alt_id']
    conda:
        'envs/tobias.yml'
    benchmark: "benchmarks/3_plot-aggregate_{comp}_{motif}.txt"
    shell: """
        TOBIAS PlotAggregate --TFBS {params.bedfiles} \
          --signals {input.bws} \
          --output {output.pdf} \
          --share_y both --plot_boundaries --signal-on-x \
          --title "{params.motname}" \
          --output-txt {output.txt}
    """