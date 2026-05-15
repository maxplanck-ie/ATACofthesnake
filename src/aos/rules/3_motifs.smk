rule clustermotifs:
  output:
    'motifs/clusteredmotifs_consensus_motifs.meme'
  params:
    motiffile = config['motif']
  conda: "envs/tobias.yml"
  threads: 2
  benchmark: "benchmarks/3_clustermotifs.txt"
  shell:'''
  TOBIAS ClusterMotifs -m {params.motiffile} -t 0.4 -a meme -p clusteredmotifs -o 'motifs' --dist_method seqcor
  '''

rule bed2fna:
  input:
    "motifs/{motif_comp}/{motif_sample}.bed"
  output:
    "motifs/{motif_comp}/{motif_sample}.fna"
  threads: 1
  conda: "envs/seqtools.yml"
  params:
    fna = config['fna']
  benchmark: "benchmarks/3_bed2fna_{motif_comp}_{motif_sample}.txt"
  shell:'''
  bedtools getfasta -fi {params.fna} -bed {input} > {output}
  '''

rule ame:
  input:
    motifs = 'motifs/clusteredmotifs_consensus_motifs.meme',
    fna = 'motifs/{motif_comp}/{motif_sample}.fna',
  output:
    tsv = 'motifs/{motif_comp}/{motif_sample}_ame/ame.tsv',
  conda: "envs/meme.yml"
  benchmark: "benchmarks/3_ame_{motif_comp}_{motif_sample}.txt"
  run:
    from pathlib import Path
    fnapath = Path(input.fna)
    bg_arg = "--control --shuffle--"
    shell(f'''
      ame --oc motifs/{wildcards.motif_comp}/{wildcards.motif_sample}_ame {bg_arg} {input.fna} {input.motifs}
    ''')

rule plotame:
  input:
    ames = get_ames,
    motiffile = 'motifs/clusteredmotifs_consensus_motifs.meme'
  output:
    out = touch('motifs/{motif_comp}/motif_enrichment.parsed')
  params:
    motifcomp = lambda wildcards: wildcards.motif_comp,
    fdr_cutoff = config['cutoffs']['fdr_cutoff']
  conda:
    'envs/gimmemotifs.yml'
  threads: 2
  benchmark: "benchmarks/3_plotame_{motif_comp}.txt"
  script:
    'scripts/plot_ame.py'

checkpoint all_plotame_done:
    input: get_all_plotame
    output: touch('.snakemake_flags/all_plotame_done')
