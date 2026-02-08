def bg_from_gr(compdic, comp, gr):
  for bg in compdic[comp]:
    if bg != gr:
      return (
        '{}/diffpeaks_{}.fna'.format(comp, bg)
      )
 
rule clustermotifs:
  output:
    'motifs_clustered/clusteredmotifs_consensus_motifs.meme'
  params:
    motiffile = config['motif']
  conda: "envs/tobias.yml"
  threads: 10
  shell:'''
  TOBIAS ClusterMotifs -m {params.motiffile} -t 0.4 -a meme -p clusteredmotifs -o 'motifs_clustered' --dist_method seqcor
  '''

rule bed2fna:
  input:
    '{comparison}/diffpeaks_{gr}.bed'
  output:
    '{comparison}/diffpeaks_{gr}.fna'
  threads: 1
  conda: "envs/seqtools.yml"
  params:
    fna = config['fna']
  shell:'''
  bedtools getfasta -fi {params.fna} -bed {input} > {output}
  '''

rule ame:
  input:
    motifs = 'motifs_clustered/clusteredmotifs_consensus_motifs.meme',
    fna = '{comparison}/diffpeaks_{gr}.fna',
    fnabg = lambda wildcards: bg_from_gr(
      config['comparison'],
      wildcards.comparison,
      wildcards.gr
    )
  output:
    html = '{comparison}/motif_{gr}/ame.html'
  params:
    of = lambda wildcards: wildcards.comparison + '/motif_' + wildcards.gr
  conda: "envs/meme.yml"
  shell:'''
  ame --oc {params.of} --control {input.fnabg} {input.fna} {input.motifs}
  '''

rule ame_shuffled:
  input:
    motifs = 'motifs_clustered/clusteredmotifs_consensus_motifs.meme',
    fna = '{comparison}/diffpeaks_{gr}.fna'
  output:
    html = '{comparison}/shuffled_motif_{gr}/ame.html'
  params:
    of = lambda wildcards: wildcards.comparison + '/shuffled_motif_' + wildcards.gr
  conda: "envs/meme.yml"
  shell:'''
  ame --oc {params.of} --control --shuffle-- {input.fna} {input.motifs}
  '''
