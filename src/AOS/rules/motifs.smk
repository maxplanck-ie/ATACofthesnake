def bg_from_gr(compdic, comp, gr):
  for bg in compdic[comp]:
    if bg != gr:
      return (
        '{}/diffpeaks_{}.fna'.format(comp, bg)
      )
 
rule clustermotifs:
  input:
    config['files']['motif']
  output:
    'motifs_clustered/clusteredmotifs_consensus_motifs.meme'
  conda: config['envs']['tobias']
  shell:'''
  TOBIAS ClusterMotifs -m {input} -t 0.4 -a meme -p clusteredmotifs -o 'motifs_clustered' --dist_method seqcor
  '''

rule bed2fna:
  input:
    '{comparison}/diffpeaks_{gr}.bed'
  output:
    '{comparison}/diffpeaks_{gr}.fna'
  threads: 1
  conda: config['envs']['seqtools']
  params:
    fna = config['files']['fna']
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
  conda: config['envs']['meme']
  shell:'''
  ame --oc {params.of} --control {input.fnabg} {input.fna} {input.motifs}
  '''