rule clustermotifs:
  input:
    config['files']['motif']
  output:
    'motifs_clustered/motifs.meme'
  conda: config['envs']['tobias']
  shell:'''
  TOBIAS ClusterMotifs -m {input} -t 0.4 -a meme -p clusteredmotifs -o 'motifs_clustered' --dist_method seqcor
  '''