rule gp_diffacc:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'gp/{comparison}/{comparison}_gp_results.tsv',
    norm_matrix = 'gp/{comparison}/{comparison}_counts_norm.tsv'
  params:
    ss = config['samplesheet'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison],
    permutations = config['cutoffs']['permutation_iterations'],
    gp_timesteps = config['cutoffs']['gp_timesteps']
  threads: 50
  conda: "envs/gp.yml"
  script:
    "scripts/gp.py"

rule gp_diffacc_interaction:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'gp/{comparison}/inttest_{comparison}_{interaction}_gp_results.tsv'
  params:
    ss = config['samplesheet'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison],
    permutations = config['cutoffs']['permutation_iterations'],
    gp_timesteps = config['cutoffs']['gp_timesteps'],
  threads: 50
  conda: 'envs/gp.yml'
  script:
      'scripts/gp_interaction.py'
