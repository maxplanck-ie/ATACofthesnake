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

rule gp_postprocessing:
  input:
    results = 'gp/{comparison}/{comparison}_gp_results.tsv',
  output:
    touch('gp/{comparison}/{comparison}_postprocess.done')
  params:
    permutation_cutoff = config['cutoffs']['permutation_cutoff'],
    comp_name = lambda wildcards: wildcards.comparison,
    min_sigpeaks = config['cutoffs']['min_sigpeaks'],
    y_pred = lambda wildcards: f"gp/{wildcards.comparison}/{wildcards.comparison}_acc_pred.tsv"
    odir = lambda wildcards: f"gp/{wildcards.comparison}"
  threads: 20
  conda: 'envs/gp.yml'
  script:
    'scripts/gp_postprocess.py'

rule gp_postprocessing_interaction:
  input:
    results = 'gp/{comparison}/inttest_{comparison}_{interaction}_gp_results.tsv',
  output:
    touch('gp/{comparison}/inttest_{comparison}_{interaction}_postprocess.done')
  params:
    permutation_cutoff = config['cutoffs']['permutation_cutoff'],
    comp_name = lambda wildcards: wildcards.comparison,
    min_sigpeaks = config['cutoffs']['min_sigpeaks'],
    y_pred = lambda wildcards: f"gp/{wildcards.comparison}/{wildcards.comparison}_acc_pred.tsv"
    odir = lambda wildcards: f"gp/{wildcards.comparison}",
    int = lambda wildcards: wildcards.interaction
  threads: 20
  conda: 'envs/gp.yml'
  script:
    'scripts/gp_postprocess.py'