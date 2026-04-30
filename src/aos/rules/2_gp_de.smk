rule gp_diffacc:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'gp/{comparison}/{comparison}_gp_results.tsv',
    norm_matrix = 'gp/{comparison}/{comparison}_counts_norm.tsv',
    acc_pred = 'gp/{comparison}/{comparison}_acc_pred.tsv',
    acc_pred_std = 'gp/{comparison}/{comparison}_acc_pred_std.tsv'
  params:
    ss = config['samplesheet'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison],
    permutations = config['cutoffs']['permutation_iterations'],
    gp_timesteps = config['cutoffs']['gp_timesteps']
  benchmark: "benchmarks/2_gp_{comparison}.txt"
  threads: 50
  conda: "envs/gp.yml"
  script:
    "scripts/gp.py"

rule gp_diffacc_interaction:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'gp/{comparison}/inttest_{comparison}_{interaction}_gp_results.tsv',
    acc_pred = 'gp/{comparison}/inttest_{comparison}_{interaction}_acc_pred.tsv',
    acc_pred_std = 'gp/{comparison}/inttest_{comparison}_{interaction}_acc_pred_std.tsv'
  params:
    ss = config['samplesheet'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison],
    permutations = config['cutoffs']['permutation_iterations'],
    gp_timesteps = config['cutoffs']['gp_timesteps'],
    interaction = lambda wildcards: wildcards.interaction
  benchmark: "benchmarks/2_gp_interaction_{comparison}_{interaction}.txt"
  threads: 50
  conda: 'envs/gp.yml'
  script:
      'scripts/gp_interaction.py'

rule gp_postprocessing:
  input:
    results = 'gp/{comparison}/{comparison}_gp_results.tsv',
    acc_pred = 'gp/{comparison}/{comparison}_acc_pred.tsv',
    acc_pred_std = 'gp/{comparison}/{comparison}_acc_pred_std.tsv'
  output:
    donefile = 'gp/{comparison}/{comparison}_postprocess.done'
  params:
    permutation_cutoff = config['cutoffs']['permutation_cutoff'],
    comp_name = lambda wildcards: wildcards.comparison,
    min_sigpeaks = config['cutoffs']['min_sigpeaks'],
    odir = lambda wildcards: f"gp/{wildcards.comparison}"
  benchmark: "benchmarks/2_gp_postprocess_{comparison}.txt"
  threads: 20
  conda: 'envs/gp.yml'
  script:
    'scripts/gp_postprocess.py'

rule gp_postprocessing_interaction:
  input:
    results = 'gp/{comparison}/inttest_{comparison}_{interaction}_gp_results.tsv',
    acc_pred = 'gp/{comparison}/inttest_{comparison}_{interaction}_acc_pred.tsv',
    acc_pred_std = 'gp/{comparison}/inttest_{comparison}_{interaction}_acc_pred_std.tsv'
  output:
    donefile = 'gp/{comparison}/inttest_{comparison}_{interaction}_postprocess.done'
  params:
    permutation_cutoff = config['cutoffs']['permutation_cutoff'],
    comp_name = lambda wildcards: wildcards.comparison,
    min_sigpeaks = config['cutoffs']['min_sigpeaks'],
    y_pred = lambda wildcards: f"gp/{wildcards.comparison}/inttest_{wildcards.comparison}_{wildcards.interaction}_acc_pred.tsv",
    odir = lambda wildcards: f"gp/{wildcards.comparison}"
  benchmark: "benchmarks/2_gp_interaction_postprocess_{comparison}_{interaction}.txt"
  threads: 20
  conda: 'envs/gp.yml'
  script:
    'scripts/gp_postprocess.py'