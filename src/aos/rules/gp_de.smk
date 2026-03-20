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

# rule gp_diffacc_interaction:
#   input:
#     mat = 'peakset/counts.tsv'
#   output:
#     table = 'gp/{comparison}/{comparison}_{interaction}_gp_results.tsv'
#   params:
#     samplesheet = config['samplesheet'],
#     interaction_covar = lambda wildcards: wildcards.interaction,
#     comparison = lambda wildcards: config['comparison'][wildcards.comparison],
#     outputfolder = lambda wildcards: f"gp/{wildcards.comparison}",
#     comparison_name = lambda wildcards: wildcards.comparison,
#     permutations = config['cutoffs']['permutation_iterations'],
#     permutation_cutoff = config['cutoffs']['permutation_cutoff'],
#     gp_timesteps = config['cutoffs']['gp_timesteps'],
#   threads: 20
#   conda: 'envs/gp.yml'
#   script:
#       'scripts/gp_interaction.py'
