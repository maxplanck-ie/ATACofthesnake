rule gp_diffacc:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'gp/{comparison}/{comparison}_gp_results.tsv',
  params:
    samplesheet = config['samplesheet'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison],
    outputfolder = lambda wildcards: f"gp/{wildcards.comparison}",
    comparison_name = lambda wildcards: wildcards.comparison,
    permutations = config['cutoffs']['permutation_iterations'],
    permutation_cutoff = config['cutoffs']['permutation_cutoff'],
    gp_timesteps = config['cutoffs']['gp_timesteps'],
    sig_clustered = lambda wildcards: f"gp/{wildcards.comparison}/{wildcards.comparison}_gp_sig_clustered.tsv"
  threads: 20
  conda: "envs/gp.yml"
  script:
    "scripts/gp.py"
