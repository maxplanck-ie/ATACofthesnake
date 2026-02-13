import yaml
import pandas as pd

BAMSAMPLES, = glob_wildcards(config['bamdir'] + "/{sample}.bam")
CRAMSAMPLES, = glob_wildcards(config['bamdir'] + "/{sample}.cram")

SAMPLES = BAMSAMPLES + CRAMSAMPLES

# if a samplesheet is provided, only include samples that are in the samplesheet.
if config['samplesheet']:
  samplesheet = pd.read_csv(config['samplesheet'], sep='\t')
  for sample in SAMPLES:
    if sample not in samplesheet['sample'].values:
      SAMPLES.remove(sample)
  

include: "rules/peaks.smk"
include: "rules/qc.smk"
include: "rules/twogroup_de.smk"
# include: "rules/motifs.smk"
# include: "rules/tobias.smk"



# def define_comparison_output():
#   outputfiles = []
#   if config['comparison']:
#     with open(config['comparison'], 'r') as f:
#       config['comparison'] = yaml.safe_load(f)
#     # define list of comparisons & list per groups to easily zip later on.
#     comps, grs = compzip(config['comparison'])
#     outputfiles.extend(expand('{comp}/{comp}_maplot.png', comp=config['comparison'].keys()))
#     outputfiles.extend(expand('{comp}/diffpeaks_{gr}.bed', zip, comp=comps, gr=grs))
#     outputfiles.extend(expand('{comp}/diffpeaks.png', comp=config['comparison'].keys()))
    
#     # if config['motif']:
#     #   outputfiles.append('motifs_clustered/clusteredmotifs_consensus_motifs.meme')
#     #   outputfiles.extend(
#     #     expand(
#     #       '{comp}/motif_{gr}/ame.html',
#     #       zip,comp=comps,
#     #       gr=grs
#     #     )
#     #   )
#     #   outputfiles.extend(expand('{comp}/shuffled_motif_{gr}/ame.html',zip,comp=comps,gr=grs))
#   return (outputfiles)

rule all:
  input:
    # Default output
    'peakset/peaks_uropa_finalhits.txt',
    'figures/PCA.png',
    expand('bw/{sample}.scalefac.bw', sample=SAMPLES),
    'figures/mitofraction.png',
    'qc/fragsize.tsv',
    'figures/alignmentsieve.png',
    'figures/fragmentsizes.png',
    'figures/fripscores.png',
    # DE output
    define_comparison_output()