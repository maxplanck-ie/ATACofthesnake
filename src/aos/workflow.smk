import yaml

SAMPLES, = glob_wildcards(config['bamdir'] + "/{sample}.bam")

include: "rules/peaks.smk"
include: "rules/qc.smk"
include: "rules/DE.smk"
include: "rules/motifs.smk"
include: "rules/tobias.smk"


def compzip(comp):
    comps = [c for c in comp for _ in comp[c]]
    grs = [gr for c in comp for gr in comp[c]]
    return comps, grs

def define_comparison_output():
  outputfiles = []
  if config['comparison']:
    with open(config['comparison'], 'r') as f:
      config['comparison'] = yaml.safe_load(f)
    # define list of comparisons & list per groups to easily zip later on.
    comps, grs = compzip(config['comparison'])
    outputfiles.extend(expand('{comp}/{comp}_maplot.png', comp=config['comparison'].keys()))
    outputfiles.extend(expand('{comp}/diffpeaks_{gr}.bed', zip, comp=comps, gr=grs))
    outputfiles.extend(expand('{comp}/diffpeaks.png', comp=config['comparison'].keys()))
    
    if config['motif']:
      outputfiles.append('motifs_clustered/clusteredmotifs_consensus_motifs.meme')
      outputfiles.extend(
        expand(
          '{comp}/motif_{gr}/ame.html',
          zip,comp=comps,
          gr=grs
        )
      )
      outputfiles.extend(expand('{comp}/shuffled_motif_{gr}/ame.html',zip,comp=comps,gr=grs))
  return (outputfiles)

rule all:
  input:
    'peakset/peaks_uropa_finalhits.txt',
    'figures/PCA.png',
    expand('bw/{sample}.scalefac.bw', sample=SAMPLES),
    'figures/mitofraction.png',
    'qc/fragsize.tsv',
    'figures/alignmentsieve.png',
    'figures/fragmentsizes.png',
    'figures/fripscores.png',
    define_comparison_output()