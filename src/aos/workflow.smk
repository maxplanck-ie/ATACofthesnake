import yaml
import pandas as pd
from pathlib import Path

BAMSAMPLES, = glob_wildcards(config['bamdir'] + "/{sample}.bam")
CRAMSAMPLES, = glob_wildcards(config['bamdir'] + "/{sample}.cram")

SAMPLES = BAMSAMPLES + CRAMSAMPLES

# if a samplesheet is provided, only include samples that are in the samplesheet.
samplesheet = None
if config['samplesheet']:
  samplesheet = pd.read_csv(config['samplesheet'], sep='\t')
  allowed_samples = set(samplesheet['sample'].values)
  SAMPLES = [s for s in SAMPLES if s in allowed_samples]

# Define output for comparisons.
def define_comparison_output():
  outputfiles = []
  sigresults = []
  if config['comparison']:
    with open(config['comparison'], 'r') as f:
      config['comparison'] = yaml.safe_load(f)
      for comp in config['comparison'].keys():
        match config['comparison'][comp]['type']:
          case 'twogroup':
            outputfiles.extend(
              [
                f"twogroup/{comp}/{comp}_diffacc_edgeR.tsv",
                f"twogroup/{comp}/relevant_samples.txt",
                f"twogroup/{comp}/{comp}_maplot.png",
                f"twogroup/{comp}/{comp}_heatmap.done"
              ]
            )
            sigresults.append(f"twogroup/{comp}/{comp}_heatmap.done")
          case 'lrt':
            outputfiles.extend(
              [
                f"lrt/{comp}/{comp}_lrt_edgeR.tsv",
                f"lrt/{comp}/relevant_samples.txt",
                f"lrt/{comp}/{comp}_heatmap.done"
              ]
            )
            sigresults.append(f"lrt/{comp}/{comp}_heatmap.done")
          case 'timecourse':
            outputfiles.extend(
              [
                f"gp/{comp}/{comp}_gp_results.tsv",
                f"gp/{comp}/{comp}_counts_norm.tsv",
                f"gp/{comp}/{comp}_postprocess.done"
              ]
            )
            sigresults.append(f"gp/{comp}/{comp}_postprocess.done")
            if 'interaction' in config['comparison'][comp]:
              if isinstance(config['comparison'][comp]['interaction'], list):
                for interaction in config['comparison'][comp]['interaction']:
                  outputfiles.extend(
                    [
                        f"gp/{comp}/inttest_{comp}_{interaction}_gp_results.tsv",
                        f"gp/{comp}/inttest_{comp}_{interaction}_postprocess.done"
                    ]
                  )
                  sigresults.append(f"gp/{comp}/inttest_{comp}_{interaction}_postprocess.done")
              elif isinstance(config['comparison'][comp]['interaction'], str):
                interaction = config['comparison'][comp]['interaction']
                outputfiles.extend(
                  [
                      f"gp/{comp}/inttest_{comp}_{interaction}_gp_results.tsv",
                      f"gp/{comp}/inttest_{comp}_{interaction}_postprocess.done"
                  ]
                )
                sigresults.append(f"gp/{comp}/inttest_{comp}_{interaction}_postprocess.done")
  return (outputfiles, sigresults)

OUTPUTFILES, SIGRESULTS = define_comparison_output()

# Define outputs for motif enrichment.
def get_ames(wildcards):
  if config['motif']:
    checkpoints.collate_sigresults.get(**wildcards)
    motif_comps, motif_samples = glob_wildcards("motifs/{motif_comp}/{sample}.bed")
    return expand(
        "motifs/{motif_comp}/{motif_sample}_ame/ame.tsv",
        zip,
        motif_comp=motif_comps,
        motif_sample=motif_samples
    )
  return []

def get_ame_plots(wildcards):
  # For rule all — one output per unique motif_comp
  if config['motif']:
    checkpoints.collate_sigresults.get(**wildcards)
    motif_comps, _ = glob_wildcards("motifs/{motif_comp}/{sample}.bed")
    return expand(
      "motifs/{motif_comp}/motif_enrichment.parsed",
      motif_comp=set(motif_comps)
    )
  return []

# Define outputs for footprinting
FPDIC = {}
if samplesheet is not None:
  group_cols = samplesheet.columns.drop("sample").tolist()
  by = group_cols[0] if len(group_cols) == 1 else group_cols
  groups = samplesheet.groupby(by).groups
  for group in groups.keys():
    # get a groupname by joining the group values with underscores
    if isinstance(group, str):
      groupname = group
    else:
      groupname = "_".join([str(g) for g in group])
    FPDIC[groupname] = samplesheet.loc[groups[group], "sample"].tolist()

def get_all_plotame(wildcards):
  checkpoints.collate_sigresults.get()
  motif_comps, _ = glob_wildcards("motifs/{motif_comp}/{sample}.bed")
  return expand(
    'motifs/{motif_comp}/motif_enrichment.parsed', motif_comp=set(motif_comps)
  )

def get_motifs_for_fp(wildcards):
  if config['motif']:
    checkpoints.all_plotame_done.get()
    enr_found, = glob_wildcards('motifs/{fpmotif}/motif_enrichment.png')
    fpmotifs, group_fnas = glob_wildcards("motifs/{fpmotif}/{group_fna}.fna")
    valid_pairs = [
        (motif, fna)
        for motif, fna in zip(fpmotifs, group_fnas)
        if motif in enr_found and not fna.endswith('_bg')
    ]
    if not valid_pairs:
        return []
    valid_motifs, valid_fnas = zip(*valid_pairs)
    
    return expand(
      'footprints/fimo/{fpmotif}/{group_fna}/fimo.tsv',
      zip,
      fpmotif=valid_motifs,
      group_fna=valid_fnas
    )
  return []

def get_motif_for_aggplot(wildcards):
  if config['motif']:
    checkpoints.parse_fimo.get(**wildcards)
    comps, motifs, _ = glob_wildcards("footprints/plotaggregate/{comp}/bedfiles/{motif}-{group}.bed")
    return expand(
        "footprints/plotaggregate/{comp}/{motif}.pdf",
        zip,
        comp=comps,
        motif=motifs
    )
  return []

include: "rules/1_peaks.smk"
include: "rules/1_qc.smk"
include: "rules/2_twogroup_de.smk"
include: "rules/2_lrt_de.smk"
include: "rules/2_gp_de.smk"
include: "rules/3_collate_sigresults.smk"
include: "rules/3_motifs.smk"
include: "rules/3_tobias.smk"

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
    OUTPUTFILES,
    # motif enrichment - AME
    get_ames,
    get_ame_plots,
    # footprinting - TOBIAS
    expand('footprints/scores/{fp_group}_scores.bw', fp_group=FPDIC.keys()),
    get_motifs_for_fp,
    get_motif_for_aggplot,
