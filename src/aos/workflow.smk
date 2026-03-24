import yaml
import pandas as pd
from pathlib import Path

BAMSAMPLES, = glob_wildcards(config['bamdir'] + "/{sample}.bam")
CRAMSAMPLES, = glob_wildcards(config['bamdir'] + "/{sample}.cram")

SAMPLES = BAMSAMPLES + CRAMSAMPLES

# if a samplesheet is provided, only include samples that are in the samplesheet.
if config['samplesheet']:
  samplesheet = pd.read_csv(config['samplesheet'], sep='\t')
  for sample in SAMPLES:
    if sample not in samplesheet['sample'].values:
      SAMPLES.remove(sample)

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

# def group_files(wildcards):
#   sigdir = checkpoints.collate_sigresults.get(**wildcards).output[0]
#   Path()

#     for group in os.listdir(checkpoint_output):
#         group_path = os.path.join(checkpoint_output, group)
#         if os.path.isdir(group_path):
#             for f in os.listdir(group_path):
#                 result[group].append(os.path.join(group_path, f))

#     return dict(result)

# if config['motifs']:
#   # do motif analysis.
#   print("motif")

include: "rules/1_peaks.smk"
include: "rules/1_qc.smk"
include: "rules/2_twogroup_de.smk"
include: "rules/2_lrt_de.smk"
include: "rules/2_gp_de.smk"
include: "rules/3_collate_sigresults.smk"
# include: "rules/motifs.smk"
# include: "rules/tobias.smk"
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
    OUTPUTFILES
    # Prep different differential calls for motif / footprinting
