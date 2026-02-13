# Create a list containing:
from aos.helper import maplot, tsv_to_bed
import os

rule twogroup_diffacc:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'twogroup/{comparison}/{comparison}_diffacc_edgeR.tsv',
  params:
    samplesheet = config['samplesheet'],
    pseudocount = config['pseudocount'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison]
    outputfolder = lambda wildcards: f"twogroup/{wildcards.comparison}"
  threads: 1
  conda: "envs/seqtools.yml"
  script:
    "scripts/edger_twogroup.R"

# rule maplot:
#   input:
#     '{comparison}/{comparison}_diffacc_edgeR.tsv'
#   output:
#     png = '{comparison}/{comparison}_maplot.png'
#   params:
#     group1 = lambda wildcards: list(
#         config['comparison'][wildcards.comparison].keys()
#     )[0],
#     group2 = lambda wildcards: list(
#         config['comparison'][wildcards.comparison].keys()
#     )[1]
#   run:
#     maplot(
#         input[0],
#         output[0],
#         params.group1,
#         params.group2
#     )

# rule bedfiles:
#   input:
#     '{comparison}/{comparison}_diffacc_edgeR.tsv'
#   output:
#     '{comparison}/diffpeaks_{gr}.bed'
#   params:
#     grint = lambda wildcards: retgroupnumber(
#       wildcards.gr,
#       config['comparison'][wildcards.comparison]
#     )
#   threads: 1
#   run:
#     tsv_to_bed(
#       input[0],
#       output[0],
#       params.grint
#     )

# rule heatmaps:
#   input:
#     maplot = '{comparison}/{comparison}_maplot.png',
#     samples = expand('bw/{sample}.scalefac.bw', sample=SAMPLES),
#   output:
#     matrix = temp('{comparison}/mat.npz'),
#     heatmap = '{comparison}/diffpeaks.png'
#   params:
#     beds = lambda wildcards: wildcards.comparison + '/*bed',
#     samples = lambda wildcards: readsamples(wildcards.comparison)
#   threads: 10
#   conda: "envs/deeptools.yml"
#   shell:'''
#   computeMatrix reference-point -R {params.beds} -S {params.samples} -o {output.matrix} \
#   --referencePoint center \
#   -b 3000 -a 3000 \
#   --missingDataAsZero \
#   --smartLabels \
#   -p {threads}
#   plotHeatmap --matrixFile {output.matrix} --colorMap 'Blues' -o {output.heatmap}
#   '''
