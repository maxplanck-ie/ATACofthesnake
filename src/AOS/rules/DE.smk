# Create a list containing:
from AOS.helper import maplot
import os

def readsamples(_f):
  samplef = os.path.join(
    _f,'samples.txt'
  )
  if not os.path.exists(samplef):
    return('')
  samples = []
  with open(samplef) as f:
    for line in f:
      samples.append(
        'bw/{}.scalefac.bw'.format(
          line.strip()
        )
      )
  return (samples)

def prefix(comparisondic):
  return (
    '_vs_'.join(
      list(comparisondic.keys())
    )
  )

def retsubset(comparisondic, group):
  gr = list(comparisondic.keys())[group]
  tlis = []
  for fact in comparisondic[gr]:
    tlis.append(r"samplesheet${}=='{}'".format(
        fact,
        comparisondic[gr][fact]
    ))
  return ('&'.join(
    tlis
  ))
  return (
    ','.join(
        list(comparisondic[gr].keys())
    )
  )

rule diffacc:
  input:
    'peakset/counts.tsv'
  output:
    table = '{comparison}/{comparison}_diffacc_edgeR.tsv',
    samples = '{comparison}/samples.txt'
  params:
    samplesheet = config['files']['samplesheet'],
    prefix = lambda wildcards: prefix(
        config['comparison'][wildcards.comparison],
    ),
    gr1_subset = lambda wildcards: retsubset(
        config['comparison'][wildcards.comparison],
        0
    ),
    gr2_subset = lambda wildcards: retsubset(
        config['comparison'][wildcards.comparison],
        1
    ),
    interaction = config['interaction'],
    outputfolder = lambda wildcards: wildcards.comparison
  threads: 1
  conda: config['envs']['seqtools']
  script:
    config['rscripts']['edger']

rule maplot:
  input:
    '{comparison}/{comparison}_diffacc_edgeR.tsv'
  output:
    png = '{comparison}/{comparison}_maplot.png'
  params:
    group1 = lambda wildcards: list(
        config['comparison'][wildcards.comparison].keys()
    )[0],
    group2 = lambda wildcards: list(
        config['comparison'][wildcards.comparison].keys()
    )[1]
  run:
    maplot(
        input[0],
        output[0],
        params.group1,
        params.group2
    )


rule heatmaps:
  input:
    maplot = '{comparison}/{comparison}_maplot.png'
  output:
    matrix = temp('{comparison}/mat.npz'),
    heatmap = '{comparison}/diffpeaks.png'
  params:
    beds = lambda wildcards: wildcards.comparison + '/*bed',
    samples = lambda wildcards: readsamples(wildcards.comparison)
  threads: 10
  conda: config['envs']['seqtools']
  shell:'''
  computeMatrix reference-point -R {params.beds} -S {params.samples} -o {output.matrix} \
  --referencePoint center \
  -b 3000 -a 3000 \
  --missingDataAsZero \
  --smartLabels \
  -p {threads}
  plotHeatmap --matrixFile {output.matrix} --colorMap 'Blues' -o {output.heatmap}
  '''
