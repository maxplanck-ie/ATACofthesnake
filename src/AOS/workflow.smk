include: config['rules']['peaks']
include: config['rules']['qc']
include: config['rules']['de']
include: config['rules']['motifs']
include: config['rules']['tobias']

def compzip(comp):
  comps = []
  grs = []
  for c in comp:
    for gr in comp[c]:
      comps.append(c)
      grs.append(gr)
  return(comps, grs)

def geto():
  # Peaks
  _f = []
  _f.append(
    'peakset/peaks_uropa_finalhits.txt'
  )
  _f.extend(
    expand(
      'bw/{sample}.scalefac.bw',
      sample=config['samples']
    )
  ),
  _f.extend(
    [
      'figures/mitofraction.png',
      'qc/fragsize.tsv',
      'figures/alignmentsieve.png',
      'figures/fragmentsizes.png'
    ]
  )
  # Differential
  if config['comparison']:
    # define list of comparisons & list per groups to easily zip later on.
    comps, grs = compzip(config['comparison'])
    _f.extend(
      expand(
        '{comparison}/{comparison}_maplot.png',
        comparison=config['comparison'].keys()
      )
    ),
    _f.extend(
      expand(
        '{comparison}/diffpeaks_{gr}.bed',
        zip,
        comparison=comps,
        gr=grs
      )
    )
    _f.extend(
      expand(
        '{comparison}/diffpeaks.png',
        comparison=config['comparison'].keys()
      )
    )
    # Motifs.
    if config['files']['motif']:
      _f.append(
        'motifs_clustered/clusteredmotifs_consensus_motifs.meme'
      )
      _f.extend(
        expand(
          '{comparison}/motif_{gr}/ame.html',
          zip,
          comparison=comps,
          gr=grs
        )
      )
  return (_f)

localrules: lnBams
rule all:
  input:
    geto()
  

