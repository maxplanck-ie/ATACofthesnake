include: config['rules']['peaks']
include: config['rules']['qc']
include: config['rules']['de']

def geto():
  _f = []
  _f.append(
    'peakset/peaks_uropa_finalhits.txt'
  )
  _f.extend(
    expand(
      'bw/{sample}.scalefac.bw',
      sample=config['samples']
    )
  )
  _f.extend(
    [
      'figures/mitofraction.png',
      'qc/fragsize.tsv',
      'figures/alignmentsieve.png',
      'figures/fragmentsizes.png'
    ]
  )
  if config['comparison']:
    _f.extend(
      expand(
        '{comparison}/{comparison}_maplot.png',
        comparison=config['comparison'].keys()
      )
    )
    _f.extend(
      expand(
        '{comparison}/diffpeaks.png',
        comparison=config['comparison'].keys()
      )
    )
  print(_f)
  return (_f)

localrules: lnBams
rule all:
  input:
    geto()
  

