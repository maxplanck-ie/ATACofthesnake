include: config['rules']['peaks']
include: config['rules']['qc']
include: config['rules']['de']

def peaks():
  return (
    [
      'peakset/peaks_uropa_finalhits.txt',
      expand(
        'bw/{sample}.scalefac.bw',
        sample=config['samples']
      ),
      expand(
        'bw/{sample}.RPKM.bw',
        sample=config['samples']
      )
    ]
  )

def qc():
  return (
    [
      'figures/mitofraction.png',
      'qc/fragsize.tsv'
    ]
  )

def de():
  if config['comparison']:
    return (
      [
        'ret'
      ]
    )
  else:
    return ('')

localrules: lnBams
rule all:
  input:
    peaks(),
    # QC
    qc(),
    # DE
    #de()

