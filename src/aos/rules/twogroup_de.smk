# Create a list containing:
from aos.helper import maplot
import os
from pathlib import Path
import glob

def get_beds(wildcards):
    ckpt = checkpoints.twogroup_bedfiles.get(comparison=wildcards.comparison)
    beddir = ckpt.output.beddir
    return sorted(glob.glob(os.path.join(beddir, "*.bed")))


rule twogroup_diffacc:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'twogroup/{comparison}/{comparison}_diffacc_edgeR.tsv',
    samples = 'twogroup/{comparison}/relevant_samples.txt'
  params:
    samplesheet = config['samplesheet'],
    pseudocount = config['pseudocount'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison],
    outputfolder = lambda wildcards: f"twogroup/{wildcards.comparison}",
    comparison_name = lambda wildcards: wildcards.comparison
  threads: 1
  conda: "envs/seqtools.yml"
  script:
    "scripts/edger_twogroup.R"

rule twogroup_maplot:
  input:
    'twogroup/{comparison}/{comparison}_diffacc_edgeR.tsv'
  output:
    png = 'twogroup/{comparison}/{comparison}_maplot.png'
  params:
    comparison = lambda wildcards: config['comparison'][wildcards.comparison].copy(),
  run:
    maplot(input[0], output[0], params.comparison)

checkpoint twogroup_bedfiles:
  input:
    'twogroup/{comparison}/{comparison}_diffacc_edgeR.tsv'
  output:
    beddir = directory("twogroup/{comparison}/bed")
  run:
    import pandas as pd
    os.makedirs(output.beddir, exist_ok=True)

    df = pd.read_csv(input[0], sep='\t', header=0)
    gr1 = df[(df['logFC'] < 0) & (df['FDR'] < 0.05)]
    gr2 = df[(df['logFC'] > 0) & (df['FDR'] < 0.05)]
    for gr in [gr1, gr2]:
      if len(gr) > 0:
        gr_name = gr['group_assignment'].iloc[0]
        with open(f"{output.beddir}/{gr_name}.bed", "w") as f:
          for peak in gr["peak_id"]:
            f.write("\t".join(peak.split("|")) + "\n")


rule twogroup_plotheatmap:
  input:
    samples = expand('bw/{sample}.scalefac.bw', sample=SAMPLES),
    bedfiles = get_beds,
    relevant_samples = "twogroup/{comparison}/relevant_samples.txt"
  output:
    touch("twogroup/{comparison}/{comparison}_heatmap.done"),
  params:
    matrix = lambda wildcards: f"twogroup/{wildcards.comparison}/{wildcards.comparison}_mat.npz",
    heatmap = lambda wildcards: f"twogroup/{wildcards.comparison}/{wildcards.comparison}_diffpeaks.png",
  conda: "envs/deeptools.yml"
  threads: 20
  run:
    if input.bedfiles:
      bedfiles_str = ' '.join(str(b) for b in input.bedfiles)
      relevant_samples = []
      with open(input.relevant_samples, 'r') as f:
        relevant_samples = [f"bw/{line.strip()}.scalefac.bw" for line in f]
      relevant_samples_str = ' '.join(relevant_samples)
      shell(f'''
        computeMatrix reference-point \
          -R {bedfiles_str} \
          -S {relevant_samples_str} \
          -o {params.matrix} \
          --referencePoint center \
          -b 3000 -a 3000 \
          --missingDataAsZero \
          --smartLabels \
          -p {threads}
      ''')
      shell(f"""
        plotHeatmap --matrixFile {params.matrix} --colorMap 'Blues' -o {params.heatmap}
      """)
    else:
      Path(output[0]).touch()
