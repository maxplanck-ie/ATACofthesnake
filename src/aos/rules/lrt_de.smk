import os
import glob
from sklearn.cluster import KMeans
import pandas as pd
from aos.helper import get_elbow

def get_beds(wildcards):
    ckpt = checkpoints.lrt_bedfiles.get(comparison=wildcards.comparison)
    beddir = ckpt.output.beddir
    return sorted(glob.glob(os.path.join(beddir, "*.bed")))

rule lrt_diffacc:
  input:
    mat = 'peakset/counts.tsv'
  output:
    table = 'lrt/{comparison}/{comparison}_lrt_edgeR.tsv',
    samples = 'lrt/{comparison}/relevant_samples.txt'
  params:
    samplesheet = config['samplesheet'],
    comparison = lambda wildcards: config['comparison'][wildcards.comparison],
    outputfolder = lambda wildcards: f"lrt/{wildcards.comparison}",
    comparison_name = lambda wildcards: wildcards.comparison
  threads: 1
  conda: "envs/seqtools.yml"
  script:
    "scripts/edger_lrt.R"

checkpoint lrt_bedfiles:
  input:
    edger = 'lrt/{comparison}/{comparison}_lrt_edgeR.tsv',
    matrix = 'peakset/counts.tsv'
  output:
    beddir = directory("lrt/{comparison}/bed")
  run:
    import pandas as pd
    os.makedirs(output.beddir, exist_ok=True)

    df = pd.read_csv(input.edger, sep='\t', header=0)
    sig = df[df['FDR'] < 0.05]
    if len(sig) > 100:
        # Write out the bed file
        with open(f"{output.beddir}/{wildcards.comparison}_LRT_sign.bed", "w") as f:
            for peak in sig["peak_id"]:
                f.write("\t".join(peak.split("|")) + "\n")
        # Next, run determine k for heatmap
        counts = pd.read_csv(input.matrix, sep='\t', header=0)
        peaks = counts[["#chr", "start", "end"]].astype(str).agg("|".join, axis=1)
        counts = counts.drop(columns=["#chr", "start", "end"])
        counts = counts[peaks.isin(sig["peak_id"])].values
        counts_scaled = counts - counts.mean()
        counts_scaled /= counts.std() + 1e-8
        k_range = range(2, 20)
        inertias = []
        for k in k_range:
            kmeans = KMeans(n_clusters=k, random_state=1337).fit(counts_scaled)
            inertias.append(kmeans.inertia_)
        k_opt = get_elbow(inertias, k_range)
        print(f"Selected K = {k_opt}")
        with open(f"lrt/{wildcards.comparison}/k_opt.txt", 'w') as f:
            f.write(str(k_opt))


rule lrt_plotheatmap:
  input:
    samples = expand('bw/{sample}.scalefac.bw', sample=SAMPLES),
    bedfiles = get_beds,
    relevant_samples = "lrt/{comparison}/relevant_samples.txt"
  output:
    touch("lrt/{comparison}/{comparison}_heatmap.done"),
  params:
    matrix = lambda wildcards: f"lrt/{wildcards.comparison}/{wildcards.comparison}_mat.npz",
    heatmap = lambda wildcards: f"lrt/{wildcards.comparison}/{wildcards.comparison}_lrtpeaks.png",
  conda: "envs/deeptools.yml"
  threads: 20
  run:
    if input.bedfiles:
      bedfiles_str = ' '.join(str(b) for b in input.bedfiles)
      relevant_samples = []
      with open(input.relevant_samples, 'r') as f:
        relevant_samples = [f"bw/{line.strip()}.scalefac.bw" for line in f]
      relevant_samples_str = ' '.join(relevant_samples)
      with open(f'lrt/{wildcards.comparison}/k_opt.txt', 'r') as f:
        k_opt = int(f.read().strip())
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
        plotHeatmap --matrixFile {params.matrix} --colorMap 'Blues' --kmeans {k_opt} -o {params.heatmap}
      """)
    else:
      Path(output[0]).touch()