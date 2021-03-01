import os
import yaml
from ATACofthesnake import misc

# Read / set variables.
with open('Parameters_downstream.yaml') as f:
    paramDic = yaml.load(f, Loader=yaml.FullLoader)

rule all:
    input:
        expand('Motif/{Comp}_up.fna', Comp=paramDic['diffComp']),
        expand('Figures/{Comp}_Diffpeak.png', Comp=paramDic['diffComp']),
        expand('Motif/{Comp}_up/ame.html', Comp=paramDic['diffComp'])

# Does it make sense to run this as a seperate workflow?

rule formatFNA:
    input:
        up = "diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv",
        down = "diffAcc_{Comp}/{Comp}_edgeR_annotated_DOWN.tsv",
        fai = paramDic['genomeFa'] + '.fai',
    output:
        up = 'Motif/{Comp}_up.fna',
        down = 'Motif/{Comp}_down.fna',
        upBed = 'Motif/{Comp}_up.bed',
        downBed = 'Motif/{Comp}_down.bed'
    threads: 1
    params:
        genomeFa = paramDic['genomeFa']
    conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
    log:
        out = 'logs/formatFNA_{Comp}.out',
        err = 'logs/formatFNA_{Comp}.err'
    shell:'''
    set +o pipefail;cut -f1 {input.up} | tail -n +2 | tr '_' '\t' > {output.upBed}
    set +o pipefail;cut -f1 {input.down} | tail -n +2 | tr '_' '\t' > {output.downBed}
    bedtools getfasta -fi {params.genomeFa} -bed {output.upBed} > {output.up}
    bedtools getfasta -fi {params.genomeFa} -bed {output.downBed} > {output.down}
    '''

rule diffHeat:
    input:
        upBed = 'Motif/{Comp}_up.bed',
        downBed = 'Motif/{Comp}_down.bed',
        bigWigs = expand('BigWigs/{sample}.bw', sample=paramDic['Samples'])
    output:
        heatmap = "Figures/{Comp}_Diffpeak.png",
        outMat = 'deepTools/{Comp}_diff.BigwigMatrix.gz'
    threads: 10
    params: 
        bigwigs = lambda wildcards: ' '.join(expand('BigWigs/{sample}.bw', sample=paramDic['Comp'][wildcards.Comp]['Samples'])),
        Comp = lambda wildcards: wildcards.Comp
    conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
    log:
        Plotout = "logs/diffplot_{Comp}.out",
        Ploterr = "logs/diffplot_{Comp}.err",
        Matout = "logs/diffMat_{Comp}.out",
        Materr = "logs/diffMat_{Comp}.err"
    shell:'''
    computeMatrix reference-point -S {params.bigwigs} -R {input.upBed} {input.downBed} --referencePoint center -a 2000 -b 2000 -out {output.outMat} -p {threads} --missingDataAsZero > {log.Matout} 2> {log.Materr}
    plotHeatmap -m {output.outMat} -out {output.heatmap} --refPointLabel center --colorMap Blues > {log.Plotout} 2> {log.Ploterr}
    '''

rule runMeme:
    input:
        up = 'Motif/{Comp}_up.fna',
        down = 'Motif/{Comp}_down.fna'
    output:
        motifUP = 'Motif/{Comp}_up/ame.html',
        motifDO = 'Motif/{Comp}_down/ame.html'
    params:
        motif = paramDic['motifLoc'],
        upOut = lambda wildcards: wildcards.Comp + "_up",
        downOut = lambda wildcards: wildcards.Comp + "_down"
    threads: 1
    conda: os.path.join(paramDic['baseDir'], 'envs','AOS_meme.yaml')
    shell:'''
    ame --control {input.down} -o Motif/{params.upOut} {input.up} {params.motif}
    ame --control {input.up} -o Motif/{params.downOut} {input.down} {params.motif}
    '''

# Todo.
rule TOBIAS:
    input:
        "MACS2/{Comp}_Merged_peaks.narrowPeak"
    output:
        "TOBIAS_{Comp}/footprints.txt"
    threads: 1
    shell:'''
    touch {output}
    '''