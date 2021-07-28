import os
from ATACofthesnake import misc

ss = misc.readss(config['sampleSheet'], config['Samples'])

rule all:
    input:
        expand(config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up.fna', Comp=config['diffComp']),
        expand(config['outDir'] + '/Figures/{Comp}_Diffpeak.png', Comp=config['diffComp']),
        expand(config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up/ame.html', Comp=config['diffComp']),
        expand(config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up_shuffle/ame.html', Comp=config['diffComp'])

rule formatFNA:
    input:
        up = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv",
        down = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_DOWN.tsv",
        fai = config['genomeFa'] + '.fai'
    output:
        up = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up.fna',
        down = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_down.fna',
        upBed = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up.bed',
        downBed = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_down.bed'
    threads: 1
    params:
        genomeFa = config['genomeFa']
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    log:
        out = config['outDir'] + '/logs/formatFNA_{Comp}.out',
        err = config['outDir'] + '/logs/formatFNA_{Comp}.err'
    shell:'''
    set +o pipefail;cut -f1 {input.up} | tail -n +2 | tr '_' '\t' > {output.upBed}
    set +o pipefail;cut -f1 {input.down} | tail -n +2 | tr '_' '\t' > {output.downBed}
    bedtools getfasta -fi {params.genomeFa} -bed {output.upBed} > {output.up}
    bedtools getfasta -fi {params.genomeFa} -bed {output.downBed} > {output.down}
    '''

rule diffHeat:
    input:
        upBed = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up.bed',
        downBed = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_down.bed',
        bigWigs = lambda wildcards: expand(config['outDir'] + '/BigWigs/{sample}.bw', sample=ss[wildcards.Comp])
    output:
        heatmap = config['outDir'] + "/Figures/{Comp}_Diffpeak.png",
        outMat = config['outDir'] + '/deepTools/{Comp}_diff.BigwigMatrix.gz'
    threads: 10
    params: 
        bigwigs = lambda wildcards: ' '.join(expand(config['outDir'] + '/BigWigs/{sample}.bw', sample=ss[wildcards.Comp])),
        Comp = lambda wildcards: wildcards.Comp
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    log:
        Plotout = config['outDir'] + "/logs/diffplot_{Comp}.out",
        Ploterr = config['outDir'] + "/logs/diffplot_{Comp}.err",
        Matout = config['outDir'] + "/logs/diffMat_{Comp}.out",
        Materr = config['outDir'] + "/logs/diffMat_{Comp}.err"
    shell:'''
    computeMatrix reference-point -S {params.bigwigs} -R {input.upBed} {input.downBed} --referencePoint center -a 5000 -b 5000 -out {output.outMat} -p {threads} --missingDataAsZero > {log.Matout} 2> {log.Materr}
    plotHeatmap -m {output.outMat} -out {output.heatmap} --refPointLabel center --colorMap "Blues" > {log.Plotout} 2> {log.Ploterr}
    '''

rule runMeme:
    input:
        up = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up.fna',
        down = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_down.fna'
    output:
        motifUP = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up/ame.html',
        motifDO = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_down/ame.html'
    params:
        motif = config['motifLoc'],
        upOut = lambda wildcards: config['outDir'] + '/diffAcc_' + wildcards.Comp + '/Motif/' + wildcards.Comp + "_up",
        downOut = lambda wildcards: config['outDir'] + '/diffAcc_' + wildcards.Comp + '/Motif/' + wildcards.Comp + "_down"
    log:
        upOut = config['outDir'] + "/logs/runMeme_up_{Comp}.out",
        upErr = config['outDir'] + "/logs/runMeme_up_{Comp}.err",
        downOut = config['outDir'] + "/logs/runMeme_up_{Comp}.out",
        downErr = config['outDir'] + "/logs/runMeme_up_{Comp}.err"
    threads: 1
    conda: os.path.join(config['baseDir'], 'envs','AOS_meme.yaml')
    shell:'''
    ame --control {input.down} -o {params.upOut} {input.up} {params.motif} > {log.upOut} 2> {log.upErr}
    ame --control {input.up} -o {params.downOut} {input.down} {params.motif} > {log.downOut} 2> {log.downErr}
    '''

rule runMeme_shuffled:
    input:
        up = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up.fna',
        down = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_down.fna'
    output:
        motifUP = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_up_shuffle/ame.html',
        motifDO = config['outDir'] + '/diffAcc_{Comp}/Motif/{Comp}_down_shuffle/ame.html'
    params:
        motif = config['motifLoc'],
        upOut = lambda wildcards: config['outDir'] + '/diffAcc_' + wildcards.Comp + '/Motif/' + wildcards.Comp + "_up_shuffle",
        downOut = lambda wildcards: config['outDir'] + '/diffAcc_' + wildcards.Comp + '/Motif/' + wildcards.Comp + "_down_shuffle"
    log:
        upOut = config['outDir'] + "/logs/runMeme_up_{Comp}.out",
        upErr = config['outDir'] + "/logs/runMeme_up_{Comp}.err",
        downOut = config['outDir'] + "/logs/runMeme_up_{Comp}.out",
        downErr = config['outDir'] + "/logs/runMeme_up_{Comp}.err"
    threads: 1
    conda: os.path.join(config['baseDir'], 'envs','AOS_meme.yaml')
    shell:'''
    ame --control --shuffle-- -o {params.upOut} {input.up} {params.motif} > {log.upOut} 2> {log.upErr}
    ame --control --shuffle-- -o {params.downOut} {input.down} {params.motif} > {log.downOut} 2> {log.downErr}
    '''