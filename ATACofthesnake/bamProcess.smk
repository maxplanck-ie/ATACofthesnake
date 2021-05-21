import os
import shutil
from ATACofthesnake import misc

# A couple arg flags define some restructuring:
# peakSet: if an external peakset is provided, only execute fragSize, alignmentSieve and idxStats (mitobleed).
# mergeBam: merge bam files before peak calling. If a sampleSheet is provided, bam files are merged per condition (within a comparison).

if not config['sampleSheet']:
    ss = {}
    # Set the comparison to the output folder name (we don't have one).
    compStr = str(config['outDir'].split('/')[-1])
    ss['Comp'] = [compStr]
    ss[compStr] = config['Samples']
else:
    ss = misc.readss(config['sampleSheet'], config['Samples'])

print(ss)

outList = []
outList.append(
    expand(config['outDir'] + "/Figures/{CompCond}.mtFrac.png", CompCond = ss['Comp']) +
    expand(config['outDir'] + "/ShortBAM/{Sample}.bam.bai", Sample=config['Samples']) +
    expand(config['outDir'] + "/deepTools/{CompCond}.raw.fragSize.tsv", CompCond = ss['Comp'])
)

if config['mergeBam']:
    if not config['sampleSheet']:
        # We need to merge bam files before peak calling
        # There is no sampleSheet --> merge all bam files.
        outList.append(
            expand(config['outDir'] + "/mergeBAM/{CompCond}.bam", CompCond=ss['Comp']) +
            expand(config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.bed", CompCond=ss['Comp'])
        )
    elif config['sampleSheet']:
        # Merge bamFiles per condition
        outList.append(
            expand(config['outDir'] + "/mergeBAM/{CompCond}.bam", CompCond=ss['CompCond']) +
            #expand(config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.bed", CompCond=ss['CompCond']) +
            expand(config['outDir'] + '/MACS2_mergeBAM/{Comp}_peaks.bed', Comp=ss['Comp'])
        )
else:
    outList.append(
        expand(config['outDir'] + "/MACS2/{CompCond}_peaks.bed", CompCond=ss['Comp'])
    )
#Common rules:
outList.append(
    # centralize peaks, FRIP plot, Uropa
    expand(config['outDir'] + "/Peaks/{Comp}_peaks.bed", Comp=ss['Comp']) +
    expand(config['outDir'] + "/Figures/{Comp}.FRIP.png", Comp=ss['Comp']) +
    expand(config['outDir'] + "/Annotation/{Comp}_uropa_finalhits.txt", Comp=ss['Comp']) +
    expand(config['outDir'] + "/Figures/{Comp}_Heatmap.png", Comp=ss['Comp']) +
    expand(config['outDir'] + "/Figures/{Comp}_PCA.png", Comp=ss['Comp']) +
    expand(config['outDir'] + "/Figures/{Comp}_plotCorr_pearson.png", Comp=ss['Comp'])
)
# differential accessibility:
if config['sampleSheet']:
    outList.append(
        expand(config['outDir'] + "/Figures/{Comp}_maPlot.png", Comp=ss['Comp']) +
        expand(config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv", Comp=ss['Comp'])
    )
print(outList)

ruleorder: MACS2_nptobed > unionMACS2_merge 
localrules: fripPlotter, idxStatPlotter, maPlot
rule all:
    input:
        outList


rule checkGenomeIndex:
	input: config['genomeFa']
	output: config['genomeFa'] + '.fai'
	log:
		out = config['outDir'] + "/logs/checkFai.out",
		err = config['outDir'] + "/logs/checkFai.err"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs', 'AOS_SeqTools.yaml')
	shell:'''
	samtools faidx {input}
	'''
rule checkIndex:
	input:
		sample = config['bamDir'] + '/{sample}.bam',
		genomeFai = config['genomeFa'] + '.fai'
	output:
		config['bamDir'] + '/{sample}.bam.bai'
	log:
		out = config['outDir'] + '/logs/checkIndex.{sample}.out',
		err = config['outDir'] + '/logs/checkIndex.{sample}.err'
	threads: 10
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	sambamba index -t {threads} {input.sample} > {log.out} 2> {log.err}
	'''

rule idxStat:
	input:
		index = config['bamDir'] + '/{sample}.bam.bai',
		bam = config['bamDir'] + '/{sample}.bam'
	output:
		config['outDir'] + "/QC/{sample}.idxstat.txt"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	set +o pipefail;
	samtools idxstats {input.bam} | cut -f1,3 > {output}
    '''

rule idxStatPlotter:
	input: 
		lambda wildcards: expand(config['outDir'] + "/QC/{sample}.idxstat.txt", sample=ss[wildcards.CompCond])
	output:
		config['outDir'] + "/Figures/{CompCond}.mtFrac.png"
	threads: 1
	run:
		misc.plotter('idxstat', input, str(output))

rule alignmentSieve:
	input:
		inBam = config['bamDir'] + '/{sample}.bam',
		index = config['bamDir'] + '/{sample}.bam.bai',
		mitoBleed = config['outDir'] + "/QC/{sample}.idxstat.txt"
	output:
		shortBam = config['outDir'] + "/ShortBAM/{sample}.bam",
		filterMetrics = config['outDir'] + "/ShortBAM/{sample}.metrics"
	params:
		blackList = config['blackList'],
		fragSize = config['fragSize']
	log:
		out = config['outDir'] + '/logs/alignmentSieve.{sample}.out',
		err = config['outDir'] + '/logs/alignmentSieve.{sample}.err'
	threads:10
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	alignmentSieve --bam {input.inBam} --outFile {output.shortBam} -p {threads} --filterMetrics {output.filterMetrics} --maxFragmentLength {params.fragSize} --minFragmentLength 0 --blackListFileName {params.blackList} > {log.out} 2> {log.err}
	'''

rule shortIndex:
	input:
		config['outDir'] + "/ShortBAM/{sample}.bam"
	output:
		index = config['outDir'] + "/ShortBAM/{sample}.bam.bai"
	log:
		out = config['outDir'] + '/logs/shortIndex.{sample}.out',
		err = config['outDir'] + '/logs/shortIndex.{sample}.err'
	threads: 10
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	sambamba index {input} > {log.out} 2> {log.err}
	'''

rule fragSize:
    input:
        lambda wildcards: expand(config['outDir'] + "/ShortBAM/{sample}.bam.bai", sample=ss[wildcards.CompCond])
    output:
        raw = config['outDir'] + '/deepTools/{CompCond}.raw.fragSize.tsv',
        table = config['outDir'] + '/deepTools/{CompCond}.fragSize.tsv'
    params:
        lambda wildcards: ' '.join(expand(config['outDir'] + "/ShortBAM/{sample}.bam", sample=ss[wildcards.CompCond]))
    log:
        out = config['outDir'] + '/logs/shortIndex.{CompCond}.out',
        err = config['outDir'] + '/logs/shortIndex.{CompCond}.err',
    threads: 5
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    shell:'''
    bamPEFragmentSize -b {params} -p {threads} --outRawFragmentLengths {output.raw} --table {output.table} > {log.out} 2> {log.err}
    '''

# Rules for Merging.
rule mergeBam:
    input:
        lambda wildcards: expand(config['outDir'] + "/ShortBAM/{sample}.bam.bai", sample=ss[wildcards.CompCond])
    output:
        config['outDir'] + "/mergeBAM/{CompCond}.bam"
    params:
        lambda wildcards: ' '.join(expand(config['outDir'] + "/ShortBAM/{sample}.bam", sample=ss[wildcards.CompCond]))
    log:
        out = config['outDir'] + '/logs/mergeBam.{CompCond}.out',
        err = config['outDir'] + '/logs/mergeBam.{CompCond}.err'
    threads: 5
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    shell:'''
    sambamba merge -t 5 {output} {params}
    '''

rule mergeBamtoBed:
    input:
        config['outDir'] + "/mergeBAM/{CompCond}.bam"
    output:
        config['outDir'] + "/mergeBAM/{CompCond}.bed"
    threads: 1
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    shell:'''
    bamToBed -i {input} > {output}
    '''
rule MACS2_mergeBam:
	input:
		config['outDir'] + "/mergeBAM/{CompCond}.bed"
	output:
		config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.narrowPeak"
	log:
		out = config['outDir'] + '/logs/MACS2_mergeBAM.{CompCond}.out',
		err = config['outDir'] + '/logs/MACS2_mergeBAM.{CompCond}.err'
	params:
		genomeSize = config['genomeSize'],
		outName = lambda wildcards: wildcards.CompCond,
		blackList = config['blackList'],
		outDir = config['outDir'] + "/MACS2_mergeBAM"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	macs2 callpeak -t {input} -f BED --nomodel --shift -75 --extsize 150 -g {params.genomeSize} -n {params.outName} -q 0.01 --outdir {params.outDir} --keep-dup all > {log.out} 2> {log.err}
	'''
rule MACS2_nptobed:
    input:
        config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.narrowPeak"
    output:
        config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.bed"
    log:
        out = config['outDir'] + '/logs/MACS2_nptobed.{CompCond}.out',
        err = config['outDir'] + '/logs/MACS2_nptobed.{CompCond}.err'
    threads: 1
    shell:'''
    cut -f1-3 {input} > {output}
    '''


rule unionMACS2_merge:
    input:
        lambda wildcards: expand(config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.narrowPeak", CompCond=ss['CompCondDic'][wildcards.Comp])
    output:
        config['outDir'] + '/MACS2_mergeBAM/{Comp}_peaks.bed'
    params:
        lambda wildcards: ' '.join(expand(config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.narrowPeak", CompCond=ss['CompCondDic'][wildcards.Comp]))
    threads: 1
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    shell:'''
    cat {params} | sort -k1,1 -k2,2n | bedtools merge > {output}
    '''

# Rules for peak per sample
rule bamtobed:
    input:
        index = config['outDir'] + "/ShortBAM/{sample}.bam.bai",
        bamFile = config['outDir'] + "/ShortBAM/{sample}.bam"
    output:
        config['outDir'] + "/ShortBAM/{sample}.bed"
    threads: 2
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    shell:'''
    bamToBed -i {input.bamFile} > {output}
    '''
rule MACS2:
	input:
		config['outDir'] + "/ShortBAM/{sample}.bed"
	output:
		config['outDir'] + "/MACS2/{sample}_peaks.narrowPeak"
	log:
		out = config['outDir'] + '/logs/MACS2.{sample}.out',
		err = config['outDir'] + '/logs/MACS2.{sample}.err'
	params:
		genomeSize = config['genomeSize'],
		outName = lambda wildcards: wildcards.sample,
		blackList = config['blackList'],
		outDir = config['outDir'] + "/MACS2"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	macs2 callpeak -t {input} -f BED --nomodel --shift -75 --extsize 150 -g {params.genomeSize} -n {params.outName} -q 0.01 --outdir {params.outDir} --keep-dup all > {log.out} 2> {log.err}
	'''

rule unionMACS2:
    input:
        lambda wildcards: expand(config['outDir'] + "/MACS2/{sample}_peaks.narrowPeak", sample=ss[wildcards.CompCond])
    output:
        config['outDir'] + '/MACS2/{CompCond}_peaks.bed'
    params:
        lambda wildcards: ' '.join(expand(config['outDir'] + "/MACS2/{sample}_peaks.narrowPeak", sample=ss[wildcards.CompCond]))
    threads: 1
    conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
    shell:'''
    cat {params} | sort -k1,1 -k2,2n | bedtools merge > {output}
    '''

# Downstream - centralize peaks.
rule centralPeak:
    input:
        config['peakSet'] if config['peakSet'] else misc.returnPeaks(config['outDir'], "{Comp}", config['mergeBam'])
    output:
        config['outDir'] +"/Peaks/{Comp}_peaks.bed"
    threads: 1
    shell:'''
    cp {input} {output}
    '''

# postPeaks - Common rules.
rule fripScore:
	input:
		bamfile = config['outDir'] + "/ShortBAM/{sample}.bam",
		Peaks = lambda wildcards: misc.returnCompfromSample(wildcards.sample,ss, config['outDir'])
	output:
		config['outDir'] + "/QC/{sample}.FRiP.txt"
	params:
		genomeSize = config['genomeSize'],
		peakFile = lambda wildcards: misc.returnCompfromSample(wildcards.sample ,ss, config['outDir']),
		sample = lambda wildcards: wildcards.sample
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	peakcount=$(cat {params.peakFile} | wc -l)
	mapped=$(samtools view -c -F 4 {input.bamfile})
	peakReads=$(samtools view -c -F 4 -L {params.peakFile} {input.bamfile})
	frip=$(bc -l <<< "$peakReads/$mapped")
	peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' {params.peakFile})
	genomecov=$(bc -l <<< "$peak_len/{params.genomeSize}")
	printf "sample\tpeakcount\tfrip\tpeak_genome_coverage\n%s\t%d\t%5.3f\t%6.4f\n" {params.sample} $peakcount $frip $genomecov > {output}
	'''

rule fripPlotter:
	input: 
		lambda wildcards: expand(config['outDir'] + "/QC/{sample}.FRiP.txt", sample=ss[wildcards.Comp])
	output:
		config['outDir'] + "/Figures/{Comp}.FRIP.png"
	threads: 1
	params: config['outDir']
	run:
		misc.plotter('frip', input, str(output), outDir=str(params))

rule uropa:
	input:
		config['outDir'] +"/Peaks/{Comp}_peaks.bed"
	output:
		config['outDir'] + "/Annotation/{Comp}_uropa_finalhits.txt"
	log:
		out = config['outDir'] + "/logs/uropa.{Comp}.out",
		err = config['outDir'] + "/logs/uropa.{Comp}.err"
	params:
		GTF = config['GTF'],
		prefix = "{Comp}_uropa",
		outDir = config['outDir'] + "/Annotation"
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	threads: 5
	shell:'''
	uropa -b {input} -g {params.GTF} --summary --feature transcript --distance 20000 10000 --internals 1 -p {params.prefix} -o {params.outDir} -t {threads} --show-attributes gene_id transcript_id gene_name gene_type transcript_type > {log.out} 2> {log.err}
	'''

rule countMat:
	input:
		config['outDir'] +"/Peaks/{Comp}_peaks.bed"
	output:
		mat = config['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat",
		matnpz = config['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.npz"
	log:
		out = config['outDir'] + "/logs/countMat.{Comp}.out",
		err = config['outDir'] + "/logs/countMat.{Comp}.err"
	params:
		blackList = config['blackList'],
		samples = lambda wildcards: ' '.join(expand(config['bamDir'] + "/{sample}.bam", sample=ss[wildcards.Comp]))
	threads: 20
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	multiBamSummary BED-file --BED {input} -bl {params.blackList} -p {threads} --outRawCounts {output.mat} -o {output.matnpz} -b {params.samples} > {log.out} 2> {log.err}
	# Remove quotes placed by multiBamSummary
	sed -i "s/'//g" {output.mat}
	# Remove .bam postfix
	sed -i 's/\.bam//g' {output.mat}
	'''

rule scaleFactors:
	input:
		config['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat",
	output:
		config['outDir'] + "/diffAcc_{Comp}/{Comp}_scaleFactors.txt"
	log:
		out = config['outDir'] + "/logs/scaleFactors.{Comp}.out",
		err = config['outDir'] + "/logs/scaleFactors.{Comp}.err"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	params:
		scriptLoc = os.path.join(config["baseDir"], "Rscripts", "edgeR_scalefactors.R")
	shell:'''
	Rscript {params.scriptLoc} {input} {output} > {log.out} 2> {log.err}
	'''

rule BigWigs:
	input:
		sf = expand(config['outDir'] + "/diffAcc_{Comp}/{Comp}_scaleFactors.txt",Comp=ss['Comp']),
		inFile = config['outDir'] + "/ShortBAM/{sample}.bam"
	output:
		config['outDir'] + "/BigWigs/{sample}.bw"
	log:
		out = config['outDir'] + "/logs/BigWigs.{sample}.out",
		err = config['outDir'] + "/logs/Bigwigs.{sample}.err"
	params:
		sampleName = '{sample}',
		blackList = config['blackList'],
		genomeSize = config['genomeSize']
	threads: 10
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	SCALEFAC=$(grep {params.sampleName} {input.sf} | cut -f2 -d ' ')
	bamCoverage --scaleFactor $SCALEFAC -b {input.inFile} -o {output} -p {threads} -bs 1 -bl {params.blackList} > {log.out} 2> {log.err}
	'''

rule computeMatrix:
	input:
		lambda wildcards: expand(config['outDir'] + '/BigWigs/{sample}.bw', sample=ss[wildcards.Comp])
	output:
		config['outDir'] + "/deepTools/{Comp}_BigwigMatrix.gz"
	log:
		out = config['outDir'] + "/logs/computeMatrix.{Comp}.out",
		err = config['outDir'] + "/logs/computeMatrix.{Comp}.err"
	params:
		bed = config['outDir'] + "/TSS.bed"
	threads: 10
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	computeMatrix reference-point -S {input} -R {params.bed} --referencePoint center -a 2000 -b 2000 -out {output} -p {threads} --missingDataAsZero > {log.out} 2> {log.err}
	'''

rule plotHeatmap:
	input:
		config['outDir'] + "/deepTools/{Comp}_BigwigMatrix.gz"
	output:
		config['outDir'] + "/Figures/{Comp}_Heatmap.png"
	log:
		out = config['outDir'] + "/logs/plotHeatmap.{Comp}.out",
		err = config['outDir'] + "/logs/plotHeatmap.{Comp}.err"
	threads: 4
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	plotHeatmap -m {input} -out {output} --refPointLabel TSS > {log.out} 2> {log.err}
	'''

rule multiBigwigSum:
	input:
		bw = lambda wildcards: expand(config['outDir'] + '/BigWigs/{sample}.bw', sample=ss[wildcards.Comp]),
		Peaks = config['outDir'] +"/Peaks/{Comp}_peaks.bed"
	output:
		config['outDir'] + '/deepTools/{Comp}_BigwigSum.npz'
	log:
		out = config['outDir'] + "/logs/multiBigWigSum.{Comp}.out",
		err = config['outDir'] + "/logs/multiBigWigSum.{Comp}.err"
	params:
		blackList = config['blackList']
	threads: 10
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	multiBigwigSummary BED-file --BED {input.Peaks} -o {output} -b {input.bw} -bl {params.blackList} -p {threads} -bs 150 > {log.out} 2> {log.err}
	'''

rule plotCorr:
	input:
		config['outDir'] + "/deepTools/{Comp}_BigwigSum.npz"
	output:
		pear = config['outDir'] + "/Figures/{Comp}_plotCorr_pearson.png",
		spear = config['outDir'] + "/Figures/{Comp}_plotCorr_spearman.png"
	log:
		outpear = config['outDir'] + "/logs/plotCorr.{Comp}.out",
		outspear = config['outDir'] + "/logs/plotCorr.{Comp}.out",
		errpear = config['outDir'] + "/logs/plotCorr.{Comp}.err",
		errspear = config['outDir'] + "/logs/plotCorr.{Comp}.err"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	plotCorrelation --corData {input} --corMethod pearson --whatToPlot heatmap --plotFile {output.pear} --skipZeros -min 0.8 -max 1 > {log.outpear} 2> {log.errpear}
	plotCorrelation --corData {input} --corMethod spearman --whatToPlot heatmap --plotFile {output.spear} --skipZeros -min 0.8 -max 1 > {log.outspear} 2> {log.errspear}
	'''

rule plotPCA:
	input:
		config['outDir'] + "/deepTools/{Comp}_BigwigSum.npz"
	output:
		config['outDir'] + "/Figures/{Comp}_PCA.png"
	log:
		out = config['outDir'] + "/logs/plotPCA.{Comp}.out",
		err = config['outDir'] + "/logs/plotPCA.{Comp}.err",
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	plotPCA --ntop 5000 --corData {input} -o {output} --transpose > {log.out} 2> {log.err}
	'''
# Differential accessibility.
rule edgeR:
	input:
		countMat = config['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat"
	output:
		sign = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.sign.tsv",
		allPeaks = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.all.tsv",
	log:
		out = config['outDir'] + "/logs/edgeR.{Comp}.out",
		err = config['outDir'] + "/logs/edgeR.{Comp}.err"
	params:
		scriptLoc = os.path.join(config["baseDir"], "Rscripts", "EdgeR.R"),
		condOrder = lambda wildcards, input: misc.conditionsfromCount(str(input.countMat) ,ss, wildcards.Comp),
		batchOrder = lambda wildcards, input: misc.batchesfromCount(str(input.countMat), ss, wildcards.Comp)
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	Rscript {params.scriptLoc} {input.countMat} {params.condOrder} {output.sign} {output.allPeaks} {params.batchOrder} > {log.out} 2> {log.err}
	'''

rule maPlot:
	input: 
		edgeR = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.all.tsv",
		countMat = config['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat"
	output: config['outDir'] + "/Figures/{Comp}_maPlot.png"
	threads: 1
	params: lambda wildcards: ss['CompCondDic'][wildcards.Comp]
	run:
		misc.plotter('maPlot',str(input.edgeR), str(output), params)

rule mergeDiff_Ann:
	input:
		annotation = config['outDir'] + "/Annotation/{Comp}_uropa_finalhits.txt",
		diffPeak = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.sign.tsv",
		nonSig = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.all.tsv"
	output:
		csvout = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated.sign.tsv",
		nonsigout = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated.all.tsv"
	threads: 1
	run:
		misc.mergeDiff_Ann(input.annotation, input.diffPeak, output.csvout)
		misc.mergeDiff_Ann(input.annotation, input.nonSig, output.nonsigout)

rule splitDiffRes:
	input:
		config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated.sign.tsv"
	output:
		UP = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv",
		DOWN = config['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_DOWN.tsv"
	threads: 1
	shell:'''
	head -1 {input} > {output.DOWN}
	# Down doesn't need to incorporate header because awk < 0 captures it.
	awk '$2 < 0' {input} >> {output.DOWN}
	awk '$2 > 0' {input} >> {output.UP}
	'''