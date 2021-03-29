import os
import yaml
from ATACofthesnake import misc

# Read / set variables.
with open('Parameters.yaml') as f:
	paramDic = yaml.load(f, Loader=yaml.FullLoader)

localrules: fripPlotter, idxStatPlotter, maPlot
rule all:
	input:
		expand(os.path.join(paramDic['Loc']['bamDir'], '{sample}.bam.bai'), sample=paramDic['Samples']),
		paramDic['genomeFa'] + '.fai',
		expand(paramDic['Loc']['outDir'] + "/QC/{sample}.idxstat.txt", sample=paramDic['Samples']),
		expand(paramDic['Loc']['outDir'] + "/Figures/{Comp}.mtFrac.png", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam", sample=paramDic['Samples']),
		expand(paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam.bai", sample=paramDic['Samples']),
		expand(paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bed", sample=paramDic['Samples']),
		expand(paramDic['Loc']['outDir'] + "/deepTools/{Comp}.fragSizes.raw.tsv", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/MACS2/{sample}_peaks.narrowPeak", sample=paramDic['Samples']),
		expand(paramDic['Loc']['outDir'] + "/MACS2/{Comp}_Merged_peaks.narrowPeak", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/QC/{sample}.FRiP.txt", sample=paramDic['Samples']),
		expand(paramDic['Loc']['outDir'] + "/Figures/{Comp}.FRIP.png", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_scaleFactors.txt", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/BigWigs/{sample}.bw", sample=paramDic['Samples']),
		expand(paramDic['Loc']['outDir'] + '/deepTools/{Comp}_BigwigSum.npz', Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + '/Figures/{Comp}_plotCorr.png', Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/Figures/{Comp}_PCA.png", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/deepTools/{Comp}_BigwigMatrix.gz",Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/Figures/{Comp}_Heatmap.png", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.sign.tsv", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/Figures/{Comp}_maPlot.png", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/Annotation/{Comp}_uropa_finalhits.txt", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated.sign.tsv", Comp=paramDic['Comp']),
		expand(paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv", Comp=paramDic['Comp'])

rule checkGenomeIndex:
	input: paramDic['genomeFa']
	output: paramDic['genomeFa'] + '.fai'
	log:
		out = paramDic['Loc']['outDir'] + "/logs/checkFai.out",
		err = paramDic['Loc']['outDir'] + "/logs/checkFai.err"
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs', 'AOS_SeqTools.yaml')
	shell:'''
	samtools faidx {input}
	'''

rule checkIndex:
	input:
		sample = paramDic['Loc']['bamDir'] + '/{sample}.bam',
		genomeFai = paramDic['genomeFa'] + '.fai'
	output:
		paramDic['Loc']['bamDir'] + '/{sample}.bam.bai'
	log:
		out = paramDic['Loc']['outDir'] + '/logs/checkIndex.{sample}.out',
		err = paramDic['Loc']['outDir'] + '/logs/checkIndex.{sample}.err'
	threads: 10
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	sambamba index -t {threads} {input.sample} > {log.out} 2> {log.err}
	'''

rule idxStat:
	input:
		index = paramDic['Loc']['bamDir'] + '/{sample}.bam.bai',
		bam = paramDic['Loc']['bamDir'] + '/{sample}.bam'
	output:
		paramDic['Loc']['outDir'] + "/QC/{sample}.idxstat.txt"
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	set +o pipefail;
	samtools idxstats {input.bam} | cut -f1,3 > {output}
	'''

rule alignmentSieve:
	input:
		inBam = paramDic['Loc']['bamDir'] + '/{sample}.bam',
		index = paramDic['Loc']['bamDir'] + '/{sample}.bam.bai',
		mitoBleed = paramDic['Loc']['outDir'] + "/QC/{sample}.idxstat.txt"
	output:
		shortBam = paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam",
		filterMetrics = paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.metrics"
	params:
		blackList = paramDic['blackList'],
		fragSize = paramDic['fragSize']
	log:
		out = paramDic['Loc']['outDir'] + '/logs/alignmentSieve.{sample}.out',
		err = paramDic['Loc']['outDir'] + '/logs/alignmentSieve.{sample}.err'
	threads:10
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	alignmentSieve --bam {input.inBam} --outFile {output.shortBam} -p {threads} --filterMetrics {output.filterMetrics} --maxFragmentLength {params.fragSize} --minFragmentLength 0 --blackListFileName {params.blackList} > {log.out} 2> {log.err}
	'''

rule shortIndex:
	input:
		paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam"
	output:
		paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam.bai"
	log:
		out = paramDic['Loc']['outDir'] + '/logs/shortIndex.{sample}.out',
		err = paramDic['Loc']['outDir'] + '/logs/shortIndex.{sample}.err',
	threads: 10
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	sambamba index {input} > {log.out} 2> {log.err}
	'''

rule fragSize:
	input:
		expand(paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam.bai", sample=paramDic['Samples'])
	output:
		raw = paramDic['Loc']['outDir'] + "/deepTools/{Comp}.fragSizes.raw.tsv",
		table = paramDic['Loc']['outDir'] + "/deepTools/{Comp}.fragSizes.metrics.tsv"
	log:
		out = paramDic['Loc']['outDir'] + '/logs/fragSize.{Comp}.out',
		err = paramDic['Loc']['outDir'] + '/logs/fragSize.{Comp}.err'
	params:
		lambda wildcards: ' '.join(expand(paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam", sample=paramDic['Comp'][wildcards.Comp]['Samples']))
	threads: 10
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	bamPEFragmentSize -b {params} -p {threads} --outRawFragmentLengths {output.raw} --table {output.table} > {log.out} 2> {log.err}
	'''

rule bamToBed:
	input:
		bai = paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam.bai",
		bam = paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam"
	output:
		outBed = paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bed"
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	bamToBed -i {input.bam} > {output.outBed}
	'''

rule MACS2:
	input:
		paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bed"
	output:
		paramDic['Loc']['outDir'] + "/MACS2/{sample}_peaks.narrowPeak"
	log:
		out = paramDic['Loc']['outDir'] + '/logs/MACS2.{sample}.out',
		err = paramDic['Loc']['outDir'] + '/logs/MACS2.{sample}.err'
	params:
		genomeSize = paramDic['genomeSize'],
		outName = lambda wildcards: wildcards.sample,
		blackList = paramDic['blackList']
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	macs2 callpeak -t {input} -f BED --nomodel --shift -75 --extsize 150 -g {params.genomeSize} -n {params.outName} -q 0.01 --outdir AOS/MACS2/ --keep-dup all > {log.out} 2> {log.err}
	'''

rule mergePeak:
	input:
		expand(paramDic['Loc']['outDir'] + "/MACS2/{sample}_peaks.narrowPeak", sample=paramDic['Samples'])
	output:
		paramDic['Loc']['outDir'] + "/MACS2/{Comp}_Merged_peaks.narrowPeak"
	params:
		lambda wildcards: ' '.join(expand(paramDic['Loc']['outDir'] +"/MACS2/{sample}_peaks.narrowPeak", sample=paramDic['Comp'][wildcards.Comp]['Samples']))
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	echo {params}
	cat {params} | sort -k1,1 -k2,2n | bedtools merge > {output}
	'''


rule fripScore:
	input:
		bamfile = paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam",
		MACS2done = expand(paramDic['Loc']['outDir'] + "/MACS2/{Comp}_Merged_peaks.narrowPeak", Comp=paramDic['Comp'])
	output:
		paramDic['Loc']['outDir'] + "/QC/{sample}.FRiP.txt"
	params:
		genomeSize = paramDic['genomeSize'],
		peakFile = lambda wildcards: misc.returnCompfromSample(wildcards.sample ,paramDic),
		sample = lambda wildcards: wildcards.sample
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
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
		fripFiles = expand(paramDic['Loc']['outDir'] + "/QC/{sample}.FRiP.txt", sample=paramDic['Samples'])
	output:
		paramDic['Loc']['outDir'] + "/Figures/{Comp}.FRIP.png"
	threads: 1
	params: lambda wildcards: paramDic['Comp'][wildcards.Comp]['Samples']
	run:
		misc.plotter('frip', params, str(output))

rule idxStatPlotter:
	input: 
		idxstatFiles = expand(paramDic['Loc']['outDir'] + "/QC/{sample}.idxstat.txt", sample=paramDic['Samples'])
	output:
		paramDic['Loc']['outDir'] + "/Figures/{Comp}.mtFrac.png"
	threads: 1
	params: lambda wildcards: paramDic['Comp'][wildcards.Comp]['Samples']
	run:
		misc.plotter('idxstat', params, str(output))

rule countMat:
	input:
		paramDic['Loc']['outDir'] + "/MACS2/{Comp}_Merged_peaks.narrowPeak"
	output:
		mat = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat",
		matnpz = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.npz"
	log:
		out = paramDic['Loc']['outDir'] + "/logs/countMat.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/countMat.{Comp}.err"
	params:
		blackList = paramDic['blackList'],
		samples = lambda wildcards: ' '.join(expand(paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam", sample=paramDic['Comp'][wildcards.Comp]['Samples']))
	threads: 20
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	multiBamSummary BED-file --BED {input} -bl {params.blackList} -p {threads} --outRawCounts {output.mat} -o {output.matnpz} -b {params.samples} > {log.out} 2> {log.err}
	# Remove quotes placed by multiBamSummary
	sed -i "s/'//g" {output.mat}
	# Remove .bam postfix
	sed -i 's/\.bam//g' {output.mat}
	'''

rule scaleFactors:
	input:
		paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat"
	output:
		paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_scaleFactors.txt"
	log:
		out = paramDic['Loc']['outDir'] + "/logs/scaleFactors.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/scaleFactors.{Comp}.err"
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	params:
		scriptLoc = os.path.join(paramDic["baseDir"], "Rscripts", "edgeR_scalefactors.R")
	shell:'''
	Rscript {params.scriptLoc} {input} {output} > {log.out} 2> {log.err}
	'''

rule BigWigs:
	input:
		sf = expand(paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_scaleFactors.txt",Comp=paramDic['Comp']),
		inFile = paramDic['Loc']['outDir'] + "/ShortBAM/{sample}.bam"
	output:
		paramDic['Loc']['outDir'] + "/BigWigs/{sample}.bw"
	log:
		out = paramDic['Loc']['outDir'] + "/logs/BigWigs.{sample}.out",
		err = paramDic['Loc']['outDir'] + "/logs/Bigwigs.{sample}.err"
	params:
		sampleName = '{sample}',
		blackList = paramDic['blackList'],
		genomeSize = paramDic['genomeSize']
	threads: 10
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	SCALEFAC=$(grep {params.sampleName} {input.sf} | cut -f2 -d ' ')
	bamCoverage --scaleFactor $SCALEFAC -b {input.inFile} -o {output} -p {threads} -bs 25 --extendReads --ignoreDuplicates -bl {params.blackList} > {log.out} 2> {log.err}
	'''

rule multiBigwigSum:
	input:
		expand(paramDic['Loc']['outDir'] + '/BigWigs/{sample}.bw', sample=paramDic['Samples']),
		Peaks = paramDic['Loc']['outDir'] + "/MACS2/{Comp}_Merged_peaks.narrowPeak"
	output:
		paramDic['Loc']['outDir'] + '/deepTools/{Comp}_BigwigSum.npz'
	log:
		out = paramDic['Loc']['outDir'] + "/logs/multiBigWigSum.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/multiBigWigSum.{Comp}.err"
	params:
		bigwigs = lambda wildcards: ' '.join(expand(paramDic['Loc']['outDir'] + '/BigWigs/{sample}.bw', sample=paramDic['Comp'][wildcards.Comp]['Samples'])),
		blackList = paramDic['blackList']
	threads: 10
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	multiBigwigSummary BED-file --BED {input.Peaks} -o {output} -b {params.bigwigs} -bl {params.blackList} -p {threads} -bs 150 > {log.out} 2> {log.err}
	'''

rule plotCorr:
	input:
		paramDic['Loc']['outDir'] + "/deepTools/{Comp}_BigwigSum.npz"
	output:
		paramDic['Loc']['outDir'] + "/Figures/{Comp}_plotCorr.png"
	log:
		out = paramDic['Loc']['outDir'] + "/logs/plotCorr.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/plotCorr.{Comp}.err"
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	plotCorrelation --corData {input} --corMethod pearson --whatToPlot heatmap --plotFile {output} --skipZeros -min 0.8 -max 1 > {log.out} 2> {log.err}
	'''

rule plotPCA:
	input:
		paramDic['Loc']['outDir'] + "/deepTools/{Comp}_BigwigSum.npz"
	output:
		paramDic['Loc']['outDir'] + "/Figures/{Comp}_PCA.png"
	params: 
		colors = lambda wildcards: " ".join(sum([[s] * n for s, n in zip(["red","blue"], [len(paramDic['Comp'][wildcards.Comp]['Cond'][paramDic['Cond'][0]]),len(paramDic['Comp'][wildcards.Comp]['Cond'][paramDic['Cond'][1]])])], []))
	log:
		out = paramDic['Loc']['outDir'] + "/logs/plotPCA.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/plotPCA.{Comp}.err",
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	plotPCA --corData {input} -o {output} --transpose --colors {params.colors} > {log.out} 2> {log.err}
	'''

rule computeMatrix:
	input:
		expand(paramDic['Loc']['outDir'] + '/BigWigs/{sample}.bw', sample=paramDic['Samples']),
	output:
		paramDic['Loc']['outDir'] + "/deepTools/{Comp}_BigwigMatrix.gz"
	log:
		out = paramDic['Loc']['outDir'] + "/logs/computeMatrix.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/computeMatrix.{Comp}.err"
	params:
		bigwigs = lambda wildcards: ' '.join(expand(paramDic['Loc']['outDir'] + '/BigWigs/{sample}.bw', sample=paramDic['Comp'][wildcards.Comp]['Samples'])),
		bed = "TSS.bed"
	threads: 10
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	computeMatrix reference-point -S {params.bigwigs} -R {params.bed} --referencePoint center -a 2000 -b 2000 -out {output} -p {threads} --missingDataAsZero > {log.out} 2> {log.err}
	'''

rule plotHeatmap:
	input:
		paramDic['Loc']['outDir'] + "/deepTools/{Comp}_BigwigMatrix.gz"
	output:
		paramDic['Loc']['outDir'] + "/Figures/{Comp}_Heatmap.png"
	log:
		out = paramDic['Loc']['outDir'] + "/logs/plotHeatmap.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/plotHeatmap.{Comp}.err"
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	plotHeatmap -m {input} -out {output} --refPointLabel TSS > {log.out} 2> {log.err}
	'''

rule edgeR:
	input:
		countMat = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat"
	output:
		sign = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.sign.tsv",
		allPeaks = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.all.tsv",
	log:
		out = paramDic['Loc']['outDir'] + "/logs/edgeR.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/edgeR.{Comp}.err"
	params:
		scriptLoc = os.path.join(paramDic["baseDir"], "Rscripts", "EdgeR.R"),
		condOrder = lambda wildcards, input: misc.conditionsfromCount(str(input.countMat) ,paramDic['Comp'][wildcards.Comp]['Cond']),
		batchOrder = lambda wildcards, input: misc.batchesfromCount(str(input.countMat), paramDic) if paramDic['batchStatus'] == 1 else 'None'
	threads: 1
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	Rscript {params.scriptLoc} {input.countMat} {params.condOrder} {output.sign} {output.allPeaks} {params.batchOrder} > {log.out} 2> {log.err}
	'''

rule maPlot:
	input: 
		edgeR = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.all.tsv",
		countMat = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_counts.mat"
	output: paramDic['Loc']['outDir'] + "/Figures/{Comp}_maPlot.png"
	threads: 1
	params: lambda wildcards: list(paramDic['Comp'][wildcards.Comp]['Cond'].keys())
	run:
		misc.plotter('maPlot',str(input.edgeR), str(output), params)

rule uropa:
	input:
		paramDic['Loc']['outDir'] + "/MACS2/{Comp}_Merged_peaks.narrowPeak"
	output:
		paramDic['Loc']['outDir'] + "/Annotation/{Comp}_uropa_finalhits.txt"
	log:
		out = paramDic['Loc']['outDir'] + "/logs/uropa.{Comp}.out",
		err = paramDic['Loc']['outDir'] + "/logs/uropa.{Comp}.err"
	params:
		GTF = paramDic['GTF'],
		prefix = "{Comp}_uropa"
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	threads: 5
	shell:'''
	uropa -b {input} -g {params.GTF} --summary --feature transcript --distance 10000 --internals 1 -p {params.prefix} -o "AOS/Annotation" -t {threads} --show-attributes gene_id transcript_id gene_name gene_type transcript_type > {log.out} 2> {log.err}
	'''

rule mergeDiff_Ann:
	input:
		annotation = paramDic['Loc']['outDir'] + "/Annotation/{Comp}_uropa_finalhits.txt",
		diffPeak = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.sign.tsv",
		nonSig = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR.all.tsv"
	output:
		csvout = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated.sign.tsv",
		nonsigout = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated.all.tsv"
	threads: 1
	run:
		misc.mergeDiff_Ann(input.annotation, input.diffPeak, output.csvout)
		misc.mergeDiff_Ann(input.annotation, input.nonSig, output.nonsigout)

rule splitDiffRes:
	input:
		paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated.sign.tsv"
	output:
		UP = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv",
		DOWN = paramDic['Loc']['outDir'] + "/diffAcc_{Comp}/{Comp}_edgeR_annotated_DOWN.tsv"
	threads: 1
	shell:'''
	head -1 {input} > {output.DOWN}
	# Down doesn't need to incorporate header because awk < 0 captures it.
	awk '$2 < 0' {input} >> {output.DOWN}
	awk '$2 > 0' {input} >> {output.UP}
	'''
