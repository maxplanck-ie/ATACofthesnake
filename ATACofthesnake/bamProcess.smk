import os
import shutil
from ATACofthesnake import misc

# Create a sample dict for merging, diff analysis etc.
#if config['sampleSheet']:
#    ss = {}
    # Define comp
    # Define config[Comp]
if not config['sampleSheet']:
    ss = {}
    # Set the comparison to the output folder name (we don't have one).
    compStr = str(config['outDir'].split('/')[-1])
    ss['Comp'] = [compStr]
    ss[compStr] = config['Samples']
print(ss)

def ruleFetcher(ss):
    outList = []
    # checkGenomeIndex, checkIndex, idxStat, idxStatPlotter
    outList.append(
        expand(config['outDir'] + "/Figures/{CompCond}.mtFrac.png", CompCond = ss['Comp'])
    )
    # alignmentSieve
    outList.append(
        expand(config['outDir'] + "/ShortBAM/{Sample}.bam", Sample=config['Samples'])
    )
    if not config['sampleSheet']:
        if config['mergeBam']:
            # No sampleSheet, merge all bamFiles before peak calling.
            outList.append(
                expand(config['outDir'] + "/MergeBAM/{CompCond}.bam", CompCond = ss['Comp'])
            )
            outList.append(
                expand(config['outDir'] + "/MergeBAM/{CompCond}.bam", CompCond = ss['Comp'])
            )
            outList.append(
                expand(config['outDir'] + "/MergeBAM/{CompCond}.bed", CompCond = ss['Comp'])
            )
            outList.append(
                expand(config['outDir'] + "/MACS2/{CompCond}_peaks.narrowPeak", CompCond=ss['Comp'])
            )
            return outList
        if not config['mergeBam']:
            print('nomerge')
    if config['sampleSheet']:
        if config['mergeBam']:
            print(mergeSS)
        if not config['mergeBam']:
            print(noMergeSS)


rule all:
    input:
        ruleFetcher(ss)


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
		lambda wildcards: expand(config['outDir'] + "/QC/{Sample}.idxstat.txt", Sample=ss[wildcards.Comp])
	output:
		config['outDir'] + "/Figures/{Comp}.mtFrac.png"
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

rule mergeBam:
	input:
		expand(config['outDir'] + "/ShortBAM/{Sample}.bam", Sample=config['Samples'])
	output:
		bam = config['outDir'] + "/MergeBAM/{CompCond}.bam"
	params:
		lambda wildcards: ' '.join(expand(config['outDir'] + "/ShortBAM/{Sample}.bam", Sample=ss[wildcards.CompCond]))
	threads: 5
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	samtools merge -@ {threads} {output.bam} {params}
	samtools index -@ {threads} {output.bam}
	'''

rule mergeBam_to_bed:
	input:
		bai = config['outDir'] + "/MergeBAM/{CompCond}.bam.bai",
		bam = config['outDir'] + "/MergeBAM/{CompCond}.bam"
	output:
		outBed = config['outDir'] + "/MergeBAM/{CompCond}.bed"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	bamToBed -i {input.bam} > {output.outBed}
	'''
rule mergeMACS2:
	input:
		config['outDir'] + "/MergeBAM/{CompCond}.bed"
	output:
		config['outDir'] + "/MACS2/{CompCond}_peaks.narrowPeak"
	log:
		out = config['outDir'] + '/logs/MACS2.{CompCond}.out',
		err = config['outDir'] + '/logs/MACS2.{CompCond}.err'
	params:
		genomeSize = config['genomeSize'],
		outName = lambda wildcards: wildcards.CompCond,
		blackList = config['blackList'],
		outDir = config['outDir'] + "/MACS2"
	threads: 1
	conda: os.path.join(config['baseDir'], 'envs','AOS_SeqTools.yaml')
	shell:'''
	macs2 callpeak -t {input} -f BED --nomodel --shift -75 --extsize 150 -g {params.genomeSize} -n {params.outName} -q 0.01 --outdir {params.outDir} --keep-dup all > {log.out} 2> {log.err}
	'''