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
            expand(config['outDir'] + "/MACS2_mergeBAM/{CompCond}_peaks.bed", CompCond=ss['CompCond']) +
            expand(config['outDir'] + '/MACS2_mergeBAM/{Comp}_peaks.bed', Comp=ss['Comp'])
        )
else:
    outList.append(
        expand(config['outDir'] + "/MACS2/{CompCond}_peaks.bed", CompCond=ss['Comp'])
    )

print(outList)

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