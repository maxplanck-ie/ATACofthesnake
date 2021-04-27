import os
import shutil
from ATACofthesnake import misc

# Set some fallbacks if we don't have a sampleSheet.

if not config['sampleSheet']:
    ss = {}
    # Set the comparison to the output folder name (we don't have one).
    compStr = str(config['outDir'].split('/')[-1])
    ss['Comp'] = [compStr]
    ss[compStr] = config['Samples']


outList = []
if not config['peakSet']:
    # Create peaks.
    outList.append(
        expand(config['outDir'] + "/Figures/{CompCond}.mtFrac.png", CompCond = ss['Comp']) +
        expand(config['outDir'] + "/ShortBAM/{Sample}.bam.bai", Sample=config['Samples']) +
        expand(config['outDir'] + "/deepTools/{CompCond}.raw.fragSize.tsv", CompCond = ss['Comp'])
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
    shell:'''
    bamPEFragmentSize -b {params} -p {threads} --outRawFragmentLengths {output.raw} --table {output.table} > {log.out} 2> {log.err}
    '''