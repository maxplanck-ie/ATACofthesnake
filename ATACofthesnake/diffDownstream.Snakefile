import os
import yaml
from ATACofthesnake import misc

# Read / set variables.
with open('Parameters.yaml') as f:
	paramDic = yaml.load(f, Loader=yaml.FullLoader)
Conditions = list(paramDic['Cond'])

# Does it make sense to run this workflow?

rule diffHeat:
	input:
		up = "diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv",
		down = "diffAcc_{Comp}/{Comp}_edgeR_annotated_DOWN.tsv",
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
	cut -f1 {input.up} | tail -n +2 | tr '_' '\t' > {params.Comp}.UP.bed
	cut -f1 {input.down} | tail -n +2 | tr '_' '\t' > {params.Comp}.DOWN.bed
	upcount=$(wc -l < {params.Comp}.UP.bed)
	downcount=$(wc -l < {params.Comp})
	computeMatrix reference-point -S {params.bigwigs} -R {params.Comp}.UP.bed {params.Comp}.DOWN.bed --referencePoint center -a 2000 -b 2000 -out {output.outMat} -p {threads} --missingDataAsZero > {log.Matout} 2> {log.Materr}
	rm {params.Comp}.UP.bed
	rm {params.Comp}.DOWN.bed
	plotHeatmap -m {output.outMat} -out {output.heatmap} --refPointLabel center --colorMap Blues > {log.Plotout} 2> {log.Ploterr}
	'''

rule formatFNA:
	input:
		up = "diffAcc_{Comp}/{Comp}_edgeR_annotated_UP.tsv",
		down = "diffAcc_{Comp}/{Comp}_edgeR_annotated_DOWN.tsv",
		fai = paramDic['genomeFa'] + '.fai'
	output:
		up = 'Motif/{Comp}_up.fna',
		down = 'Motif/{Comp}_down.fna'
	threads: 1
	params:
		upBed = "Motif/{Comp}_up.bed",
		downBed = "Motif/{Comp}_down.bed",
		genomeFa = paramDic['genomeFa']
	conda: os.path.join(paramDic['baseDir'], 'envs','AOS_SeqTools.yaml')
	log:
		out = 'logs/formatFNA_{Comp}.out',
		err = 'logs/formatFNA_{Comp}.err'
	shell:'''
	set +o pipefail;cut -f1 {input.up} | tail -n +2 | tr '_' '\t' > {params.upBed}
	set +o pipefail;cut -f1 {input.down} | tail -n +2 | tr '_' '\t' > {params.downBed}
	bedtools getfasta -fi {params.genomeFa} -bed {params.upBed} > {output.up}
	bedtools getfasta -fi {params.genomeFa} -bed {params.downBed} > {output.down}
	'''

#rule runMeme:


rule TOBIAS:
	input:
		"MACS2/{Comp}_Merged_peaks.narrowPeak"
	output:
		"TOBIAS_{Comp}/footprints.txt"
	threads: 1
	shell:'''
	touch {output}
	'''