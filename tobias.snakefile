from os.path import join
import sys
import os
import pandas as pd
import yaml
import glob

include: join("rules/init.smk")

localrules: all

rule all:
	input:
		expand(join(WORKDIR, "bams", "{condition}.bam"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "coverage", "{condition}_coverage.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_uncorrected.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_bias.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_expected.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_corrected.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "footprinting", "{condition}_footprints.bw"), condition=CONDITION_IDS),
		# expand(join(WORKDIR, "overview", "all_{condition}_bound.bed"), condition=CONDITION_IDS),

rule condition_bam:
	input:
		lambda wildcards: config["data"][wildcards.condition]
	output:
		bam = join(WORKDIR,"bams","{condition}.bam"),
	envmodules: TOOLS["sambamba"]["version"]
	threads: 16
	message:
		"Running {rule} for condition: {wildcards.condition}"
	params:
		outdir=join(WORKDIR,"bams"),
		mem=MEMORY,
		condition="{condition}",
	shell:"""
sorted_bams=""
for bam in {input}
do
bn=$(basename $bam)
sambamba sort --memory-limit={params.mem}G --tmpdir=/lscratch/${{SLURM_JOBID}} --nthreads={threads} --out=/lscratch/${{SLURM_JOBID}}/${{bn%.*}}.sorted.bam $bam
sorted_bams="$sorted_bams /lscratch/${{SLURM_JOBID}}/${{bn%.*}}.sorted.bam"
done
sambamba merge --nthreads={threads} {output.bam} $sorted_bams
"""

rule coverage_bw:
	input:
		bam = rules.condition_bam.output.bam
	output:
		coveragebw = join(WORKDIR, "coverage", "{condition}_coverage.bw")
	params:
		condition = "{condition}"
	message:
		"Running {rule} for condition: {wildcards.condition} ({input.bam})"
	threads: 4
	envmodules: TOOLS["bedtools"]["version"],TOOLS["samtools"]["version"], TOOLS["ucsc"]["version"]
	shell:"""
bedtools genomecov -ibam {input.bam} -bg >  /lscratch/${{SLURM_JOBID}}/{params.condition}.bg
bedSort /lscratch/${{SLURM_JOBID}}/{params.condition}.bg /lscratch/${{SLURM_JOBID}}/{params.condition}.bg
samtools view -H {input.bam}|grep '^@SQ'|cut -f2,3|sed 's/SN:\|LN://g' > /lscratch/${{SLURM_JOBID}}/{params.condition}.chrom.sizes
bedGraphToBigWig /lscratch/${{SLURM_JOBID}}/{params.condition}.bg /lscratch/${{SLURM_JOBID}}/{params.condition}.chrom.sizes {output.coveragebw}
"""

rule atacorrect:
	input:
		bam = rules.condition_bam.output.bam, 	#os.path.join(OUTPUTDIR, "mapping", "{condition}.bam"),
		peaks = PEAKS, 	#os.path.join(OUTPUTDIR, "peak_calling", "all_merged.bed"),
		genome = REFFA,
	output:
		uncorrected = join(WORKDIR, "bias_correction", "{condition}_uncorrected.bw"),
		bias = join(WORKDIR, "bias_correction", "{condition}_bias.bw"),
		expected = join(WORKDIR, "bias_correction", "{condition}_expected.bw"),
		corrected = join(WORKDIR, "bias_correction", "{condition}_corrected.bw"),
	params: 
		blacklist = BLACKLIST,
		outdir = join(WORKDIR, "bias_correction"),
		condition = "{condition}",
		autocorrect_extra_params = AUTOCORRECT_EXTRA_PARAMS
	threads: 
		56 	#unless there are more than 99 cores, this rule will run on max threads
	message: 
		"Running {rule} for condition: {wildcards.condition} ({input.bam})"
	container: "docker://nciccbr/ccbr_tobias:latest"
	shell:"""
TOBIAS ATACorrect \
	-b {input.bam} \
	-g {input.genome} \
	-p {input.peaks} \
	--cores {threads} \
	--prefix {params.condition} \
	--outdir {params.outdir} \
	{params.blacklist} {params.autocorrect_extra_params}
"""

rule footprinting:
	input: 
		signal = rules.atacorrect.output.corrected, 	#os.path.join(OUTPUTDIR, "bias_correction", "{condition}_corrected.bw"),
		regions = PEAKS		#os.path.join(OUTPUTDIR, "peak_calling", "all_merged.bed")
	output: 
		footprints = join(WORKDIR, "footprinting", "{condition}_footprints.bw"),
	params:
		footprinting_extra_params = FOOTPRINTING_EXTRA_PARAMS
	threads: 
		56
	message: 
		"Running {rule} for condition: {wildcards.condition} ({input.signal})"
	container: "docker://nciccbr/ccbr_tobias:latest"
	shell:"""
TOBIAS FootprintScores \
	--signal {input.signal} \
	--regions {input.regions} \
	--output {output.footprints} \
	--cores {threads} {params.footprinting_extra_params}
"""

#Estimate bound sites from scored file
checkpoint bindetect:
	input: 
		motifs = MOTIFS, 		
		footprints = expand(rules.footprinting.output.footprints, condition=CONDITION_IDS),
		genome = REFFA,
		peaks = config["annotatedpeaks"],
		peak_header = config["annotatedpeaksheader"]
	output:
		directory(join(WORKDIR, "TFBS"))
	threads: 56	
	params:
		conditions="--cond_names " + " ".join(CONDITION_IDS),	#comma inserts space between elements
		bindetect_extra_params=BINDETECT_EXTRA_PARAMS,
		workdir=WORKDIR
	message: 
		"Running {rule}"
	container: "docker://nciccbr/ccbr_tobias:latest"
	shell:"""
TOBIAS BINDetect \
	--motifs {input.motifs} \
	--signals {input.footprints} \
	--genome {input.genome} \
	--peaks {input.peaks} \
	--peak_header {input.peak_header} \
	--cores {threads} \
	--outdir {output} \
	{params.conditions}
mkdir -p {params.workdir}/overview
cp {output}/*.txt {params.workdir}/overview/
cp {output}/*.xlsx {params.workdir}/overview/
cp {output}/*.pdf {params.workdir}/overview/
"""
#--------------------------------------------------------------------------------------------------------#

def get_TF_ids(wildcards):
	bindetect_output = checkpoints.bindetect.get(**wildcards).output[0]
	TF_IDS = glob_wildcards(os.path.join(bindetect_output, "{TF}", "beds")).TF	#wildcard for each TF from dir name
	return(TF_IDS)

#--------------------------------------------------------------------------------------------------------#

#Join bound estimates per condition
rule join_bound:
	input:
		lambda wildcards: expand(join(WORKDIR, "TFBS", "{TF}", "beds", "{TF}_{{condition}}_bound.bed"), TF=get_TF_ids(wildcards)),
	output:
		unsorted = temp(join(WORKDIR, "overview", "all_{condition}_bound.tmp")),
		final = join(WORKDIR, "overview", "all_{condition}_bound.bed")
	threads: 56
	shell:"""
touch {output.unsorted}
touch {output.final}
"""
	# run:
	# 	n = 100 	#chunks of 100 files
		
	# 	shell("> " + output[0])
	# 	file_chunks = [input[i:i+n] for i in range(0,len(input),n)]
	# 	for chunk in file_chunks:
	# 		shell("cat {0} | bedtools sort >> {1}".format(" ".join(chunk), output.unsorted))
	# 	shell("sort -k 1,1 -k2,2 -n {output.unsorted} > {output.final}")
	# 	shell("igvtools index {output.final}")