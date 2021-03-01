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
		# expand(os.path.join(WORKDIR, "footprinting", "{condition}_footprints.bw"), condition=CONDITION_IDS),
		# expand(os.path.join(WORKDIR, "overview", "all_{condition}_bound.bed"), condition=CONDITION_IDS),

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
samtools view -H {input.bam}|grep '^@SQ'|cut -f2,3|sed 's/[SN:\|LN:]//g' > /lscratch/${{SLURM_JOBID}}/{params.condition}.chrom.sizes
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