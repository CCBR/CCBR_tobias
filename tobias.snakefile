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
		expand(os.path.join(WORKDIR, "coverage", "{condition}_coverage.bw"), condition=CONDITION_IDS),
		# expand(os.path.join(WORKDIR, "footprinting", "{condition}_footprints.bw"), condition=CONDITION_IDS),
		# expand(os.path.join(WORKDIR, "overview", "all_{condition}_bound.bed"), condition=CONDITION_IDS),

rule coverage_bw:
	input:
		lambda wildcards: config["data"][wildcards.condition]
	output:
		bam = join(WORKDIR,"bams","{condition}.bam"),
		coveragebw = join(WORKDIR, "coverage", "{condition}_coverage.bw")
	envmodules: TOOLS["samtools"]["version"], TOOLS["sambamba"]["version"], TOOLS["ucsc"]["version"]
	threads: 4
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
bedtools genomecov -ibam {output.bam} -bg | sort -k1,1 -k2,2n >  /lscratch/${{SLURM_JOBID}}/${{bn%.*}}.sorted.bg
samtools view -H D4_Neural_Nkx1_2_Dox_1.qsorted.bam|grep '^@SQ'|cut -f2,3|sed 's/[SN:\|LN:]//g' > /lscratch/${{SLURM_JOBID}}/${{bn%.*}}.sorted.chrom.sizes
bedGraphToBigWig /lscratch/${{SLURM_JOBID}}/${{bn%.*}}.sorted.bg /lscratch/${{SLURM_JOBID}}/${{bn%.*}}.sorted.chrom.sizes {output.coveragebw}
"""
