from os.path import join
import sys
import os
import pandas as pd
import yaml
import glob

include: join("rules/init.smk")
print(CONTRASTS)

def get_bound_bed_list(wildcards) :
	files=dict()
	for cont,conds in CONTRASTS2CONDITIONS.items():
		for cond in conds:
			files[cont+"."+cond]=join(WORKDIR, "overview_"+cont, "all_"+cond+"_bound.bed.gz")
	return files

def get_bound_bed_list_per_contrast_per_condition_per_TF(wildcards) :
	files=dict()
	for cont,conds in CONTRASTS2CONDITIONS.items():
		for cond in conds:
			for tf in TFs:
				files[cont+"."+cond+"."+TF]=join(WORKDIR, "TFBS_"+cont, tf, "beds", tf+"_"+cond+"_bound.bed")
	return files

localrules: all, plot_heatmaps_aggregates

rule all:
	input:
		expand(join(WORKDIR, "bams", "{condition}.bam"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "coverage", "{condition}_coverage.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_uncorrected.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_bias.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_expected.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "bias_correction", "{condition}_corrected.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "footprinting", "{condition}_footprints.bw"), condition=CONDITION_IDS),
		expand(join(WORKDIR, "TFBS_{contrast}", "bindetect_figures.pdf"),contrast=CONTRASTS),
		expand(join(WORKDIR,"TFBS_{contrast}","bound_beds_list.tsv"),contrast=CONTRASTS),
		# unpack(get_bound_bed_list),
		# unpack(get_bound_bed_list_per_contrast_per_condition_per_TF),
		expand(join(WORKDIR, "overview_{cont}", "all_{cond}_bound.bed.gz"),zip,cont=CC1,cond=CC2),
		expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.heatmap.pdf"),contrast=CONTRASTS,TF=TFs),
		expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.aggregate.pdf"),contrast=CONTRASTS,TF=TFs),
		# join(WORKDIR, "footprinting","dummy"),
		# unpack(get_bound_bed_list),	
		# expand(join(WORKDIR, "overview", "all_{condition}_bound.bed.gz"), condition=CONDITION_IDS),



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

def c2c(wildcards):
	files=[]
	for c in CONTRASTS2CONDITIONS[wildcards.contrast]:
		files.append(join(WORKDIR, "footprinting", c+"_footprints.bw"))
	return files

rule bindetect:
	input:		
		lambda wildcards: unpack(c2c(wildcards))
	output:
		pdf=join(WORKDIR,"TFBS_{contrast}","bindetect_figures.pdf"),
		# beds=expand(join(WORKDIR, "TFBS_{{contrast}}", "{TF}", "beds", "{TF}_{{condition}}_bound.bed"), TF=TFs)
	params:
		contrast="{contrast}",
		motifs = MOTIFS, 		
		genome = REFFA,
		peaks = config["annotatedpeaks"],
		peak_header = config["annotatedpeaksheader"],
		bindetect_extra_params=BINDETECT_EXTRA_PARAMS,
	threads: 56
	container: "docker://nciccbr/ccbr_tobias:latest"
	message: 
		"Running {rule} for {wildcards.contrast}"
	shell:"""
c1=$(echo {params.contrast}|awk -F\"_vs_\" '{{print $1}}')
c2=$(echo {params.contrast}|awk -F\"_vs_\" '{{print $2}}')
outdir=$(dirname "{output.pdf}")
TOBIAS BINDetect \
	--motifs {params.motifs} \
	--signals {input} \
	--genome {params.genome} \
	--peaks {params.peaks} \
	--peak-header {params.peak_header} \
	--cores {threads} \
	--cond-names $c1 $c2 \
	--outdir $outdir \
	{params.bindetect_extra_params}
"""

rule create_cont_cond_tf_join_filelist:
	input: 
		rules.bindetect.output.pdf
	output: 
		tsv=join(WORKDIR,"TFBS_{contrast}","bound_beds_list.tsv"),
	run:
		cont=wildcards.contrast
		conds=CONTRASTS2CONDITIONS[cont]
		out=open(output[0], "w")
		for cond in conds:
			for tf in TFs:
				f=join(WORKDIR, "TFBS_"+cont, tf, "beds", tf+"_"+cond+"_bound.bed")
				out.write("%s\t%s\n"%(cond,f))
		out.close()

# Join bound estimates per condition
rule join_bound:
	input:
		rules.create_cont_cond_tf_join_filelist.output.tsv
	output:
		bedgz = join(WORKDIR, "overview_{contrast}", "all_{cond}_bound.bed.gz")
	threads: 2
	envmodules: TOOLS["ucsc"]["version"], TOOLS["samtools"]["version"]
	params:
		cont="{contrast}",
		cond="{cond}"
	shell:"""
x="{output.bedgz}"
xbed="${{x%.*}}"
while read a b;do 
	if [ "$a" == "{params.cond}" ];then
		cat $b
	fi
done < {input} > $xbed
bedSort $xbed $xbed
bgzip -f --threads {threads} $xbed
tabix -f -p bed $x
"""

rule plot_heatmaps_aggregates:
	input:
		rules.bindetect.output.pdf
	output:
		heatmap = join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.heatmap.pdf"),
		aggregate = join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.aggregate.pdf")
	params:
		contrast = "{contrast}",
		tf = "{TF}",
		workdir = WORKDIR
	container: "docker://nciccbr/ccbr_tobias:latest"
	message: 
		"Running {rule} for Contrast:{wildcards.contrast} and TF:{wildcards.TF}"
	shell:"""
contrast="{params.contrast}"
TF="{params.tf}"
workdir="{params.workdir}"
c1=$(echo {params.contrast}|awk -F\"_vs_\" '{{print $1}}')
c2=$(echo {params.contrast}|awk -F\"_vs_\" '{{print $2}}')
x=$(echo "${{TF}}"|awk -F"_" '{{print NF/2}}')
tf=$(echo "${{TF}}"|sed "s/_/\\t/g"|cut -f1,$x|sed "s/\\t/_/g")
TOBIAS PlotHeatmap \
 --TFBS ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c1}}_bound.bed ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c1}}_unbound.bed  \
 --TFBS-labels "bound" "unbound" \
 --TFBS ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c2}}_bound.bed ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c2}}_unbound.bed \
 --TFBS-labels "bound" "unbound" \
 --signals ${{workdir}}/bias_correction/${{c1}}_corrected.bw ${{workdir}}/bias_correction/${{c2}}_corrected.bw \
 --output {output.heatmap} \
 --share_colorbar --sort_by -1 --signal_labels $c1 $c2 --title "${{tf}}"
TOBIAS PlotAggregate \
 --TFBS ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c1}}_bound.bed ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c2}}_bound.bed \
 --TFBS-labels "${{c1}}_bound" "${{c2}}_bound" \
 --signals ${{workdir}}/bias_correction/${{c1}}_corrected.bw ${{workdir}}/bias_correction/${{c2}}_corrected.bw \
 --signal_labels $c1 $c2 \
 --output {output.aggregate} \
 --share_y both --plot_boundaries --title "${{tf}}"
"""
