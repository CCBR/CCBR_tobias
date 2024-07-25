#######################################################
from os.path import join
import sys
import os
import pandas as pd
import yaml
import glob
#######################################################

#######################################################
# INCLUDE OTHER .smk FILES
#######################################################
include: 'rules/init.smk'
#######################################################

#######################################################
# INPUT FUNCTIONS
#######################################################

def get_replicate_bams_for_condition(wildcards):
    """
    Return a list of bams for this condition
    """
    return CONDITION2BAMS[wildcards.condition]

#######################################################

def get_condition_footprint_bigwigs_from_contrast(wildcards):
    files=[]
    for c in CONTRASTS2CONDITIONS[wildcards.contrast]:
        files.append(join(WORKDIR, "footprinting", c+"_footprints.bw"))
    return files

#######################################################

def get_bound_bed_list(wildcards) :
    files=dict()
    for cont,conds in CONTRASTS2CONDITIONS.items():
        for cond in conds:
            files[cont+"."+cond]=join(WORKDIR, "overview_"+cont, "all_"+cond+"_bound.bed.gz")
    return files

#######################################################

def get_bound_bed_list_per_contrast_per_condition_per_TF(wildcards) :
    files=dict()
    for cont,conds in CONTRASTS2CONDITIONS.items():
        for cond in conds:
            for tf in TFs:
                files[cont+"."+cond+"."+TF]=join(WORKDIR, "TFBS_"+cont, tf, "beds", tf+"_"+cond+"_bound.bed")
    return files

#######################################################

#######################################################
# LOCAL RULES
# These will not be submitted to the job scheduler but run locally
#######################################################

# localrules: all, plot_heatmaps_aggregates
localrules: all, plot_heatmaps_aggregates_init, reformat_annotated_bed, map_motifs_to_genes

#######################################################

#######################################################
# RULE ALL
#######################################################
rule all:
    input:
        # merge replicate BAMs to per-condition BAM
        expand(join(WORKDIR, "bams", "{condition}.bam"), condition=CONDITION_IDS),
        # create per-condition coverage BIGWIG
        expand(join(WORKDIR, "coverage", "{condition}_coverage.bw"), condition=CONDITION_IDS),
        # correct for Tn5 bias
        expand(join(WORKDIR, "bias_correction", "{condition}_uncorrected.bw"), condition=CONDITION_IDS),
        expand(join(WORKDIR, "bias_correction", "{condition}_bias.bw"), condition=CONDITION_IDS),
        expand(join(WORKDIR, "bias_correction", "{condition}_expected.bw"), condition=CONDITION_IDS),
        expand(join(WORKDIR, "bias_correction", "{condition}_corrected.bw"), condition=CONDITION_IDS),
        # do footprinting
        expand(join(WORKDIR, "footprinting", "{condition}_footprints.bw"), condition=CONDITION_IDS),
        # run Bindetect for each contrast
        expand(join(WORKDIR, "TFBS_{contrast}", "bindetect_figures.pdf"),contrast=CONTRASTS),
        expand(join(WORKDIR,"TFBS_{contrast}","bound_beds_list.tsv"),contrast=CONTRASTS),
        expand(join(WORKDIR, "overview_{cont}", "all_{cond}_bound.bed.gz"),zip,cont=CC1,cond=CC2),
        # plots
        expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.bash"),contrast=CONTRASTS,TF=TFs),
        expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.heatmap.pdf"),contrast=CONTRASTS,TF=TFs),
        expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.aggregate.pdf"),contrast=CONTRASTS,TF=TFs),
        # network
        # only needed for first condition in each contrast
        [f"{WORKDIR}/network/{contrast}/{conditions_list[0]}/edges.txt" for contrast, conditions_list in CONTRASTS2CONDITIONS.items()],
#######################################################

#######################################################

rule tobias_download_data:
    """ 
    Download example data from Tobias
    """
    output:
        directory('data-tobias-2020')
    container: CONTAINERS["tobias"]
    shell: """
        TOBIAS DownloadData --bucket data-tobias-2020
        """


rule condition_bam:
# """ 
# Merge all replicate BAMs into a single "condition" BAM file
# lambda wildcards: config["data"][wildcards.condition]
# """
    input:
        get_replicate_bams_for_condition
    output:
        bam = join(WORKDIR,"bams","{condition}.bam"),
    envmodules: TOOLS["sambamba"]["version"]
    threads: getthreads("condition_bam")
    message:
        "Running {rule} for condition: {wildcards.condition}"
    params:
        outdir=join(WORKDIR,"bams"),
        mem=getmemG("condition_bam"),
        condition="{condition}",
    shell: """
        if [ -w "/lscratch/${{SLURM_JOB_ID}}" ]; then tmpdir="/lscratch/${{SLURM_JOB_ID}}"; else tmpdir="/dev/shm"; fi
        sorted_bams=""
        for bam in {input}
        do
        bn=$(basename $bam)
        sambamba sort --memory-limit={params.mem} --tmpdir=${{tmpdir}} --nthreads={threads} --out=/lscratch/${{SLURM_JOBID}}/${{bn%.*}}.sorted.bam $bam
        sorted_bams="$sorted_bams ${{tmpdir}}/${{bn%.*}}.sorted.bam"
        done
        sambamba merge --nthreads={threads} {output.bam} $sorted_bams
        rm -rf ${{tmpdir}}/*
        """

#######################################################

rule coverage_bw:
# """
# Convert the merged (per condition) BAM file to a unnormalized BIGWIG file
# """
    input:
        bam = rules.condition_bam.output.bam
    output:
        coveragebw = join(WORKDIR, "coverage", "{condition}_coverage.bw")
    params:
        condition = "{condition}"
    message:
        "Running {rule} for condition: {wildcards.condition} ({input.bam})"
    threads: getthreads("coverage_bw")
    envmodules: TOOLS["bedtools"]["version"],TOOLS["samtools"]["version"], TOOLS["ucsc"]["version"]
    shell: """
        if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then tmpdir="/lscratch/${{SLURM_JOB_ID}}";else tmpdir="/dev/shm";fi
        bedtools genomecov -ibam {input.bam} -bg >  ${{tmpdir}}/{params.condition}.bg
        bedSort ${{tmpdir}}/{params.condition}.bg ${{tmpdir}}/{params.condition}.bg
        samtools view -H {input.bam}|grep '^@SQ'|cut -f2,3|sed 's/SN:\|LN://g' > ${{tmpdir}}/{params.condition}.chrom.sizes
        bedGraphToBigWig ${{tmpdir}}/{params.condition}.bg ${{tmpdir}}/{params.condition}.chrom.sizes {output.coveragebw}
        rm -rf ${{tmpdir}}/*
        """

#######################################################

rule atacorrect:
# """
# Correcting the per-condition merged BAM file for Tn5 Bias
# Output 4 BIGWIG files per-condition: uncorrected, corrected, bias, expected
# """
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
    threads: getthreads("atacorrect")
    message: 
        "Running {rule} for condition: {wildcards.condition} ({input.bam})"
    container: CONTAINERS["tobias"]
    shell:  """
        TOBIAS ATACorrect \
            -b {input.bam} \
            -g {input.genome} \
            -p {input.peaks} \
            --cores {threads} \
            --prefix {params.condition} \
            --outdir {params.outdir} \
            {params.blacklist} {params.autocorrect_extra_params}
        """

#######################################################

rule footprinting:
# """
# Create a per-condition footprinting BIGWIG from a per-condition bias-corrected BIGWIG
# """
    input: 
        signal = rules.atacorrect.output.corrected, 	#os.path.join(OUTPUTDIR, "bias_correction", "{condition}_corrected.bw"),
        regions = PEAKS,		#os.path.join(OUTPUTDIR, "peak_calling", "all_merged.bed")
        annotatedpeaks = join(WORKDIR,"peaks","annotated_peaks"),   # annotated_peaks are not required by this rule, but this makes sure that uropa is run upstream from footprinting, which is in turn upstream from bindetect which requires annotated peaks
        annotatedpeaksheader = join(WORKDIR,"peaks","annotated_peaks_header")
    output: 
        footprints = join(WORKDIR, "footprinting", "{condition}_footprints.bw"),
    params:
        footprinting_extra_params = FOOTPRINTING_EXTRA_PARAMS
    threads: getthreads("footprinting")
    message: 
        "Running {rule} for condition: {wildcards.condition} ({input.signal})"
    container: CONTAINERS["tobias"]
    shell: """
        TOBIAS FootprintScores \
            --signal {input.signal} \
            --regions {input.regions} \
            --output {output.footprints} \
            --cores {threads} {params.footprinting_extra_params}
        """

#######################################################

rule uropa:
# """
# Annotate peaks with Uropa and create the header file required later
# """
    input:
        peaks = PEAKS		# BED format
    output:
        annotatedpeaks = join(WORKDIR,"peaks","annotated_peaks"),
        annotatedpeaksheader = join(WORKDIR,"peaks","annotated_peaks_header")
    params:
        outdir = join(WORKDIR,"peaks"),
        uropa_base_config = UROPABASEYAML,
        gtf = GTF,
        create_uropa_json_script = join(SCRIPTSDIR,"create_uropa_config.py")
    threads: getthreads("uropa")
    envmodules: TOOLS["uropa"]["version"], TOOLS["python"]["version"]
    message: 
        "Running {rule}:"
    shell: """
        cd {params.outdir}
        cp {input.peaks} .
        python {params.create_uropa_json_script} \
            -j {params.uropa_base_config} \
            -b {input.peaks} \
            -g {params.gtf} > uropa.json
        uropa -i uropa.json -t {threads}
        head -n1 peaks_finalhits.txt > {output.annotatedpeaksheader}
        tail -n +2 peaks_finalhits.txt > {output.annotatedpeaks}
        """

#######################################################

# lambda wildcards: unpack(c2c(wildcards))
rule bindetect:
    input:
        signals = get_condition_footprint_bigwigs_from_contrast,
        peaks = rules.uropa.output.annotatedpeaks,
        peak_header = rules.uropa.output.annotatedpeaksheader,	
    output:
        pdf=join(WORKDIR,"TFBS_{contrast}","bindetect_figures.pdf"),
    params:
        motifs = MOTIFS, 		
        genome = REFFA,
        bindetect_extra_params=BINDETECT_EXTRA_PARAMS,
    threads: getthreads("bindetect")
    container: CONTAINERS["tobias"]
    message: 
        "Running {rule} for {wildcards.contrast}"
    shell: """
        c1=$(echo {wildcards.contrast} | awk -F\"_vs_\" '{{print $1}}')
        c2=$(echo {wildcards.contrast} | awk -F\"_vs_\" '{{print $2}}')
        outdir=$(dirname "{output.pdf}")
        TOBIAS BINDetect \
            --motifs {params.motifs} \
            --signals {input.signals} \
            --genome {params.genome} \
            --peaks {input.peaks} \
            --peak-header {input.peak_header} \
            --cores {threads} \
            --cond-names $c1 $c2 \
            --outdir $outdir \
            {params.bindetect_extra_params}
        sleep 600
        """

rule reformat_annotated_bed:
    """
    Strip ensembl gene version numbers, shorten motif IDs, and combine into one bed file.
    Annotated bed files from uropa contain version numbers,
    but gene names in map_motifs_to_genes do not.
    """
    input:
        bindetect=rules.bindetect.output.pdf, # actually using annotated TFBS files, but checkpoints are cumbersome to deal with
        py=f"{SCRIPTSDIR}/reformat_annotated_bed.py"
    output:
        bed=f"{WORKDIR}/TFBS_{{contrast}}_reformatted/{{condition}}.bed",
        tmp=temp(f"{WORKDIR}/TFBS_{{contrast}}_reformatted/{{condition}}.bed.tmp")
    container: CONTAINERS["base"]
    shell: """
        bed_files=$(dirname {input.bindetect})/*/beds/*{wildcards.condition}_bound.bed
        cat $bed_files > {output.tmp}
        python {SCRIPTSDIR}/reformat_annotated_bed.py {output.tmp} {output.bed}
        """


rule map_motifs_to_genes:
    input:
        gtf=GTF,
        pfm=config["pfm"],
        py=f"{SCRIPTSDIR}/motif2gene_mapping.py"
    output:
        txt=f"{WORKDIR}/motif2gene/{GENOME}.motif2gene_mapping.txt"
    container: CONTAINERS["base"]
    script:
        "scripts/motif2gene_mapping.py"

rule create_network:
    input:
        origin=rules.map_motifs_to_genes.output.txt,
        bed=rules.reformat_annotated_bed.output.bed
    output:
        adjacency=f"{WORKDIR}/network/{{contrast}}/{{condition}}/adjacency.txt",
        edges=f"{WORKDIR}/network/{{contrast}}/{{condition}}/edges.txt"
    container: CONTAINERS["tobias"]
    shell: """
        TOBIAS CreateNetwork \
            --TFBS {input.bed} \
            --origin {input.origin} \
            --outdir $(dirname {output.edges})
        """

#######################################################

rule create_cont_cond_tf_join_filelist:
    input: 
        rules.bindetect.output.pdf
    output: 
        tsv=temp(join(WORKDIR,"TFBS_{contrast}","bound_beds_list.tsv.tmp")),
    run:
        cont=wildcards.contrast
        conds=CONTRASTS2CONDITIONS[cont]
        out=open(output[0], "w")
        for cond in conds:
            for tf in TFs:
                f=join(WORKDIR, "TFBS_"+cont, tf, "beds", tf+"_"+cond+"_bound.bed.gz")
                out.write("%s\t%s\n"%(cond,f))
        out.close()

#######################################################

rule bgzip_beds:
    input:
        rules.create_cont_cond_tf_join_filelist.output.tsv,
        [f"{WORKDIR}/network/{contrast}/{conditions_list[0]}/edges.txt" for contrast, conditions_list in CONTRASTS2CONDITIONS.items()] # run create_network first
    output:
        tsv=join(WORKDIR,"TFBS_{contrast}","bound_beds_list.tsv")
    params:
        workdir = WORKDIR,
        script = join(SCRIPTSDIR,"_bgzip_all_beds.bash")
    threads: getthreads("gzip_beds")
    envmodules: TOOLS["ucsc"]["version"], TOOLS["samtools"]["version"], TOOLS["parallel"]["version"]
    shell: """
        cd {params.workdir}
        awk '{{print $2}}' {input} | xargs dirname | sort | uniq | xargs -I % echo bash {params.script} % > do_bgzip
        parallel -j {threads} < do_bgzip
        # rm -f do_bgzip
        cp {input} {output}
        """

#######################################################

# Join bound estimates per condition
rule join_bound:
    input:
        rules.create_cont_cond_tf_join_filelist.output.tsv
    output:
        bedgz = join(WORKDIR, "overview_{contrast}", "all_{cond}_bound.bed.gz")
    threads: getthreads("join_bound")
    envmodules: TOOLS["ucsc"]["version"], TOOLS["samtools"]["version"]
    params:
        cont="{contrast}",
        cond="{cond}"
    shell: """
        x="{output.bedgz}"
        xbed="${{x%.*}}"
        while read a b;do 
            if [ "$a" == "{params.cond}" ];then
                zcat $b
            fi
        done < {input} > $xbed
        bedSort $xbed $xbed
        bgzip -f --threads {threads} $xbed
        tabix -f -p bed $x
        """

#######################################################

rule plot_heatmaps_aggregates_init:
    input:
        rules.bindetect.output.pdf,
        rules.bgzip_beds.output.tsv
    output:
        plotcmd = temp(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.bash"))
    params:
        contrast = "{contrast}",
        tf = "{TF}",
        workdir = WORKDIR,
        heatmap = join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.heatmap.pdf"),
        aggregate = join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.aggregate.pdf")
    # message: 
    #     "Running {rule} for Contrast:{wildcards.contrast} and TF:{wildcards.TF}"
    shell: """
        tmpdir="/dev/shm"
        chars="abcdefghijklmnopqrstuvwxyz"
        randomtxt=$(for i in 1 2 3 4 5 6 7 8; do echo -n "${{chars:RANDOM%${{#chars}}:1}}";done)
        randomtxt=$(echo "$randomtxt$RANDOM")
        mkdir ${{tmpdir}}/${{randomtxt}}
        tmpdir="${{tmpdir}}/${{randomtxt}}"

        contrast="{params.contrast}"
        TF="{params.tf}"
        workdir="{params.workdir}"
        c1=$(echo {params.contrast}|awk -F\"_vs_\" '{{print $1}}')
        c2=$(echo {params.contrast}|awk -F\"_vs_\" '{{print $2}}')
        x=$(echo "${{TF}}"|awk -F"_" '{{print NF/2}}')
        tf=$(echo "${{TF}}"|sed "s/_/\\t/g"|cut -f1,$x|sed "s/\\t/_/g")
        cat > {output.plotcmd} << EOF

        if [ ! -w "$tmpdir" ];then mkdir -p $tmpdir;fi

        zcat ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c1}}_bound.bed.gz > ${{tmpdir}}/${{TF}}_${{c1}}_bound.bed
        zcat ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c1}}_unbound.bed.gz > ${{tmpdir}}/${{TF}}_${{c1}}_unbound.bed
        zcat ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c2}}_bound.bed.gz > ${{tmpdir}}/${{TF}}_${{c2}}_bound.bed
        zcat ${{workdir}}/TFBS_${{contrast}}/${{TF}}/beds/${{TF}}_${{c2}}_unbound.bed.gz > ${{tmpdir}}/${{TF}}_${{c2}}_unbound.bed

        TOBIAS PlotHeatmap \
        --TFBS ${{tmpdir}}/${{TF}}_${{c1}}_bound.bed ${{tmpdir}}/${{TF}}_${{c1}}_unbound.bed  \
        --TFBS-labels "bound" "unbound" \
        --TFBS ${{tmpdir}}/${{TF}}_${{c2}}_bound.bed ${{tmpdir}}/${{TF}}_${{c2}}_unbound.bed \
        --TFBS-labels "bound" "unbound" \
        --signals ${{workdir}}/bias_correction/${{c1}}_corrected.bw ${{workdir}}/bias_correction/${{c2}}_corrected.bw \
        --output {params.heatmap} \
        --share_colorbar --sort_by -1 --signal_labels $c1 $c2 --title "${{tf}}"
        TOBIAS PlotAggregate \
        --TFBS ${{tmpdir}}/${{TF}}_${{c1}}_bound.bed ${{tmpdir}}/${{TF}}_${{c2}}_bound.bed \
        --TFBS-labels "${{c1}}_bound" "${{c2}}_bound" \
        --signals ${{workdir}}/bias_correction/${{c1}}_corrected.bw ${{workdir}}/bias_correction/${{c2}}_corrected.bw \
        --signal_labels $c1 $c2 \
        --output {params.aggregate} \
        --share_y both --plot_boundaries --title "${{tf}}"

        rm -rf ${{tmpdir}}

        EOF
        """

#######################################################

rule plot_heatmaps_aggregates:
    input:
        expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.bash"),contrast=CONTRASTS,TF=TFs),
    output:
        expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.heatmap.pdf"),contrast=CONTRASTS,TF=TFs),
        expand(join(WORKDIR,"TFBS_{contrast}","{TF}","plots","{TF}_{contrast}.aggregate.pdf"),contrast=CONTRASTS,TF=TFs),
    params:
        workdir = WORKDIR
    container: CONTAINERS["tobias"]
    threads: getthreads("plot_heatmaps_aggregates")
    message: 
        "Generating plots"
    shell: """
        set -exo pipefail

        for i in {input};
        do 
        echo "bash $i"
        done > {params.workdir}/do_plots
        #parallel -j {threads} < {params.workdir}/do_plots
        bash {params.workdir}/do_plots

        # cleanup
        # while read a b;do rm -f $b;done < {params.workdir}/do_plots
        rm -f {params.workdir}/do_plots
        """
