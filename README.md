## CCBR TOBIAS snakemake pipeline

TOBIAS or "Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal" is a framework of tools for investigating transcription factor binding from ATAC-seq signal. The analysis involves numerous sequential steps (or tasks) to be executed in order to successfully predict TF occupancy footprint from deduplicated alignment BAM files of ATACseq raw data (fastq files). Here we use Snakemake to automate the sequential execution on any HPC. Most tools used by the pipeline are completely containerized in docker format and can be invoked using singularity on the HPC. The minimum requirements for running this pipeline are:

* Python (>=3.5)
* Snakemake(>=5.24.1)
* Singularity(>=3.7.4)

This pipeline was built using [the CCBR_SnakemakePipelineCookiecutter](https://github.com/CCBR/CCBR_SnakemakePipelineCookiecutter).

Please visit the following pages for more details directly from the authors of TOBIAS:

* https://github.molgen.mpg.de/pages/loosolab/www/software/TOBIAS/
* https://github.com/loosolab/TOBIAS
* https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake

### Quick start instructions for running CCBR_tobias on [Biowulf](https://hpc.nih.gov/)

Various version of the pipeline have been checked out at `/data/CCBR_Pipeliner/Pipelines/CCBR_tobias` on biowulf. You can get help about running the pipeline using:

```bash
%  bash /data/CCBR_Pipeliner/Pipelines/CCBR_tobias/v0.2/run_tobias.bash --help
Pipeline Dir: /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/CCBR_tobias/v0.2
Git Commit/Tag: 6c8726023269ace0fd8fe886a1213859b363f9fd	v0.2
/data/CCBR_Pipeliner/Pipelines/CCBR_tobias/v0.2/run_tobias.bash: run CCBR TOBIAS workflow for ATAC seq data
USAGE:
  bash /data/CCBR_Pipeliner/Pipelines/CCBR_tobias/v0.2/run_tobias.bash -m/--runmode=<MODE> -w/--workdir=<path_to_workdir>
Required Arguments:
1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) run : run with slurm
    *) reset : DELETE workdir dir and re-init it
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
    *) runlocal : run without submitting to sbatch
2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
```

The pipeline requires only 2 arguments:

* Runmode
* Working dir

Generally, we anticipate CCBR_tobais to be run in 3 steps:

#### 1. Initialize

```bash
% bash /data/CCBR_Pipeliner/Pipelines/CCBR_tobias/dev/run_tobias.bash -m=init -w=/path/to/outfolder
```

This creates the output folder, so it should not exists before running `init`. Along with other scripts and files, `init` copies `config.yaml` and `cluster.json`  to the output folder, which can then be edited by the user. Some key input values that need to be edited before running the pipeline are as follows:

* `data`: points to the CCBR_ATACseq `dedup.bam` replicate files per sample. The sample names should match those later used in `contrasts`

* `contrasts`: which contrasts to perform using TOBIAS. The 2 groups should be already defined under `data`

* `peaks`: areas of interests to query for differential foot printing. This should be manually curated before running CCBR_tobias pipeline

* `genome`: currently supports **mm10** for mouse with Gencode M21 annotation and **hg38** for human Gencode v30 annotation.

* `motifs`: motif database to use for analysis. The choices are:

| database      | organism    | version            |
| ------------- | ----------- | ------------------ |
| HOCOMOCO\_v11 | Human       | Core               |
| HOCOMOCO\_v11 | Human       | Full               |
| HOCOMOCO\_v11 | Mouse       | Core               |
| HOCOMOCO\_v11 | Mouse       | Full               |
| HOCOMOCO\_v11 | Human+Mouse | Core               |
| HOCOMOCO\_v11 | Human+Mouse | Full               |
| JASPAR2020    | \-          | core\_nonredundant |
| JASPAR2020    | \-          | core\_redundant    |
| JASPAR2020    | vertebrate  | core\_nonredundant |
| JASPAR2020    | vertebrate  | core\_redundant    |

#### 2. Dryrun

```bash
% bash /data/CCBR_Pipeliner/Pipelines/CCBR_tobias/dev/run_tobias.bash -m=dryrun -w=/path/to/outfolder
```

Running the above command ensures that 

* output folder exists and contains the required files
* examples the `config.yaml` files and makes sure that we have appropriate permissions to the input files and output locations
* runs snakemake in `dry-run` mode using the `cluster.json` to enlist a table of rules/tasks to be run 

#### 3. Run

After successfully running `dryrun` , the user can run the same command with `-m=run` option to submit jobs to the slurm job scheduler on biowulf. By default, the `norm` partition is used to running jobs, but that and other job parameters can be changed by editing the `cluster.json` file in the output folder.

### Expected Outputs:

The following folders are expected upon successful completion.

#### bams

Individual replicate alignment BAMs are merged together and pre-sorted. This folder will contains the merged BAMs

#### coverage

The merged BAMs are converted to normalized bigwigs for visualization with IGV. The bigwigs can be found here.

#### bias_correction

The merged BAMs from the `bams` folder are corrected for Tn5 insertion bias. 4 separate bigwigs are expected as output on a per-condition basis:

* uncorrected bigwig: The uncorrected cutsite signal representing observed reads in basepair resolution. This track is normalized for sequencing depth but not corrected in terms of Tn5 bias.
* corrected bigwig: This is the corrected cutsite signal and will contain both positive and negative values for positions respectively more or less cut than expected. Remember, bigwigs cannot have positive and negative values at the same coordinate.
* biased bigwig: The raw bias score against the PWM/DWM bias matrix. This is purely based on sequence.
* expected bigwig: Knowing the cutsite preferences of the Tn5 enzyme the expected cutsite signal is reported here given the influence of bias. It is the raw bias score scaled towards the sum of cuts in the region, and can be directly compared to the uncorrected signal.
* pdf: Plot showing the observed Tn5 bias before and after correction can be found here.

#### footprinting

Using the bias corrected corrected bigwig a per-condition footprinting bigwig is created limited to the "regions of interest" defined by the `peaks` in the `config.yaml`.

#### peaks

Supplied `peaks` are annotated using *UROPA* and annotations are stored here.

#### TFBS_{contrast}

One TFBS folder is create for each contrast. There are created by running `bindetect`. Each TFBS folder contains numerous (100s) subfolders, one for each motif in the motif database selected using `motifs` parameter in `config.yaml`. Each of these per-TF-motif subfolder also has a standard folder structure including a subfolder name `beds`. This contains:

* a bed file for the TFBS for the motif in consideration which fall within the "regions of interest" as declared by the `peaks` parameter in `config.yaml`
* the TFBS sites in the above bed file are split into "bound" and "unbound" sites for each contrast separately resulting into a total of 4 bed files.

More more details see https://github.com/loosolab/TOBIAS/wiki/BINDetect

**<u>Caution</u>** This folder has a large digital footprint. Approximately, each contrast produces files amounting to about 40-60 GB. Hence, only run those contrasts that are interesting. DO NOT RUN ALL JUST BECAUSE YOU CAN!

This folder also contains:

* bindetect_results.txt
* bindetect_figures.pdf

which are the key results for this contrast as a table and as plots.

#### overview_{contrast}

All "bound" bed for all the TF motifs considered are concatenated together to be reported here as 2 sorted and indexed bed files. As these are indexed they can be easily loaded in a IGV session for visual inspection.





