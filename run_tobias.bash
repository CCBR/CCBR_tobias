#!/usr/bin/env bash
# Author: Vishal Koparde, Ph.D.
# CCBR, NCI
# (c) 2021
#
# wrapper script to run the TOBIAS snakemake workflow
# run tobias
# https://github.com/loosolab/tobias/
# ## clone the pipeline to a folder
# ## git clone https://github.com/loosolab/TOBIAS.git

set -eo pipefail
module purge

SINGULARITY_BINDS="-B /data/CCBR_Pipeliner/:/data/CCBR_Pipeliner/"

function get_git_commitid_tag() {
  cd $1
  gid=$(git rev-parse HEAD)
  tag=$(git describe --tags $gid 2>/dev/null)
  echo -ne "$gid\t$tag"
}

# ## setting PIPELINE_HOME
PIPELINE_HOME=$(readlink -f $(dirname "$0"))
echo "Pipeline Dir: $PIPELINE_HOME"
# set snakefile
SNAKEFILE="${PIPELINE_HOME}/tobias.snakefile"

# get github commit tag
GIT_COMMIT_TAG=$(get_git_commitid_tag $PIPELINE_HOME)
echo "Git Commit/Tag: $GIT_COMMIT_TAG"

function usage() { cat << EOF
run_tobias.sh: run TOBIAS for ATAC seq data
USAGE:
  bash run_tobias.sh <MODE>  <path_to_workdir>
Required Positional Argument:
  MODE: [Type: Str] Valid options:
    a) init <path_to_workdir> : initialize workdir
    b) run <path_to_workdir>: run with slurm
    c) reset <path_to_workdir> : DELETE workdir dir and re-init it
    e) dryrun <path_to_workdir> : dry run snakemake to generate DAG
    f) unlock <path_to_workdir> : unlock workdir if locked by snakemake
    g) runlocal <path_to_workdir>: run without submitting to sbatch
EOF
}

function err() { cat <<< "
#
#
#
  $@
#
#
#
" && usage && exit 1 1>&2; }

function init() {

if [ "$#" -eq "1" ]; then err "init needs an absolute path to the working dir"; fi
if [ "$#" -gt "2" ]; then err "init takes only one more argument"; fi
WORKDIR=$2
x=$(echo $WORKDIR|awk '{print substr($1,1,1)}')
if [ "$x" != "/" ]; then err "working dir should be supplied as an absolute path"; fi
echo "Working Dir: $WORKDIR"
if [ -d $WORKDIR ];then err "Folder $WORKDIR already exists!"; exit 1; fi
mkdir -p $WORKDIR
sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $WORKDIR/config.yaml

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

echo "Done Initializing $WORKDIR. You can now edit $WORKDIR/config.yaml"

}

function runcheck(){
  if [ "$#" -eq "1" ]; then err "absolute path to the working dir needed"; usage; exit 1; fi
  if [ "$#" -gt "2" ]; then err "too many arguments"; usage; exit 1; fi
  WORKDIR=$2
  echo "Working Dir: $WORKDIR"
  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!"; exit 1; fi
  module load python/3.7
  module load snakemake/5.24.1
  SINGULARITY_BINDS="$SINGULARITY_BINDS -B ${PIPELINE_HOME}:${PIPELINE_HOME} -B ${WORKDIR}:${WORKDIR}"
}

function dryrun() {
  runcheck "$@"
  run "--dry-run"
}

function unlock() {
  runcheck "$@"
  run "--unlock"  
}

function runlocal() {
  runcheck "$@"
  if [ "$SLURM_JOB_ID" == "" ];then err "runlocal can only be done on an interactive node"; exit 1; fi
  module load singularity
  run "local"
}

function runslurm() {
  runcheck "$@"
  run "slurm"
}

function preruncleanup() {
  echo "Running..."

  cd $WORKDIR
  ## check if initialized
  for f in config.yaml; do
    if [ ! -f $WORKDIR/$f ]; then err "Error: '${f}' file not found in workdir ... initialize first!";usage && exit 1;fi
  done
  ## Archive previous run files
  if [ -f ${WORKDIR}/snakemake.log ];then 
    modtime=$(stat ${WORKDIR}/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/snakemake.log ${WORKDIR}/stats/snakemake.${modtime}.log
    if [ -f ${WORKDIR}/snakemake.log.HPC_summary.txt ];then 
      mv ${WORKDIR}/snakemake.log.HPC_summary.txt ${WORKDIR}/stats/snakemake.${modtime}.log.HPC_summary.txt
    fi
    if [ -f ${WORKDIR}/snakemake.stats ];then 
      mv ${WORKDIR}/snakemake.stats ${WORKDIR}/stats/snakemake.${modtime}.stats
    fi
  fi
  nslurmouts=$(find ${WORKDIR} -maxdepth 1 -name "slurm-*.out" |wc -l)
  if [ "$nslurmouts" != "0" ];then
    for f in $(ls ${WORKDIR}/slurm-*.out);do gzip -n $f;mv ${f}.gz ${WORKDIR}/logs/;done
  fi

}

function postrun() {
  bash ${PIPELINE_HOME}/scripts/gather_cluster_stats.sh ${WORKDIR}/snakemake.log > ${WORKDIR}/snakemake.log.HPC_summary.txt
}

function run() {


  if [ "$1" == "local" ];then

  preruncleanup

  snakemake -s ${PIPELINE_HOME}/$SNAKEFILE \
  --directory $WORKDIR \
  --printshellcmds \
  --use-singularity \
  --singularity-args "$SINGULARITY_BINDS" \
  --use-envmodules \
  --latency-wait 120 \
  --configfile ${WORKDIR}/config.yaml \
  --cores all \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  if [ "$?" -eq "0" ];then
    snakemake -s ${PIPELINE_HOME}/$SNAKEFILE \
    --report ${WORKDIR}/runlocal_snakemake_report.html \
    --directory $WORKDIR \
    --configfile ${WORKDIR}/config.yaml 
  fi

  postrun

  elif [ "$1" == "slurm" ];then
  
  preruncleanup

  cat > ${WORKDIR}/submit_script.sbatch << EOF
#!/bin/bash
#SBATCH --job-name="tobias"
#SBATCH --mem=10g
#SBATCH --partition="ccr,norm"
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2

module load python/3.7
module load snakemake/5.24.1
module load singularity

cd \$SLURM_SUBMIT_DIR

snakemake -s $SNAKEFILE \
--directory $WORKDIR \
--use-singularity \
--singularity-args "$SINGULARITY_BINDS" \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile ${WORKDIR}/config.yaml \
--cluster-config ${PIPELINE_HOME}/config/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} --qos=ccrcore" \
-j 500 \
--rerun-incomplete \
--keep-going \
--stats ${WORKDIR}/snakemake.stats \
2>&1|tee ${WORKDIR}/snakemake.log

if [ "\$?" -eq "0" ];then
  snakemake -s $SNAKEFILE \
  --directory $WORKDIR \
  --report ${WORKDIR}/runslurm_snakemake_report.html \
  --configfile ${WORKDIR}/config.yaml 
fi

bash ${PIPELINE_HOME}/scripts/gather_cluster_stats.sh ${WORKDIR}/snakemake.log > ${WORKDIR}/snakemake.log.HPC_summary.txt

EOF

  sbatch ${WORKDIR}/submit_script.sbatch

  else

snakemake $1 -s $SNAKEFILE \
--directory $WORKDIR \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile ${WORKDIR}/config.yaml \
--cluster-config ${PIPELINE_HOME}/config/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
-j 500 \
--rerun-incomplete \
--keep-going \
--stats ${WORKDIR}/snakemake.stats

  fi

}

function reset() {
if [ "$#" -eq "1" ]; then err "cleanup needs an absolute path to the existing working dir"; usage; fi
if [ "$#" -gt "2" ]; then err "cleanup takes only one more argument"; usage; fi
WORKDIR=$2
echo "Working Dir: $WORKDIR"
if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!";fi
echo "Deleting $WORKDIR"
rm -rf $WORKDIR
echo "Re-Initializing $WORKDIR"
init "$@"
}


function main(){

  if [ $# -eq 0 ]; then usage; exit 1; fi

  case $1 in
    init) init "$@" && exit 0;;
    dryrun) dryrun "$@" && exit 0;;
    unlock) unlock "$@" && exit 0;;
    run) runslurm "$@" && exit 0;;
    runlocal) runlocal "$@" && exit 0;;
    reset) reset "$@" && exit 0;;
    -h | --help | help) usage && exit 0;;
    -* | --*) err "Error: Failed to provide mode: <init|run>."; usage && exit 1;;
    *) err "Error: Failed to provide mode: <init|run>. '${1}' is not supported."; usage && exit 1;;
  esac
}

main "$@"





