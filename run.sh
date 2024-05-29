
samplesheet="./assets/samplesheet_test.csv"
workdir="./work-short"
outdir="./results"
binddir="/projects/"

#----------------------------------------------------------------------------------------------------

SLURM_SUBMIT_DIR=`pwd`

#sbatch \
submit.sb $workflow "$samplesheet" "$workdir" "$outdir" "$binddir"