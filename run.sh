#!/bin/bash
#SBATCH --job-name=omp_runner
#SBATCH --output=./omprunner_%j.out
#SBATCH --error=./omprunner_%j.err
#SBATCH --nodes=1
#SBATCH --time=00:02:00
#SBATCH --partition=az4-n4090
#SBATCH --exclusive

if [ $# -lt 6 ]; then
        echo "!!!!ERROR!!!! usage : sbatch $0 <tipo> <chunk> <ncore> <n> <i> <impl>"
        exit -1
fi

echo "Executing [$SLURM_JOB_ID - $SLURM_JOB_NAME] in $SLURM_JOB_NODELIST"

tipo=$1
chunk=$2
ncore=$3
nbod=$4
its=$5
impl=$6

export OMP_NUM_THREADS=$ncore
export OMP_SCHEDULE="$tipo,$chunk"
echo "$typo, $chunk - $ncore"
./build/bin/murb -n $nbod -i $its --im $impl --nv --gf
