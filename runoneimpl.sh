#!/bin/bash

if [ $# -lt 1 ]; then
        echo "!!!!fuck you!!!! usage : bash <impl>"
        exit -1
fi

BINARY="./build/bin/murb"
OUTPUT_FILE="plot_results.csv"
ITERATIONS=1000

NB_BODIES=(1000 2000 5000 10000 15000 20000)

export OMP_NUM_THREADS=16
export OMP_SCHEDULE="static,16"

impl=$1

for n in "${NB_BODIES[@]}"
do
        sbatch ./run.sh static 16 16 $n $ITERATIONS $impl
        echo $BINARY -i $ITERATIONS -n $n -im $impl
done
echo ""
