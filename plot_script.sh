#!/bin/bash

#srun -p az4-mixed --exclusive ./bin/murb -i 10 -n 1000 --im cpu+simd --gf --nv > results.txt

BINARY="./build/bin/murb"
OUTPUT_FILE="plot_results.csv"
ITERATIONS=1000

NB_BODIES=(1000 2000 5000 10000 15000 20000)

IMPLNAME=("cpu+naive" "cpu+optim" "cpu+simd" "cpu+barneshut" "cpu+barneshut+omp" "cpu+omp")
export OMP_NUM_THREADS=16
export OMP_SCHEDULE="static,16"

echo "implementation,iterations,nb_bodies,time,fps,gflops" > $OUTPUT_FILE

for impl in "${IMPLNAME[@]}"
do
	for n in "${NB_BODIES[@]}"
	do
		output=$(sbatch ./run.sh static 16 16 $n $ITERATIONS $impl)
		echo $BINARY -i $ITERATIONS -n $n -im $impl --gf --nv
		
		info_line=$(echo "$output" | grep "Entire simulation took")

		time=$(echo "$info_line" | awk '{print $4}')
		fps=$(echo "$info_line" | awk '{print $6}' | tr -d '(')
		gflops=$(echo "$info_line" | awk '{print $8}')
		echo "$impl,$ITERATIONS,$n,$time,$fps,$gflops" >> $OUTPUT_FILE

	done
	echo ""
done
