for f in $(ls *.out); do
    OUTPUT_FILE=$f-parsed.txt
    info_line=$(cat $f | grep "Entire simulation took")
    impl=$(cat $f | grep implementation | awk '{print $5}')
    ITERATIONS=$(cat $f | grep iterations | awk '{print $7}')
    n=$(cat $f | grep bodies | awk '{print $7}')
    time=$(echo "$info_line" | awk '{print $4}')
    fps=$(echo "$info_line" | awk '{print $6}' | tr -d '(')
    gflops=$(echo "$info_line" | awk '{print $8}')
    echo "$impl,$ITERATIONS,$n,$time,$fps,$gflops" >> $OUTPUT_FILE
done

cat *-parsed.txt > data.csv