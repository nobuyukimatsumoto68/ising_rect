#!/bin/bash -l

module load gcc/13.2.0

app='wolff.o'
# app='wolff2.o'
make ${app}

ns=(4 6 8 10) # 1e5, omp 4
len=${#ns[@]}
for (( i=0; i<$len; i++ ));
do
    n=${ns[$i]}
    dir=mult${n}
    mkdir -p ${dir}

    cp run.sh ${dir}
    cp ${app} ${dir}
    cd ${dir}
    qsub -N ${app}-sq-${n} -v mult=$n -v binsize=100000 -v app=$app -j y run.sh
    echo $n
    cd ..
done


ns=(12 14 16) # 1e4, omp 8
len=${#ns[@]}
for (( i=0; i<$len; i++ ));
do
    n=${ns[$i]}
    dir=mult${n}
    mkdir -p ${dir}

    cp run.sh ${dir}
    cp ${app} ${dir}
    cd ${dir}
    qsub -N ${app}-sq-${n} -v mult=$n -v binsize=10000 -v app=$app -j y run.sh
    echo $n
    cd ..
done
