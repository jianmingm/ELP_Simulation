#!/bin/bash

if [ "$#" -lt 3 ]; then
    echo "Error: Insufficient arguments provided."
    echo "Usage: $0 <m> <n> <save> <save_dir>"
    exit 1
fi

m=$1
n=$2
save=$3
save_dir=$4

## The residues
set1="G T S W Y H E Q D N K R"
set2="I V L F C M A"

for x in $set1; do
    for y in $set2; do
        python main.py -x1 $x -m $m -x2 $y -n $n -save_dir ${save_dir}/m${m}n${n} $save 
    done
done

