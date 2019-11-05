#!/bin/bash

mkdir -p data/fake
n=0
while [ $n -lt 1500 ]
do
    n=$(($n + 1))
    nm=$(echo ${n} | md5sum | cut -c 1-12)
    fn="data/fake/${nm}.tsv"
    s=$((17 + $n))
    # 150x coverage, 0.2% errors, 1% SNPs
    pypy fakeCounts.py some-genes.bed $s 150 0.002 0.01 ~/data/hg19 > ${fn}
done
