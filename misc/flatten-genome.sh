#!/bin/bash

out_dir=$1
shift;

mkdir -p ${out_dir}
for f in $@
do
    b=$(basename $f .fa.gz)
    g=${out_dir}/${b}.txt
    gzcat ${f} | tail -n +2 | tr -d ' \r\n\t\v' > ${g}
    z=$(wc -c < ${g}| tr -d ' \t')
    echo "$b	$z	$b.txt" >> ${out_dir}/toc.txt
done
