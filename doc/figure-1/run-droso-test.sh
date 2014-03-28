#!/bin/bash

file=/local/bigdata/SRR332538.fasta

cd ../../
make deb=1 
cd /local/bigdata/

for mem in 1 10 100 1000
do
    echo "mem=$mem"
    for disk_ratio in 0.1 0.25 0.5 1 2 4
    do
        read_size=`du $file -b | cut -f1`
        disk=`echo "$disk_ratio * $read_size / 1024 / 1024" | bc`
        time ~/dsk/paper/figure-1/monitor-memory.sh ~/dsk/dsk $file 21 -m $mem -d $disk -vv
    done
done
