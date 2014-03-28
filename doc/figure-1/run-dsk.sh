#!/bin/bash

cd ../../
make deb=1 
cd /local/bigdata/

for mem in 10 50 100 500 1000 2000
do
    echo "mem=$mem"
    for disk_ratio in 0.1 0.5 1 2 4 6
    do
        read_size=`du SRR001665.fasta -b | cut -f1`
        disk=`echo "$disk_ratio * $read_size / 1024 / 1024" | bc`
        time ~/dsk/paper/figure-1/monitor-memory.sh ~/dsk/dsk SRR001665.fasta 21 -m $mem -d $disk -vv
    done
done
