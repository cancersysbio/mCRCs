#!/bin/bash

if [ $# -lt 2 ]
then
    echo "missing arguments: <file.bam> <ref_file.fai>"
    exit -1
fi

NAME=`echo $1 | sed 's/\.bam$//g'`
REF=$2

samtools sort -o $NAME.sorted.bam -O bam -T $NAME $NAME.bam
mv $NAME.sorted.bam $NAME.bam
samtools index $NAME.bam
