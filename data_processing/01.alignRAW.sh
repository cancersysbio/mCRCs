#!/bin/bash

if [ $# -lt 3 ]
then
    echo "missing arguments 1.fq 2.fq reference.fa readgroup"
    exit -1
fi

SAMPLE=`basename $1 .fastq.gz`
READGROUP="@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA"

if [ $# -eq 4 ]
then
    READGROUP=$4
fi


bwa mem -R $READGROUP $3 $1 $2 > $SAMPLE.sam
samtools view -bS $SAMPLE.sam > $SAMPLE.bam
rm $SAMPLE.sam
