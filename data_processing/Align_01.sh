#!/bin/bash

module load samtools/1.2
module load bwa/0.7.10
module load java/7u03
module load gatk/3.4.0
module load picard-tools/1.111

if [ $# -lt 4 ]
then
    echo "missing arguments: <fastq folder> <fastq name> <bam name> <bed file>"
    exit -1
fi

REFDIR=$mCRC_DIR/ref_genomes
REF='broad.Homo_sapiens_assembly19.fa'
SCRIPTDIR=$mCRC_DIR/scripts/data_processing
DATADIR=$1
FASTQ=$2
SAMPLENAME=$3
BED=$4

echo $FASTQ 
echo $SAMPLENAME 

for i in $DATADIR/${FASTQ}_*R1_001.fastq.gz
do
    PAIR=`echo $i | sed s/R1/R2/`
    echo $i $PAIR
    SM=$SAMPLENAME
    ID=`basename $i .fastq.gz`
    
    echo "Align"
    $SCRIPTDIR/01.alignRAW.sh $i $PAIR $REFDIR/$REF "@RG\tID:$ID\tSM:$SM\tPL:ILLUMINA\tLB:$SM"
    SAMPLE=`basename $i _R1_001.fastq.gz`
    mv ${SAMPLE}_R1_001.bam $SAMPLE.bam

    echo "Sort bam files"
    $SCRIPTDIR/02.sam2sortedbam.sh $SAMPLE.bam $REFDIR/$REF.fai

    echo "Realignment"
    $SCRIPTDIR/03.local_realign.sh $SAMPLE.bam $REFDIR/$REF

    rm $SAMPLE.bam
    samtools view -b -x BD -x BI -o $SAMPLE.bam $SAMPLE.recal.bam
    rm $SAMPLE.recal.bam
    rm $SAMPLE.recal.bai
    rm $SAMPLE.recal_data.table
    rm $SAMPLE.target_intervals.list
    samtools index $SAMPLE.bam
done

echo "Merge"
ls -1 ${FASTQ}_*bam > $SAMPLENAME.list
samtools merge -b $SAMPLENAME.list $SAMPLENAME.bam
samtools index $SAMPLENAME.bam
rm ${FASTQ}_*bam*

SAMPLE=$SAMPLENAME

echo "Realignment"
$SCRIPTDIR/03b.local_realign_norecal.sh $SAMPLE.bam $REFDIR/$REF
rm $SAMPLE.bam
rm $SAMPLE.rmdup.bam
rm $SAMPLE.rmdup*bai
mv $SAMPLE.realigned.bam $SAMPLE.bam
rm $SAMPLE.realigned*bai
rm $SAMPLE.target_intervals.list
samtools index $SAMPLE.bam
