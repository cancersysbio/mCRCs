#!/bin/bash

##Redo alignment for bam files generated at Broad

module load samtools/1.3.1
module load bwa/0.7.10
module load java/7u03
module load bedtools/2.23.0
module load gatk/3.4.0

if [ $# -lt 2 ]
then
    echo "missing arguments: <old bam file> <bed file>"
    exit -1
fi

SCRIPTDIR=$mCRC_DIR/scripts/data_processing
REFDIR=$mCRC_DIR/ref_genomes
REF='broad.Homo_sapiens_assembly19.fa'
BAMFILE=$1
BAM=`basename $BAMFILE .bam`
BED=$2

echo $BAM

mkdir -p $BAM
cd $BAM

samtools split $BAMFILE

for i in ${BAM}_[0-9]*.bam
do
    echo $i
    NAME=`basename $i .bam`
    RG=`samtools view -H $i | grep ^@RG`
    RG2=`echo $RG | sed s/\ /'\\\t'/g`
    echo $RG2
    samtools sort -n -o $NAME.sorted.bam $i 
    rm $i
    samtools fastq -n -O -1 $NAME.1.fastq -2 $NAME.2.fastq $NAME.sorted.bam
    rm $NAME.sorted.bam
    gzip $NAME.1.fastq
    gzip $NAME.2.fastq

    echo "Align"
    $SCRIPTDIR/01.alignRAW.sh $NAME.1.fastq.gz $NAME.2.fastq.gz $REFDIR/$REF "$RG2"
    mv $NAME.1.bam $NAME.ra.bam
    rm $NAME.1.fastq.gz
    rm $NAME.2.fastq.gz

    echo "Sort bam files"
    $SCRIPTDIR/02.sam2sortedbam.sh $NAME.ra.bam $REFDIR/$REF.fai

    echo "Realignment"
    $SCRIPTDIR/03.local_realign.sh $NAME.ra.bam $REFDIR/$REF

    rm $NAME.ra.bam
    rm $NAME.ra.target_intervals.list
    samtools view -b -x BD -x BI -o $NAME.ra.bam $NAME.ra.recal.bam
    rm $NAME.ra.recal.bam
    rm $NAME.ra.recal.bai
    rm $NAME.ra.recal_data.table
    samtools index $NAME.ra.bam
done

echo "Merge"
ls -1 ${BAM}_*.ra.bam > $BAM.list
samtools merge -b $BAM.list ${BAM}_RA.bam
samtools index ${BAM}_RA.bam
rm ${BAM}_*.ra.bam
rm ${BAM}_*.ra.bam.bai

echo "Realignment"
$SCRIPTDIR/03b.local_realign_norecal.sh ${BAM}_RA.bam $REFDIR/$REF
rm ${BAM}_RA.bam
rm ${BAM}_RA.bam.bai
rm ${BAM}_RA.rmdup.bam
rm ${BAM}_RA.rmdup*bai
rm ${BAM}_RA.target_intervals.list
rm ${BAM}_RA.realigned*bai
mv ${BAM}_RA.realigned.bam ${BAM}_RA.bam
samtools index ${BAM}_RA.bam
