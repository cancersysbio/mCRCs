#!/bin/bash

if [ $# -lt 2 ]
then
    echo "missing arguments: <file.bam> <reference.fa>"
    exit -1
fi

SAMPLE=`basename $1 .bam`

INDEL=$mCRC_DIR/ref_genomes/indels.vcf
DBSNP=$mCRC_DIR/ref_genomes/dbsnp_132_b37.leftAligned.vcf

java -Xmx3g -jar $PICARD/MarkDuplicates.jar I=$1 O=$SAMPLE.rmdup.bam M=$SAMPLE.rmdup.metric VALIDATION_STRINGENCY=LENIENT
samtools index $SAMPLE.rmdup.bam


java -Xmx3g -jar $GATKROOT/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $2 -I $SAMPLE.rmdup.bam -known $INDEL -o $SAMPLE.target_intervals.list

java -Xmx3g -jar $GATKROOT/GenomeAnalysisTK.jar -T IndelRealigner -R $2 -I $SAMPLE.rmdup.bam -targetIntervals $SAMPLE.target_intervals.list -known $INDEL -o $SAMPLE.realigned.bam

rm $SAMPLE.rmdup.bam
rm $SAMPLE.rmdup*bai


java -Xmx3g -jar $GATKROOT/GenomeAnalysisTK.jar -T BaseRecalibrator -R $2 -I $SAMPLE.realigned.bam -knownSites $DBSNP -knownSites $INDEL -o $SAMPLE.recal_data.table 

java -Xmx3g -jar $GATKROOT/GenomeAnalysisTK.jar -T PrintReads -R $2 -I $SAMPLE.realigned.bam -BQSR $SAMPLE.recal_data.table --emit_original_quals -o $SAMPLE.recal.bam

rm $SAMPLE.realigned.bam
rm $SAMPLE.realigned*bai
