#!/bin/bash

if [ $# -lt 2 ]
then
    echo "missing arguments: <patient name> <tumorbam.list>"
    exit -1
fi

REFDIR=$mCRC_DIR/ref_genomes
REF=$REFDIR/broad.Homo_sapiens_assembly19.fa

PATIENT=$1
SNV=$PATIENT.combined.varscanfilter

cat $SNV | awk -F"\t" '{print $1"\t"$2-1"\t"$2}' > $PATIENT.tmp.pos
sort -u $PATIENT.tmp.pos | sort -n -k2,2 | sort -n -s -k1,1 > $PATIENT.tmp.pos2
mv $PATIENT.tmp.pos2 $PATIENT.tmp.pos

for i in `cat $2`
do
    BAM=$i
    TUMOR=`basename $i`
    samtools mpileup -B -R -q 40 -Q 20 -l $PATIENT.tmp.pos -f $REF $BAM > $TUMOR.pileup
    cat $SNV | awk -F"\t" '{print $1"_"$2}' > $TUMOR.tmp.pos
    cat $TUMOR.pileup | awk -F"\t" '{print $1"_"$2}' > $TUMOR.tmp.pileup.pos
    paste $TUMOR.tmp.pileup.pos $TUMOR.pileup > $TUMOR.tmp.pileup
    fgrep -f $TUMOR.tmp.pos $TUMOR.tmp.pileup > $TUMOR.tmp.pileup2
    cut -f 2- $TUMOR.tmp.pileup2 > $TUMOR.pileup
    rm $TUMOR.tmp*
done

rm $PATIENT.tmp.pos
