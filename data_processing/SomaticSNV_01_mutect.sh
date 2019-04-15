#!/bin/bash

if [ $# -lt 2 ]
then
    echo "missing arguments: <normal.bam> <tumor.bam>"
    exit -1
fi

REFDIR=$mCRC_DIR/ref_genomes
REF=$REFDIR/broad.Homo_sapiens_assembly19.fa
COSMIC=$REFDIR/COSMICv75_GRCh37_sorted.vcf
DBSNP=$REFDIR/dbsnp_138.b37.vcf

BAMNAME=`basename $2`

##MuTect
java -Xmx2g -jar $MUTECT -rf BadCigar --analysis_type MuTect --reference_sequence $REF --dbsnp $DBSNP --cosmic $COSMIC --input_file:normal $1 --input_file:tumor $2 --out $BAMNAME.mutect
