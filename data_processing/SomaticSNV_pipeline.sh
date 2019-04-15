#!/bin/bash

module load samtools/1.2
module load java/7u03
module load mutect/1.1.7
module load blast/2.2.30+
module load bam-readcount/0.7.4
module load varscan2/2.3.9
module load annovar/20150617

if [ $# -lt 2 ]
then
    echo "missing arguments: <normal.bam> <tumor.bam>"
    exit -1
fi

SCRIPTDIR=$mCRC_DIR/scripts/data_processing
REFDIR=$mCRC_DIR/ref_genomes
REF=$REFDIR/broad.Homo_sapiens_assembly19.fa
COSMIC=$REFDIR/COSMICv75_GRCh37_sorted.vcf
DBSNP=$REFDIR/dbsnp_138.b37.vcf

$SCRIPTDIR/SomaticSNV_01_mutect.sh $1 $2
$SCRIPTDIR/SomaticSNV_02_mutect_filter.sh $2
$SCRIPTDIR/SomaticSNV_03_varscan_filter.sh $2
$SCRIPTDIR/SomaticSNV_04_annovar.sh $2
