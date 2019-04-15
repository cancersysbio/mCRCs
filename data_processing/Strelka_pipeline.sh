#!/bin/bash

module load samtools/1.2
module load annovar/20150617
module load vcftools/0.1.13

SCRIPTDIR=$mCRC_DIR/scripts/data_processing
REFDIR=$mCRC_DIR/ref_genomes
REF='broad.Homo_sapiens_assembly19.fa'
STRELKA_INSTALL_DIR='scripts/strelka'

if [ $# -lt 2 ]
then
    echo "missing arguments: <normal.bam> <tumor.bam>"
    exit -1
fi

NORMALBAM=$1
TUMORBAM=$2

echo $NORMALBAM $TUMORBAM

$SCRIPTDIR/Strelka_01.sh $NORMALBAM $TUMORBAM
$SCRIPTDIR/Strelka_02_annovar.sh $NORMALBAM $TUMORBAM
