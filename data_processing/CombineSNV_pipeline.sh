#!/bin/bash

##Combine SNV called from multiple samples of the same patient

module load samtools/1.2
module load java/7u03
module load varscan2/2.3.9
module load bam-readcount/0.7.4

if [ $# -lt 2 ]
then
    echo "missing arguments: <patient name> <tumorbam.list>"
    exit -1
fi

SCRIPTDIR=$mCRC_DIR/scripts/data_processing
REFDIR=$mCRC_DIR/ref_genomes
REF=$REFDIR/broad.Homo_sapiens_assembly19.fa

$SCRIPTDIR/CombineSNV_01_varscan_filter_multi_sample.sh $1 $2
$SCRIPTDIR/CombineSNV_02_pileup.sh $1 $2
$SCRIPTDIR/CombineSNV_03_generate_table.R $1 $2
