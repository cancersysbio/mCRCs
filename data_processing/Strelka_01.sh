#!/bin/bash

if [ $# -lt 2 ]
then
    echo "missing arguments: <normal.bam> <tumor.bam>"
    exit -1
fi

REFDIR=$mCRC_DIR/ref_genomes
REF=broad.Homo_sapiens_assembly19.fa

NORMALBAM=$1
TUMORBAM=$2

TUMOR=`basename $TUMORBAM .bam`
OUTPUTDIR=${TUMOR}_strelka 

echo Call $TUMOR

rm -rf $OUTPUTDIR

$STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl --normal=$NORMALBAM --tumor=$TUMORBAM --ref=$REFDIR/$REF --config=$STRELKA_INSTALL_DIR/config.ini --output-dir=$OUTPUTDIR

cd $OUTPUTDIR
make -j 1

cp $OUTPUTDIR/results/passed.somatic.indels.vcf $TUMOR.strelka.indels.vcf
cp $OUTPUTDIR/results/passed.somatic.snvs.vcf $TUMOR.strelka.snvs.vcf
