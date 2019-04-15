#!/bin/bash

##Allow the use of two normal bam files, one matched normal for calling of germline heterozygous SNVs, one possibly pooled normal for calculating depth ratios

if [ $# -lt 3 ]
then
    echo "missing arguments normal.bam tumor.bam ref.fa dbsnp.pos [target.bed]"
    exit -1
fi

module load samtools/1.2
module load bcftools/1.2
module load java/7u03
module load r/3.1.1
module load snpeff/3.6

SCRIPTDIR1=$mCRC_DIR/tools
SCRIPTDIR2=$mCRC_DIR/scripts/data_processing
SNPSIFT=$SnpSiftJar
##NORMALBAM=NORMALBAM1|NORMALBAM2
##NORMALBAM1 for calling SNVs
##NORMALBAM2 for calculating depth ratios
##If no "|", use the same bam file for both purpose
NORMALBAM=$1
TUMORBAM=$2
REF=$3
DBSNP=$4
BED=$5

NORMALBAM1=`echo $NORMALBAM | cut -d "|" -f 1`
NORMALBAM2=`echo $NORMALBAM | cut -d "|" -f 2`

NORMAL1=`basename $NORMALBAM1 | cut -d "." -f 1`
NORMAL2=`basename $NORMALBAM2 | cut -d "." -f 1`

TUMOR=`basename $TUMORBAM | cut -d "." -f 1`

if [ ! -s $NORMAL1.titancna.vcf ]
then
    samtools mpileup -R -q 1 -uv -I -f $REF -l $DBSNP $NORMALBAM1 | bcftools call -v -c - | java -Xmx2g -jar $SNPSIFT filter "isHet( GEN[0] ) & ( DP >= 6  ) & ( DP4[2] > 0 ) & ( DP4[3] > 0 )" > $NORMAL1.titancna.vcf
fi

if [ ! -s $NORMAL2.readcount.wig ]
then
    $SCRIPTDIR1/HMMcopy/bin/readCounter -w 1000 -q 1 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $NORMALBAM2 > $NORMAL2.readcount.wig
fi

if [ ! -s $TUMOR.readcount.wig ]
then
    $SCRIPTDIR1/HMMcopy/bin/readCounter -w 1000 -q 1 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $TUMORBAM > $TUMOR.readcount.wig
fi

$SCRIPTDIR2/TitanCNA_01.R 2.0 $NORMALBAM $TUMORBAM $BED
