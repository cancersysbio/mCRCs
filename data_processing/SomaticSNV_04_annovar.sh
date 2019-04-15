#!/bin/bash

if [ $# -lt 1 ]
then
    echo "missing arguments: <tumor.bam>"
    exit -1
fi

REFDIR=$mCRC_DIR/ref_genomes
REF=$REFDIR/broad.Homo_sapiens_assembly19.fa
COSMIC=$REFDIR/COSMICv75_GRCh37_sorted.vcf
DBSNP=$REFDIR/dbsnp_138.b37.vcf

BAMNAME=`basename $1`

##Annovar
MUTECTNAME=$BAMNAME.mutect.filtered.nopara.txt
INPUT=$MUTECTNAME.varscan.filtered

cat $INPUT | sed 1d | awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4}' > $INPUT.tmp.input

table_annovar.pl $INPUT.tmp.input $ANNOVAR/humandb -buildver hg19 -out $INPUT.tmp.input.out -remove -protocol refGene,snp138,1000g2015aug_all,cosmic70,ljb26_all -operation g,f,f,f,f -nastring NA

mv $INPUT.tmp.input.out.hg19_multianno.txt $INPUT.annovar
rm $INPUT.tmp.input
