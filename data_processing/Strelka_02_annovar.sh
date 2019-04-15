#!/bin/bash

if [ $# -lt 2 ]
then
    echo "missing arguments: <normal.bam> <tumor.bam>"
    exit -1
fi

NORMALBAM=$1
TUMORBAM=$2

TUMOR=`basename $TUMORBAM .bam`
VCF=$TUMOR.strelka.indels.vcf
VCF2=$TUMOR.strelka.indels

echo Annotate $TUMOR

grep ^#CHROM $VCF > $VCF.avinput
grep -v ^# $VCF | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9":GT\t"$10":0/1\t"$11":0/1\t"}' >> $VCF.avinput
convert2annovar.pl -format vcf4 $VCF.avinput -allsample -includeinfo -outfile $VCF
mv $VCF.TUMOR.avinput $VCF.avinput
rm $VCF.NORMAL.avinput

table_annovar.pl $VCF.avinput $ANNOVAR/humandb/ -buildver hg19 -out $VCF.annovar -remove -protocol refGene,snp138,1000g2015aug_all,cosmic70 -operation g,f,f,f -nastring NA
mv $VCF.annovar.hg19_multianno.txt $VCF2.annovar
rm $VCF.avinput

echo -e "Tumor_cover\tTumor_ref_count\tTumor_alt_count\tTumor_alt_freq" > $VCF2.tmp2
vcf-query -c TUMOR -f '[%DP\t%TAR\t%TIR\t]\n' $VCF | sed s/,[0-9]*//g | awk '{print $0$3/($1)}' >> $VCF2.tmp2 

echo -e "Normal_cover\tNormal_alt_freq" > $VCF2.tmp22
vcf-query -c NORMAL -f '[%DP\t%TAR\t%TIR\t]\n' $VCF | sed s/,[0-9]*//g | awk '{print $1"\t"$3/($1)}' >> $VCF2.tmp22

cut -f 1-5 $VCF2.annovar > $VCF2.tmp1
cut -f 6- $VCF2.annovar > $VCF2.tmp3
paste $VCF2.tmp1 $VCF2.tmp22 $VCF2.tmp2 $VCF2.tmp3 > $VCF2.annovar
rm $VCF2.tmp*
