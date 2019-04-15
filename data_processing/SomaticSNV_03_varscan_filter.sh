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

##Varscan2 filter

BAM=$1
MUTECTNAME=$BAMNAME.mutect.filtered.nopara.txt
VARSCAN=$VARSCAN2ROOT/VarScan.v2.3.9.jar

cat $MUTECTNAME | sed 1d | awk '{print $1"\t"$2"\t"$2}' > $MUTECTNAME.tmp.pos
cat $MUTECTNAME | sed 1d | awk '{print $1"\t"$2-1"\t"$2}' > $MUTECTNAME.tmp.bed
cat $MUTECTNAME | sed 1d | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $MUTECTNAME.tmp.varscan

samtools view -b -o $MUTECTNAME.tmp.bam -L $MUTECTNAME.tmp.bed -q 0 -F 1536 -x BQ $BAM
samtools index $MUTECTNAME.tmp.bam

bam-readcount -q 0 -w 3 -l $MUTECTNAME.tmp.pos -f $REF $MUTECTNAME.tmp.bam > $MUTECTNAME.tmp.readcount

java -Xmx2g -jar $VARSCAN fpfilter $MUTECTNAME.tmp.varscan $MUTECTNAME.tmp.readcount --keep-failures 1 --output-file $MUTECTNAME.varscanfilter --min-var-count 3 --min-var-freq 0.01 --min-ref-basequal 20 --min-var-basequal 20 --max-mmqs-diff 100

java -Xmx2g -jar $VARSCAN fpfilter $MUTECTNAME.tmp.varscan $MUTECTNAME.tmp.readcount --keep-failures 1 --output-file $MUTECTNAME.tmp.varscanfilter --min-var-count 3 --min-var-freq 0.01 --min-ref-basequal 20 --min-var-basequal 20 --max-var-mmqs 100 --max-mmqs-diff 70 --min-ref-mapqual 40 --min-var-mapqual 40 --max-mapqual-diff 20

cut -f 21 $MUTECTNAME.tmp.varscanfilter > $MUTECTNAME.tmp2.varscanfilter
paste $MUTECTNAME.varscanfilter $MUTECTNAME.tmp2.varscanfilter > $MUTECTNAME.tmp3.varscanfilter
mv $MUTECTNAME.tmp3.varscanfilter $MUTECTNAME.varscanfilter

cut -f 22 $MUTECTNAME.varscanfilter > $MUTECTNAME.tmp.vf
head -n 1 $MUTECTNAME > $MUTECTNAME.tmp1
cat $MUTECTNAME | sed 1d > $MUTECTNAME.tmp2
paste $MUTECTNAME.tmp.vf $MUTECTNAME.tmp2 | grep ^PASS | cut -f 2- >> $MUTECTNAME.tmp1
mv $MUTECTNAME.tmp1 $MUTECTNAME.varscan.filtered

rm $MUTECTNAME.tmp*
