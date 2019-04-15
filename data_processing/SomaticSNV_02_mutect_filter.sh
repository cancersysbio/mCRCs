#!/bin/bash

if [ $# -lt 1 ]
then
    echo "missing arguments: <tumor.bam>"
    exit -1
fi

SCRIPTDIR=$mCRC_DIR/scripts/data_processing
REFDIR=$mCRC_DIR/ref_genomes
REF=$REFDIR/broad.Homo_sapiens_assembly19.fa
COSMIC=$REFDIR/COSMICv75_GRCh37_sorted.vcf
DBSNP=$REFDIR/dbsnp_138.b37.vcf

BAMNAME=`basename $1`

##Filters
COVERAGE=10
MINMUT=3

echo -ne "chr\tpos\tref\tvar\tnormal_reads\ttumor_reads\tnormal_varfrac\ttumor_varfrac\n" > $BAMNAME.mutect.filtered.txt
cat $BAMNAME.mutect | sed 's/chr//g' | awk -F"\t" '{ if( ( ( $1 >= 1 && $1 <= 22) || $1 == "X" || $1 == "Y" ) && $38+$39 >= '$COVERAGE' && $26+$27 >= '$COVERAGE' && $27 >= '$MINMUT' && $10 == "COVERED" && $51 == "KEEP" && $4 != $5) {print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $38+$39 "\t" $26+$27 "\t" $36 "\t" $22 } }' >> $BAMNAME.mutect.filtered.txt

cat /dev/null > $BAMNAME.mutect.blast.fna

echo -ne "Checking paralogous reads..."

for j in `seq 1 22`
do
    echo -ne " chr$j"
    for i in `cat $BAMNAME.mutect.filtered.txt | sed '1d' | awk -F"\t" '{ if($1 == "'$j'") print $1 "_" $2 }'`
    do 
	CHR=$j
	POS=`echo $i | awk -F_ '{ print $2 }'`
	echo ">${CHR}_${POS}" >> $BAMNAME.mutect.blast.fna
	PILEUP=$REFDIR/broad.Homo_sapiens_assembly19.pileup.chr${j}.pileup_filtered.pileup
	cat $PILEUP | egrep -A 20 -B 19 "^${POS}[[:space:]]" | awk '{ printf "%s", $2 } END { print "" }' >> $BAMNAME.mutect.blast.fna
    done
done

echo ""

BLASTREF=$REFDIR/blastDB.2.2.30/broad.Homo_sapiens_assembly19
$SCRIPTDIR/TOOLS.blastSEQ.sh $BAMNAME.mutect.blast.fna $BLASTREF

cat /dev/null > $BAMNAME.mutect.filtered.para.txt
head -n 1 $BAMNAME.mutect.filtered.txt > $BAMNAME.mutect.filtered.para.txt

for i in `cat $BAMNAME.mutect.blast.fna.out | sed '1d' | awk '{ print $1 }' | sort | uniq`
do
    N=`cat $BAMNAME.mutect.blast.fna.out | grep $i | wc -l`
    if [ $N -gt 1 ]
    then
	CHR=`echo $i | awk -F_ '{ print $1 }'`
	POS=`echo $i | awk -F_ '{ print $2 }'`

	cat $BAMNAME.mutect.filtered.txt | awk '{ if($1 == "'$CHR'" && $2 == "'$POS'") print $0 }' >> $BAMNAME.mutect.filtered.para.txt
    fi
done

cat $BAMNAME.mutect.filtered.txt > $BAMNAME.mutect.tmp
for i in `cat $BAMNAME.mutect.filtered.para.txt | sed '1d' | awk '{ print $1 "_" $2 }'`
do
    cat $BAMNAME.mutect.tmp | awk -F"\t" '{ if($1 "_" $2 != "'$i'") print $0 }' > $BAMNAME.mutect.tmp2
    mv $BAMNAME.mutect.tmp2 $BAMNAME.mutect.tmp
done

mv $BAMNAME.mutect.tmp $BAMNAME.mutect.filtered.nopara.txt
