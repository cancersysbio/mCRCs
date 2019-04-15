#!/bin/bash

if [ $# -lt 2 ]
then
    echo "missing arguments: <patient name> <tumorbam.list>"
    exit -1
fi

REFDIR=$mCRC_DIR/ref_genomes
REF=$REFDIR/broad.Homo_sapiens_assembly19.fa
VARSCAN=$VARSCAN2ROOT/VarScan.v2.3.9.jar

PATIENT=$1

cat /dev/null > $PATIENT.tmp.pos
for i in `cat $2`
do
    TUMORBAM=$i
    TUMOR=`basename $TUMORBAM`

    MUTECT=$TUMOR.mutect
    MUTECT_F=$TUMOR.mutect.filtered.nopara.txt.varscan.filtered

    cat $MUTECT_F | sed 1d |  awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4}' >> $PATIENT.tmp.pos
done

cut -f 1-4 $PATIENT.tmp.pos | sort -u | sort -n -k2,2 | sort -n -s -k1,1 > $PATIENT.tmp.pos2
mv $PATIENT.tmp.pos2 $PATIENT.tmp.pos



cat $PATIENT.tmp.pos | awk '{print $1"\t"$2"\t"$2}' > $PATIENT.tmp.pos2
cat $PATIENT.tmp.pos | awk '{print $1"\t"$2-1"\t"$2}' > $PATIENT.tmp.bed

cat /dev/null > $PATIENT.tmp.bam.list

for i in `cat $2`
do
    TUMORBAM=$i
    TUMOR=`basename $TUMORBAM`

    MUTECT=$TUMOR.mutect
    MUTECT_F=$TUMOR.mutect.filtered.nopara.txt.varscan.filtered

    echo $TUMOR $TUMORBAM

    samtools view -b -q 0 -F 1536 -x BQ -o $PATIENT.tmp.$TUMOR.bam -L $PATIENT.tmp.bed $TUMORBAM
    samtools index $PATIENT.tmp.$TUMOR.bam

    echo "$PATIENT.tmp.$TUMOR.bam" >> $PATIENT.tmp.bam.list
done
    
if [ `wc -l $PATIENT.tmp.bam.list | cut -d " " -f 1` -gt 1 ]
then
    samtools merge -b $PATIENT.tmp.bam.list $PATIENT.tmp.bam
else
    cp $PATIENT.tmp.$TUMOR.bam $PATIENT.tmp.bam
fi
samtools index $PATIENT.tmp.bam
    
bam-readcount -q 0 -w 3 -l $PATIENT.tmp.pos2 -f $REF $PATIENT.tmp.bam > $PATIENT.tmp.readcount
    
java -Xmx2g -jar $VARSCAN fpfilter $PATIENT.tmp.pos $PATIENT.tmp.readcount --keep-failures 1 --output-file $PATIENT.combined.varscanfilter --min-var-count 3 --min-var-freq 0.0 --min-ref-basequal 20 --min-var-basequal 20 --max-mmqs-diff 100
    
java -Xmx2g -jar $VARSCAN fpfilter $PATIENT.tmp.pos $PATIENT.tmp.readcount --keep-failures 1 --output-file $PATIENT.tmp.varscanfilter --min-var-count 3 --min-var-freq 0.0 --min-ref-basequal 20 --min-var-basequal 20 --max-var-mmqs 100 --max-mmqs-diff 70 --min-ref-mapqual 40 --min-var-mapqual 40 --max-mapqual-diff 20
    
cut -f 21 $PATIENT.tmp.varscanfilter > $PATIENT.tmp2.varscanfilter
paste $PATIENT.combined.varscanfilter $PATIENT.tmp2.varscanfilter > $PATIENT.tmp3.varscanfilter
mv $PATIENT.tmp3.varscanfilter $PATIENT.combined.varscanfilter

rm $PATIENT.tmp*
