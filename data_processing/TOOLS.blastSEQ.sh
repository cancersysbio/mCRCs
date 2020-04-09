#!/bin/bash

#output less aligned sequences (-max_target_seqs 100)

if [ $# -lt 2 ]
then
    echo "need arg <query.fna> <ref>"
    exit -1
fi

module load blast

OPTIONS="-perc_identity 94.5 -max_target_seqs 100 -outfmt 6 -strand plus -num_threads 1"
WSIZEMID=7
WSIZELOC=6
HEADER="query\tsubject\tpercidentity\talignlen\tmismatches\tgaps\tq.start\tq.end\ts.start\ts.end\tevalue\tscore"
MINLEN=40
MINMATCH=95

echo -e $HEADER > $1.out.tmp

split -a 3 -l 100 -d $1 $1.split.
for i in $1.split.*
do
    blastn -query $i -db $2 -out $i.out $OPTIONS -word_size $WSIZEMID
    cat $i.out | awk '{ if(int($3) >= '$MINMATCH' && $4 >= '$MINLEN') print $0 }' >> $1.out.tmp
    rm $i
    rm $i.out
done

mv $1.out.tmp $1.out
