#!/bin/bash

if [ $# -lt 1 ]
then
    echo "missing arguments: <patient name>"
    exit -1
fi

module load java/7u03

LICHEE=$mCRC_DIR/scripts/LICHeE/release/lichee.jar

PATIENT=$1
INPUT=${PATIENT}_SNV_Lichee_Input.txt

LOGFILE=${PATIENT}_Lichee.log

echo $PATIENT > $LOGFILE
echo "java -Xmx2g -jar $LICHEE -build -i $INPUT -o ${PATIENT}_Lichee.trees.txt -maxVAFAbsent 0.05 -minVAFPresent 0.3 -n 0 -s 1 -v -e 0.1" >> $LOGFILE
java -Xmx2g -jar $LICHEE -build -i $INPUT -o ${PATIENT}_Lichee -maxVAFAbsent 0.05 -minVAFPresent 0.3 -n 0 -s 1 -v -e 0.1 -showTree 1 &>> $LOGFILE
