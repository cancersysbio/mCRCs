#!/bin/bash

source ~/venv/bin/activate

PyClone='scripts/PyClone-0.12.7/PyClone'

rm -f *_config.yaml

for i in *mut
do
    $PyClone build_mutations_file --ref_prior normal_variant --var_prior parental_copy_number $i $i.yaml
    PATIENT=`echo $i | awk -F"_" '{print $1}'`
    SAMPLE=`basename $i _SNV.mut`
    PURITY=`grep $SAMPLE mCRC_All_Purity.txt | cut -f 2`
    if [ ! -e ${PATIENT}_config.yaml ] 
    then
	cp config_base.yaml ${PATIENT}_config.yaml
	echo -e "trace_dir: $PATIENT\n" >> ${PATIENT}_config.yaml
	echo "samples:" >> ${PATIENT}_config.yaml
    fi
    echo "  $SAMPLE:" >> ${PATIENT}_config.yaml
    echo "    mutations_file: $i.yaml" >> ${PATIENT}_config.yaml
    echo "    tumour_content:" >> ${PATIENT}_config.yaml
    echo "      value: $PURITY" >> ${PATIENT}_config.yaml
    echo "    error_rate: 0.001" >> ${PATIENT}_config.yaml
done
