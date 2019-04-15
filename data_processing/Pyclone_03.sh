#!/bin/bash

if [ $# -lt 1 ]
then
    echo "missing arguments: <patient name>"
    exit -1
fi

source $mCRC_DIR/tools/venv/bin/activate

PyClone=$mCRC_DIR/tools/PyClone-0.12.7/PyClone

PATIENT=$1

YAML=${PATIENT}_config.yaml
echo $PATIENT $YAML

$PyClone analyse $YAML
$PyClone build_table $YAML $PATIENT.pyclone.out
