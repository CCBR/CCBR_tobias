#!/bin/bash

FOLDER=$1
cd $FOLDER
for f in `ls *.bed`;do
	bedSort $f $f
	bgzip -f --threads 2 $f
	tabix -f -p bed ${f}.gz
done