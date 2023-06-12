#!/bin/bash

FOLDER=$1
sleep 10
cd $FOLDER
nbeds=$(ls *.bed 2>/dev/null | wc -l)
if [ "$nbeds" != "0" ];then
for f in `ls *.bed`;do
        bedSort $f $f
        bgzip -f --threads 2 $f
        tabix -f -p bed ${f}.gz
done
fi