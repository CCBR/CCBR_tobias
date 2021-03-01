#!/bin/bash
module load uropa
#uropa --input uropa.config --prefix merged --threads 2 --log log
cut -f 1-4,7-13,16-19 merged_finalhits.txt > merged_finalhits_sub
head -n 1 merged_finalhits_sub > merged_header
tail -n +2 merged_finalhits_sub > merged_peaks
