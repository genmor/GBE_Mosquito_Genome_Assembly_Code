#!/bin/bash
mkdir 1_longest_iso_gff

cd 0_original

gff=*.gff

for i in $gff
do 
agat_sp_keep_longest_isoform.pl --gff ${i} -o ../1_longest_iso_gff/${i%.gff}.gff
done