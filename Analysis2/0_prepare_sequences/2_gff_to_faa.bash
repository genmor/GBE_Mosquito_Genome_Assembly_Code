#!/bin/bash

mkdir 2_gff_to_faa

cd 1_longest_iso_gff
gff=*.gff

for i in ${gff}
	do
		agat_sp_extract_sequences.pl -g $i -f ../fna/${i%.gff}.fna -p -o ../2_gff_to_faa/${i%.gff}.faa
done

##note that for some genomes (i.e.., Bf05), the fna will need to be folded prior
