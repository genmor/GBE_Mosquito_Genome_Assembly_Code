#!/bin/bash
mkdir 4_simple_faa_headers
cd 4_simple_faa_headers
cp ../4_no_stop/*.faa .
sed -i '/^>/s/ .*//' *.faa

###append filename to each sequence header
###do this if you're running orthofinder
#faa=*.faa
# for i in ${faa}
	# do
		# sed -i "s/^>/>${i%.*}_/" "$i"
# done

