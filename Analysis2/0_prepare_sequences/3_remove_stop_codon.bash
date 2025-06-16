#!/bin/bash
mkdir 3_no_stop
cd 3_no_stop
cp ../2_gff_to_faa/*.faa .
sed -i "s/*$//g" *.faa
sed -E -i '/>/!s/\*/X/g' *.faa