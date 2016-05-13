#!/bin/bash

# take a vcf file
input=$1

# strip it to a 
# generate a 4 column file
output=$2


awk 'length($4) == 1 && length($5) == 1 && $5 != "\." { print $1,$2,$4,$5,$10 }' ${input} |\
	# only take lines with single reference and single alt (other than .)
       	sed 's/\:.*//' |\
	# strip coverage etc. information from genotype (last) column
	sed 's/\(Osj\)\?[Cc]hr0\?//' |\
	# strip chromosome keyword from the first column
	awk '{if ($5 == "1/1") print $1,$2,$3,$4}' > ${output}.dat
	# only print homozygous SNPs
