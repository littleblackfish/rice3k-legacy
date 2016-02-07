#! /bin/bash

#generate regionlist from annotation
grep $'\tCDS\t' ../MSU7/all.chrs.gff3 | cut -f 1,4,5 |sed 's/Chr//' | sort -n -k1 | uniq | awk '{print $1,$2,$3,NR}'> cds_regions


#filter snps by cds
plink --make-set cds_regions --gene-all --file $1 --make-bed --recode --out filtered

gzip filtered.ped
gzip filtered.map

rm filtered.???
