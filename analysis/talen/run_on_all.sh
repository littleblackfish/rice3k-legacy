#!/bin/bash

#run for reference samples 

for file in /DATA/GROUP/huangcy/Rice_3000_Assembly/VCF/*snp.vcf ;
	do python2 identify.py $file >> reference.results;
done
   


#run for unknown samples

for file in /DATA/GROUP/huangcy/rice_genomes_BY/Illumina8/VCF/*snp.vcf ;
	do python2 identify.py $file >> unknown.results;
done
