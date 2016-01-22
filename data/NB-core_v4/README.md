data from 

http://oryzasnp-atcg-irri-org.s3-website-ap-southeast-1.amazonaws.com/

990K CoreSNP dataset, called vs Nipponbare MSU7/IRGSP1.0 genome, PLINK ped file
https://s3.amazonaws.com/3kricegenome/reduced/NB-core_v4.ped.gz

990K CoreSNP dataset, called vs Nipponbare MSU7/IRGSP1.0 genome, PLINK map file
https://s3.amazonaws.com/3kricegenome/reduced/NB-core_v4.map.gz

https://s3.amazonaws.com/3kricegenome/reduced/990k_3krg-snp-README.txt

Core SNP set (v.4) README
======================
3kRG core SNP set v.4 is a subset of the 3kRG filtered SNP set v4 obtained using 2-step LD pruning procedure:

1) LD pruning with window size 10kb, window step: 1 snp, R2 threshold: 0.8, followed by

2) LD pruning with window size 500 snps, window step 1 snp, R2 threshold 0.8

Total SNPs: 996,009
samples : 3023

39d8c264feaa75c9ee4c9b59e7d38c074c8d08f1  NB-core_v4.map.gz
17238669e38bc13831d8d78f3c1f7deca48c229d  NB-core_v4.ped.gz
