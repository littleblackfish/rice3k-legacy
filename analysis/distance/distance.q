#PBS -k o
#PBS -N distmat-b60
#PBS -l nodes=1:ppn=16,walltime=4:00:00:00

projdir=/N/dc2/scratch/muroztur/rice3k
data=${projdir}/data

gff3=${data}/MSU7/all.chrs.gff3
fasta=${data}/MSU7/all.chrs.con
map=${data}/3krg_filt_snp_v4/coding.map.gz
ped=${data}/3krg_filt_snp_v4/coding.ped.gz


cd ${projdir}/analysis/distance/blosum60

python2 ${projdir}/scripts/distance_parallel.py ${gff3} ${fasta} ${map} ${ped}

