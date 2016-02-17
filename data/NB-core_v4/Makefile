all : nipponbare

base = NB-core_v4
url = https://s3.amazonaws.com/3kricegenome/reduced/

raw : ${base}.map.gz  ${base}.ped.gz  shasum
	sha1sum -c shasum

${base}.%.gz : 
	wget -N ${url}/$(@F)

coding : raw ../MSU7/cds_regions
	which plink
	gunzip 	NB-core_v4.map.gz -c > tmp.map
	gunzip  NB-core_v4.ped.gz -c > tmp.ped
	plink --make-set ../MSU7/cds_regions --gene-all --file tmp --make-bed --recode --out coding
	gzip coding.ped 
	gzip coding.map
	rm tmp.???
