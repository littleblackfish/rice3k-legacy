all : coding

base = 3krg_filt_snp_v4
url = https://s3.amazonaws.com/3kricegenome/reduced

raw : ${base}.map.gz  ${base}.ped.gz  shasum
	sha1sum -c shasum

${base}.%.gz : 
	wget -nv -N ${url}/$(@F)

coding : raw ../MSU7/cds_regions
	which plink
	gunzip 	${base}.map.gz -c > tmp.map
	gunzip  ${base}.ped.gz -c > tmp.ped
	plink --allow-no-sex --make-set ../MSU7/cds_regions --gene-all --file tmp --make-bed --recode --out coding
	gzip coding.ped 
	gzip coding.map
	rm tmp.???
