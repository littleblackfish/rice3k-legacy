all : nipponbare

nipponbare : MSU7/README MSU7/all.chrs.con MSU7/all.chrs.gff3 

MSU7/README :
	mkdir -p $(@D)
	cd $(@D); wget -N ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/README
MSU7/all.chrs.con : 
	mkdir -p $(@D)
	cd $(@D); wget -N ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.chrs.con


MSU7/all.gff3 : 
	mkdir -p $(@D)
	cd $(@D); wget -N ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3


MSU7/all.chrs.gff3 : MSU7/all.gff3
	cd $(@D); grep -v 'ChrSy' all.gff3 |grep -v 'ChrUn' |gt gff3 -sort -tidy -retainids -o all.chrs.gff3

MSU7/cds_regions : MSU7/all.chrs.gff3 
	grep -P '\tCDS\t' $< | cut -f 1,4,5 |sed 's/Chr//' | sort -n -k1 | uniq | awk '{print $$1,$$2,$$3,NR}' > $(@D)/cds_regions

NB-core_v4/NB-core_v4.%.gz : 
	mkdir -p $(@D)
	cd $(@D); wget -N https://s3.amazonaws.com/3kricegenome/reduced/$(@F)

NB-core_v4/core-filtered.%.gz :  NB-core_v4/NB-core_v4.map.gz  NB-core_v4/NB-core_v4.ped.gz  MSU7/cds_regions
	cd $(@D);\
       		gunzip 	NB-core_v4.map.gz -c > tmp.map;\
		gunzip  NB-core_v4.ped.gz -c > tmp.ped;\
		plink --make-set ../MSU7/cds_regions --gene-all --file tmp --make-bed --recode --out core-filtered; \
		gzip core-filtered.ped ; \
		gzip core-filtered.map ;\
		rm tmp.???

3krg_filt_snp_v4/3krg_filt_snp_v4.%.gz : 
	cd $(@D); wget -N https://s3.amazonaws.com/3kricegenome/reduced/$(@F)
	
%/filtered.%.gz :  NB-core_v4/NB-core_v4.map.gz  NB-core_v4/NB-core_v4.ped.gz  MSU7/cds_regions
	cd $(@D);\
       		gunzip 	$(@D).map.gz -c > tmp.map;\
		gunzip  $(@D).ped.gz -c > tmp.ped;\
		plink --make-set ../MSU7/cds_regions --gene-all --file tmp --make-bed --recode --out filtered; \
		gzip filtered.ped ; \
		gzip filtered.map ;\
		rm tmp.???
