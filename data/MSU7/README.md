MSU rice genome assembly and annotation (MSU7/IRGSP1.0)
release 7


main page 

http://rice.plantbiology.msu.edu/annotation_pseudo_current.shtml

getting the genome

Chromosome-ony version omitting ChrSy and ChrUn segments
ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.chrs.con

full versions
ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.con

proteins
ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.cdna


getting the annotation

ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3

This will remove ChrSy and ChrUn which are unplaced pseudomolecules and sanitize the annotation 

grep -v 'ChrSy' all.gff3 |grep -v 'ChrUn' |gt gff3 -sort -tidy -retainids -o all.chrs.gff3

