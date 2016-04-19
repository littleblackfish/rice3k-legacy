from parsers import fasta_reference,ped_stats
from sys import argv
import gzip


print '# fasta file :', argv[1]
print '# plink file : ', argv[2]

reference = fasta_reference(argv[1], argv[2]+'.map.gz')

stats = ped_stats(gzip.open(argv[2]+'.ped.gz'), reference, nrows = 10)

print 'name\tchanging\thomo\thetero\tmissing'
for name in stats :
    this= stats[name]
    print '{}\t{:d}\t{:d}\t{:d}\t{:d}'.format(name,this['changing'], this['homo'], this['hetero'], this['missing'])




