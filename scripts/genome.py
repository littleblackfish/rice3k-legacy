from parsers import *


class genome :
    def __init__ (self, gff3fname, fastafname=None, mapfname=None) :
        # a genome definitely has an annotation
        self.feat = gff3_parser(gff3fname)

        # populate a gene dict for fast lookup of genes
        self.geneindex = dict()
        for c in self.feat : 
            chromosome = self.feat[c]
            for gene in chromosome.iterkeys() :
                self.geneindex[gene]=self.feat[c][gene]

        self.genelist = self.geneindex.keys() 

        # it can also have an actual sequence
        if fastafname :
            self.sequence = fasta_parser(fastafname)
            for c in self.feat :
                self.feat[c]['sequence']=self.sequence[c]
            
            print '# Extracting sequences for genes...'
            for gene in self.geneindex.itervalues() : 
                for mrna in gene['mRNA'].itervalues() : 
                    mrna['seq']=[]
                    for interval in mrna['CDS'] : 
                        mrna['seq'].append(self.sequence[mrna['seqid']][interval[0]-1:interval[1]])

            self.hasSequence=True
        else :
            self.hasSequence = False

        # and a MAP file listing SNPs
        if mapfname :
            self.map_snps(mapfname)
            self.hasSNPs = True
        else : 
            self.hasSNPs = False

    def __getitem__(self,key) : 
        if isinstance(key,str) :
            #assert key in geneindex.iterkeys()
            return self.geneindex[key]
        elif isinstance(key,int) :
            return self.feat[key]
        else : raise KeyError('Use int as chromosome no or string as gene name')

    def __iter__(self) : 
        for gene in self.geneindex.itervalues() :
            yield gene 
    
    def __len__ (self) : return len(self.genelist)

    def __repr__ (self) : return 'Genome with {} genes in {} chromosomes'.format(len(self.geneindex), len(self.feat))

    # assembles a product from a given list of cds
    def cds_assemble(self,cdslist) :
        cds = Seq('', generic_dna)
        for s in cdslist :
            cds += s
        return cds

            
    def protein_iter(self) :
        for gene in self.geneindex.itervalues() :
            for protein in gene['mRNA'].itervalues() :
                tmp = self.cds_assemble (protein['seq']) 
                if gene['strand'] == '-' :
                    tmp=tmp.reverse_complement()
                yield tmp

    # maps SNPs from a MAP file to CDS segments
    # this way we know exactly where to look in the PED file for each CDS

    def map_snps(self, mapfname):
        mapdict = map_parser(mapfname)
        print '# Mapping snps to genes...' 
        for gene in self.geneindex.itervalues() :
            for mrna in gene['mRNA'].itervalues() :
                mrna['SNPind']=[]
                mrna['SNPpos']=[]
                mrna['SNPpos-old']=[]
        
                nCDS = len(mrna['CDS']) # number of ranges to append
                
                seqCDS = Seq('',generic_dna) # empty sequence
                

                for i in range(nCDS) : 
                    # interval for this cds segment
                    interval = mrna['CDS'][i]

                    # find SNPs in this interval
                    indices, positions = map_find_loci(mapdict, mrna['seqid'], interval)

                    # PED indices remain unchanged 
                    mrna['SNPind'] += indices
                    
                    # offset and append positions (on the sequence)
                    mrna['SNPpos'] += [ p - interval[0] + len(seqCDS) for p in positions]
                    mrna['SNPpos-old'].append([p-interval[0] for p in positions])

                    # append cds sequence
                    seqCDS += mrna['seq'][i]

                mrna['seq-full'] = seqCDS 
