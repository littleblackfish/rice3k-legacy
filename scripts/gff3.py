def gff3_iterator(fname) :
    f = open(fname, 'r')
    print 'parsing file :', fname

    for line in f :
        if line[0] == '#' :
            continue
        else : 

            tmp = line.strip().split("\t")
            try :
                attrstring = tmp[8] 

            except : 
                print line
                print 'WTF!!'

            # unpack attributes to dictionary
            attrstring = attrstring.split(";")
            attrdict = dict()

            for attr in attrstring :
                key, value = attr.split("=") 
            #    value= value.split(",")
                attrdict[key]=value

#            feature = {'seqid':tmp[0], 'source':tmp[1], 'type':tmp[2], 'range':(int(tmp[3]),int(tmp[4])), 'score':tmp[5], 'strand':tmp[6], 'phase':tmp[7], 'attr':attrdict }
            
            # don't extract source, score and phase
            feature = {'seqid':tmp[0], 'type':tmp[2], 'range':(int(tmp[3]),int(tmp[4])), 'strand':tmp[6], 'attr':attrdict }

            yield feature

def gff3_parser(fname) :
    genome=dict()
    for f in gff3_iterator(fname) :

        # we don't use source or score for this application
        #f.pop('source')
        #f.pop('score')
        
        typ = f.pop('type')
        sid = f.pop('seqid')
        fid = f['attr'].pop('ID', None)
        
        if typ == 'gene' :
            #f.pop('phase') 

            if sid not in genome :
                genome[sid] = dict()
            
            f['mRNA']=dict()

            genome[sid][fid] =  f
        else : 
            
            parent = f['attr'].pop('Parent')

            if typ == 'mRNA' :
                #f.pop('phase')
                f['CDS'] = list()
                f['exon'] = list()
                genome[sid][parent]['mRNA'][fid.split('.')[1]] = f

            else :
                if typ == 'CDS' :
                    gene,no = parent.split('.')
                    genome[sid][gene]['mRNA'][no][typ].append(f['range'])

                if typ == 'exon' : 
                    pass
                if typ == 'five_prime_UTR' :
                    pass
                if typ == 'three_prime_UTR' :
                    pass


    return genome
            

    
