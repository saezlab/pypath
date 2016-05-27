class ReferenceList(object):
    
    def __init__(self,nameType,typ,tax,inFile,**kwargs):
        self.infile = inFile
        self.nameType = nameType
        self.typ = typ
        self.tax = tax
        self.kwargs = kwargs
    
    def load(self):
        if hasattr(dataio, self.infile):
            toCall = getattr(dataio, self.infile)
            lst = toCall(**self.kwargs)
        else:
            f = codecs.open(self.infile,encoding='utf-8',mode='r')
            lst = []
            for l in f:
                lst.append(l.strip())
            f.close()
        self.lst = set(lst)

def get_reflists():
    return [
        ReferenceList('uniprot', 'protein', 9606, 'all_uniprots')
    ]
