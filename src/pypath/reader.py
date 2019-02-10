import os

import pypath.common as common
import pypath.dataio as dataio
import curl


class ReaderBase(object):


    def __init__(self, settings):
        self.settings = settings


    def __iter__(self):
        for i, row in enumerate(self.resource):
            if i < self.settings.header or not row:
                continue

            if isinstance(row, common.basestring):
                if row[0] == "#":
                    continue

                row = row.strip('\n\r ').split(self.settings.separator)

            yield row


    def setup_resource(self):
        self.input = self.settings.inFile

        if callable(self.input):
            self.resource = self.input(**self.settings.inputArgs)

        elif isinstance(self.input,common.basestring):
            if hasattr(dataio, self.input):
                self.resource = getattr(dataio, self.input)(
                    **self.settings.inputArgs
                )
            elif (
                os.path.exists(self.input) or
                curl.is_url(self.input)
            ):
                c = curl.Curl(self.input, **self.settings.curlArgs)
                self.resource = c.result

        elif hasattr(self.input, '__iter__'):
            self.resource = self.input

        else:
            self.resource = []




class Reader(ReaderBase):
    """
    Reads network data from an iterable according to
    an :py:class:``pypath.input_formats.ReadSettings`` object.
    """
    def __init__(self, settings):
        """
        :arg input_formats.ReadSettings settings:
            ReadSettings instance
        """
        ReaderBase.__init__(self, settings)


    def __iter__(self):
        for row in ReaderBase.__iter__(self):
            yield Row(row)


class FieldProcessor(object):
    def __init__(self, field, method = None):
        """
        :arg field:
            Field processing definition
        :arg callable method:
            Custome method which accepts a row as an argument and returns
            the processed field value
        """
        self.field = field
        self.method = method
        self.failed = False
        self.setup_method()


    def __iter__(self):

        pass


    def setup_method(self):
        if callable(self.method):
            self._method = self.method

        else:


    def process_commom(self):
        pass




class Row(object):
    __slots__ = ["row", "failed"]

    def __init__(self, row):
        self.row = row
        self.failed = False


    def __getitem__(self, i):
        """
        :arg int i:
            Index of the field
        """
        try:
            return self.row[i]

        except IndexError:
            self.failed = True
            return ""





class EdgeAttribute(FieldProcessor):
    def __init__(self,
                 separator=None,
                 isDirected=False,
                 sign=False,
                 references=False,
                 must_have_references=True,
                 extraEdgeAttrs={},
                 interactionType='PPI',
                 positiveFilters=[],
                 negativeFilters=[],
                 inFile=None,
                 header=False
                 ):

        Reader.__init__(self,
                 separator=None,
                 isDirected=False,
                 sign=False,
                 references=False,
                 must_have_references=True,
                 extraEdgeAttrs={},
                 interactionType='PPI',
                 positiveFilters=[],
                 negativeFilters=[],
                 inFile = None,
                 header = False
                        )


        self.isDirected = isDirected
        self.extraEdgeAttrs = extraEdgeAttrs
        self.separator = separator
        self.refs = references
        self.must_have_references = must_have_references and bool(references)
        self.sign = sign
        self.intType = interactionType
        self.positiveFilters = positiveFilters
        self.negativeFilters = negativeFilters
        self.inFile = inFile
        self.header = header



class NodeAttribute(FieldProcessor):
    def __init__(self,
                 name="unknown",
                 separator=None,
                 nameColA=0,
                 nameColB=1,
                 nameTypeA="uniprot",
                 nameTypeB="uniprot",
                 typeA="protein",
                 typeB="protein",
                 inFile=None,
                 extraNodeAttrsA={},
                 extraNodeAttrsB={},
                 header=False,
                 taxonA=9606,
                 taxonB=9606,
                 ncbiTaxId=False,
                 mark_source=None,
                 mark_target=None,
                 inputArgs={},
                 huge=False,
                 resource=None,
                 positiveFilters=[],
                 negativeFilters=[]
                 ):

        self.typeA = typeA
        self.typeB = typeB
        self.nameColA = nameColA
        self.nameColB = nameColB
        self.nameTypeA = nameTypeA
        self.nameTypeB = nameTypeB
        self.extraNodeAttrsA = extraNodeAttrsA
        self.extraNodeAttrsB = extraNodeAttrsB
        self.name = name
        self.separator = separator
        self.taxonA = taxonA
        self.taxonB = taxonB
        self.ncbiTaxId = ncbiTaxId
        self.positiveFilters = positiveFilters
        self.negativeFilters = negativeFilters
        self.inputArgs = inputArgs
        self.huge = huge
        self.resource = self.name if resource is None else resource
        self.mark_source = mark_source
        self.mark_target = mark_target
        self.inFile = inFile
        self.header = header


class MergeTable():
    def __init__(self):



class CreateIgraphObject():
    def __init__(self):











