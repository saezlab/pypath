"""
Parse BRENDA data and emit Entity records.

This module converts enzyme and ligand information from BRENDA into Entity
records using the declarative schema pattern.
"""

import re
from functools import partial
from collections import defaultdict

from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.macdb import iter_associations
from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    OntologyDataset
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.internals.ontology_builder import (
    OntologyBuilder,
    RelationshipBuilder
)
from pypath.internals.ontology_schema import OntologyTypedef
from pypath.internals.cv_terms import (
    EntityTypeCv,
    ProteinFunctionalClassCv,
    InteractionTypeCv,
    OntologyCv,
    DiseaseAnnotationCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeSubtypeCv,
    MoleculeAnnotationsCv,
    AssayAnnotationsCv,
)

# =================================== SET-UP ===================================

URL = 'https://www.brenda-enzymes.org/download.php'
# XXX: As specified in inputs v1
RECORDS_ENABLED = {
    'ID', # EC-class
    'PR', # protein 	sequence identifier and name of source in {...}
    'AC', # activating compound
    'IN', # inhibitors
    'CF', # cofactor
    'KI', # Ki-value	inhibitor in {...}
    'KM', # KM-value	substrate in {...}
    'RF', # references
}

ID = re.compile(r"^#([\d,]+)#")
PR_ORGANISM_NAME = re.compile(r"#\d+# '?([-\w\s\.\[\]]+[^\s#\{<\('])")
PR_IDENTIFIER = re.compile(r".*\{(.*?)\;.*")
REFERENCE = re.compile(r".*\<([\d,]+)\>$")
AC_COMPOUND = re.compile(r"^#[\d,]+# ([^\<#]+) [\(\<].*")
AC_DETAILS = re.compile(r".*\((#.*\>)\).*")

config = ResourceConfig(
    id=ResourceCv.BRENDA,
    name='BRENDA',
    url='https://www.brenda-enzymes.org/index.php',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='41206471',
    primary_category='enzsub',
    description=(
        'BRENDA is the main collection of enzyme functional data available to '
        'the scientific community.'
    ),
)

# ================================== DOWNLOAD ==================================

def process_entry(entry):

    record = defaultdict(list)

    for line in entry:

        # Controlling for empty records
        if not line:

            break

        else:

            key, value = line

            if key == 'ID':

                # Controlling for deleted, transferred or dummy records
                if '(' in value:

                    break

                else:

                    record[key] = value

            else:

                record[key].append(value)

    return record if len(record) > 1 else {}

def process_record(record):

    eid = record['ID']

    # Processing organisms and identifiers

    proc = defaultdict(lambda: defaultdict(set))

    for pr in record.get('PR', []):

        pid = ID.match(pr).group(1)
        org = PR_ORGANISM_NAME.match(pr).group(1)
        upids = g.group(1) if (g := PR_IDENTIFIER.match(pr)) else ''
        # Dealing with cases of multiple or empty identifiers
        upids = [i for i in upids.split(' AND ') if i]
        refs = REFERENCE.match(pr).group(1)

        proc[pid]['EC'] = eid
        proc[pid]['UniProt'].update(upids)
        proc[pid]['#'] = pid
        proc[pid]['Organism'].add(org)
        proc[pid]['Refs'].update(refs.split(','))

    for ac in record.get('AC', []):

        compound = AC_COMPOUND.match(ac).group(1)
        aux = x.group(1) if (x := AC_DETAILS.match(ac)) else ''

        # Making sure all IDs are accounted for (not all described inside
        # parentheses)
        initial_pids = set(ID.match(ac).group(1).split(','))
        initial_refs = set(REFERENCE.match(ac).group(1).split(','))

        for entry in aux.split('; '):

            if not entry:

                continue

            pids = ID.match(entry).group(1).split(',')
            refs = REFERENCE.match(entry).group(1).split(',')

            for pid in pids:

                proc[pid]['Activator'].add(compound)
                proc[pid]['Refs'].update(refs)

                initial_pids -= {pid}

            initial_refs -= set(refs)

        # Adding remaining annotations for non-described entries
        for pid in initial_pids:

            proc[pid]['Activator'].add(compound)
            proc[pid]['Refs'].update(initial_refs) # Assuming a bit here

    for inh in record.get('IN', []):

        compound = AC_COMPOUND.match(inh).group(1)
        aux = x.group(1) if (x := AC_DETAILS.match(inh)) else ''

        # Making sure all IDs are accounted for (not all described inside
        # parentheses)
        initial_pids = set(ID.match(inh).group(1).split(','))
        initial_refs = set(REFERENCE.match(inh).group(1).split(','))

        for entry in aux.split('; '):

            if not entry:

                continue

            pids = ID.match(entry).group(1).split(',')
            refs = REFERENCE.match(entry).group(1).split(',')

            for pid in pids:

                proc[pid]['Inhibitor'].add(compound)
                proc[pid]['Refs'].update(refs)

                initial_pids -= {pid}

            initial_refs -= set(refs)

        # Adding remaining annotations for non-described entries
        for pid in initial_pids:

            proc[pid]['Inhibitor'].add(compound)
            proc[pid]['Refs'].update(initial_refs) # Assuming a bit here

    return proc


def parser(opener, **kwargs):

    keys = sorted(opener.result)
    # Serialized raw string
    result = opener.result[keys[0]].read().decode('utf-8')

    entries = [
        [
            split
            for line in entry.split('\n')
            if (split := line.split('\t', maxsplit=1))[0] in RECORDS_ENABLED
        ]
        # NOTE: Continuation lines start with \t character -> merged from
        # serialized string by replacing \n\t with empty string before splitting
        # the lines. Each entry spans multiple lines, separated by "///"
        for entry in result.replace('\n\t', '').split('///')[1:] # Ignore first
    ]

    records = [process_record(pe) for e in entries if (pe := process_entry(e))]

    yield from records


download = Download(
    url=URL,
    filename='brenda.txt.tar.gz',
    subfolder='brenda',
    large=True,
    ext='tar.gz',
    default_mode='rb',
    download_kwargs={
        'post': True,
        'query': {'dlfile': 'dl-textfile', 'accept-license': '1'},
    }
)

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

#schema = EntityBuilder(
#    entity_type=EntityTypeCv.PROTEIN,
#    identifiers=IdentifiersBuilder(
#        CV(term=IdentifierNamespaceCv.UNIPROT, value=f(#Extract uniprot)),
#    ),
#    annotations=AnnotationsBuilder(
#        CV(term=IdentifierNamespaceCv.EC, value=f('ID')),
#        CV(
#            term=EntityTypeCv.ORGANISM,
#            value=f('PR', delimiter='||', extract='pr')
#        ),
#        # Put EC as ontology
#        #CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('Gene Ontology IDs', delimiter=';')),
#    ),
#)

# ================================= RESOURCE ===================================

#resource = Resource(
#    config=config,
#    associations=Dataset(
#        download=download,
#        mapper=schema,
#        raw_parser=parser,
#    ),
#)

# ================================= REFERENCE ==================================
