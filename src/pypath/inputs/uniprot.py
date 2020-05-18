#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems

import re
import time
import datetime
import collections
import itertools
import timeloop
timeloop.app.logging.disable(level = 9999)

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.session as session_mod
import pypath.share.settings as settings

_logger = session_mod.Logger(name = 'uniprot_input')

db = {}
_cleanup_period = settings.get('mapper_cleanup_interval')
_lifetime = 300
_last_used = {}


_redatasheet = re.compile(r'([A-Z\s]{2})\s*([^\n\r]+)[\n\r]+')

# regex for matching UniProt AC format
# from https://www.uniprot.org/help/accession_numbers
reac = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
)
_rename = re.compile(r'Name=([\w\(\)-]+)\W')
_retaxid = re.compile(r'=(\d+)[^\d]')


def _all_uniprots(organism = 9606, swissprot = None):

    swissprot = _swissprot_param(swissprot)
    rev = '' if swissprot is None else ' AND reviewed: %s' % swissprot
    url = urls.urls['uniprot_basic']['url']
    get = {
        'query': 'organism:%s%s' % (str(organism), rev),
        'format': 'tab',
        'columns': 'id',
    }
    c = curl.Curl(url, get = get, silent = False)
    data = c.result

    return {
        l.strip() for l in data.split('\n')[1:] if l.strip()
    }


def all_uniprots(organism = 9606, swissprot = None):

    return get_db(organism = organism, swissprot = swissprot)


def all_swissprots(organism = 9606):

    return get_db(organism = organism, swissprot = True)


def all_trembls(organism = 9606):

    return get_db(organism = organism, swissprot = False)


def init_db(organism = 9606, swissprot = None):

    swissprot = _swissprot_param(swissprot)
    _logger._log(
        'Loading list of all UniProt IDs for '
        'organism `%u` (only SwissProt: %s).' % (
            organism,
            str(swissprot),
        )
    )

    key = (organism, swissprot)

    globals()['db'][key] = _all_uniprots(
        organism = organism,
        swissprot = swissprot,
    )
    globals()['_last_used'][key] = time.time()


def get_db(organism = 9606, swissprot = None):

    swissprot = _swissprot_param(swissprot)
    key = (organism, swissprot)

    if key not in globals()['db']:

        init_db(organism = organism, swissprot = swissprot)

    globals()['_last_used'][key] = time.time()

    return globals()['db'][key]


def _swissprot_param(swissprot):

    return (
        'yes'
            if swissprot in {'yes', 'YES', True} else
        'no'
            if swissprot in {'no', 'NO', False} else
        None
    )


def valid_uniprot(name):
    """
    Checks if ``name`` fits the format requirements for UniProt accession
    numbers.
    """

    return bool(reac.match(name))


def is_uniprot(name, organism = 9606, swissprot = None):
    """
    Tells if ``name`` is a UniProt ID of ``organism``.
    If ``swissprot`` is None then both SwissProt and TrEMBL IDs will be
    considered.
    """

    return name in get_db(organism = organism, swissprot = swissprot)


def is_swissprot(name, organism = 9606):
    """
    Tells if ``name`` is a SwissProt ID of ``organism``.
    For TrEMBL IDs returns False.
    """

    return is_uniprot(name, organism = organism, swissprot = True)


def is_trembl(name, organism = 9606):
    """
    Tells if ``name`` is a TrEMBL ID of ``organism``.
    For SwissProt IDs returns False.
    """

    return is_uniprot(name, organism = organism, swissprot = False)


_cleanup_timeloop = timeloop.Timeloop()

@_cleanup_timeloop.job(
    interval = datetime.timedelta(
        seconds = _cleanup_period
    )
)
def _cleanup():

    keys = list(globals()['db'].keys())

    for key in keys:

        if time.time() - globals()['_last_used'][key] > _lifetime:

            _remove(key)

_cleanup_timeloop.start(block = False)


def _remove(key):

    if key in globals()['db']:

        _logger._log(
            'Removing UniProt ID list for '
            'organism `%u` (SwissProt: %s)' % (
                key[0],
                str(key[1]),
            )
        )
        del globals()['db'][key]

    if key in globals()['_last_used']:

        del globals()['_last_used'][key]


def protein_datasheet(identifier):

    url = urls.urls['uniprot_basic']['datasheet'] % identifier.strip()

    datasheet =  _protein_datasheet(url)

    if not datasheet:

        _logger._log(
            'UniProt ID `%s` returns empty response, it might be and an old '
            'ID which has been deleted from the database. Attempting to '
            'find its history and retrieve either an archived version or '
            'the find the new ID which replaced this one.' % identifier
        )
        return uniprot_history_recent_datasheet(identifier)

    else:

        return datasheet


def deleted_uniprot_genesymbol(identifier):
    """
    Retrieves the archived datasheet for a deleted UniProt ID and returns
    the Gene Symbol and the NCBI Taxonomy ID from the datasheet.
    """

    datasheet = uniprot_history_recent_datasheet(identifier)
    genesymbol = None
    ncbi_tax_id = None

    for tag, line in datasheet:

        if tag == 'GN':

            m = _rename.search(line.strip())

            if m:

                genesymbol = m.groups()[0]

        if tag == 'OX':

            ncbi_tax_id = int(_retaxid.search(line).groups()[0])
            break

    return genesymbol, ncbi_tax_id


def _protein_datasheet(url):

    cache = True

    for a in range(3):

        c = curl.Curl(
            url,
            silent = True,
            large = False,
            cache = cache,
            connect_timeout = (
                settings.get('uniprot_datasheet_connect_timeout')
            ),
            timeout = settings.get('uniprot_datasheet_timeout'),
        )

        if not c.result or c.result.startswith('<!DOCTYPE'):

            cache = False

        else:

            break

    if not c.result:

        _logger._log(
            'Could not retrieve UniProt datasheet by URL `%s`.' % url
        )

    return _redatasheet.findall(c.result) if c.result else []


def uniprot_history_recent_datasheet(identifier):

    recent_version = uniprot_recent_version(identifier)

    if recent_version:

        if recent_version.replaced_by:

            new = recent_version.replaced_by.split(';')[0]
            url = urls.urls['uniprot_basic']['datasheet'] % new
            _logger._log(
                'UniProt ID `%s` is obsolete, has been replaced by '
                '`%s`: `%s`.' % (
                    identifier,
                    new,
                    url,
                )
            )
            return protein_datasheet(new)

        else:

            version = int(recent_version.entry_version)
            url = '%s?version=%u' % (
                urls.urls['uniprot_basic']['datasheet'] % identifier,
                version,
            )
            _logger._log(
                'UniProt ID `%s` is obsolete, downloading archived '
                'version %u: `%s`.' % (
                    identifier,
                    version,
                    url,
                )
            )
            c = curl.Curl(url, silent = True, large = False)
            return _protein_datasheet(url)

    return []


UniprotRecordHistory = collections.namedtuple(
    'UniprotRecordHistory',
    [
        'entry_version',
        'sequence_version',
        'entry_name',
        'database',
        'number',
        'date',
        'replaces',
        'replaced_by',
    ],
)


def uniprot_history(identifier):
    """
    Retrieves the history of a record.
    Returns a generator iterating over the history from most recent to the
    oldest.
    """

    if valid_uniprot(identifier):

        url_history = urls.urls['uniprot_basic']['history'] % identifier
        c_history = curl.Curl(
            url_history,
            silent = True,
            large = True,
        )

        if c_history.result:

            line0 = next(c_history.result)

            if not line0.startswith('<!DOCTYPE'):

                for line in c_history.result:

                    if line:

                        yield UniprotRecordHistory(
                            *(
                                field.strip() for field in line.split('\t')
                            )
                        )


def uniprot_recent_version(identifier):

    for version in uniprot_history(identifier):

        if (
            (
                version.entry_version != '0' and
                version.entry_name != 'null'
            ) or version.replaced_by
        ):

            return version


def uniprot_deleted(confirm = True):

    return swissprot_deleted() | trembl_deleted(confirm = confirm)


def _uniprot_deleted(swissprot = True, confirm = True):

    if not swissprot and confirm:

        resp = input(
            'Loading the list of deleted TrEMBL IDs requires '
            '>5GB memory. Do you want to proceed [y/n] '
        )

        if not resp or resp[0].lower() != 'y':

            return set()

    key = 'deleted_%s' % ('sp' if swissprot else 'tr')
    url = urls.urls['uniprot_basic'][key]
    c = curl.Curl(url, silent = False, large = True)

    result = set()

    for line in c.result:

        m = reac.match(line.strip())

        if m:

            result.add(m.groups()[0])

    return result


def swissprot_deleted():

    return _uniprot_deleted(swissprot = True)


def trembl_deleted(confirm = True):

    return _uniprot_deleted(swissprot = False, confirm = True)


def get_uniprot_sec(organism = 9606):
    """
    Downloads and processes the mapping between secondary and
    primary UniProt IDs.

    Yields pairs of secondary and primary UniProt IDs.

    :param int organism:
        NCBI Taxonomy ID of the organism.
    """

    if organism is not None:
        proteome = all_uniprots(organism=organism)
        proteome = set(proteome)

    sec_pri = []
    url = urls.urls['uniprot_sec']['url']
    c = curl.Curl(url, silent = False, large = True, timeout = 2400)

    for line in filter(
        lambda line:
            len(line) == 2 and (organism is None or line[1] in proteome),
            map(
                lambda i:
                    i[1].split(),
                filter(
                    lambda i: i[0] >= 30,
                    enumerate(c.result)
                )
            )
        ):

        yield line


_uniprot_fields = {
    'function': 'comment(FUNCTION)',
    'activity_regulation': 'comment(ACTIVITY REGULATION)',
    'tissue_specificity': 'comment(TISSUE SPECIFICITY)',
    'developmental_stage': 'comment(DEVELOPMENTAL STAGE)',
    'induction': 'comment(INDUCTION)',
    'intramembrane': 'feature(INTRAMEMBRANE)',
    'signal_peptide': 'feature(SIGNAL)',
    'subcellular_location': 'comment(SUBCELLULAR LOCATION)',
    'transmembrane': 'feature(TRANSMEMBRANE)',
    'comment': 'comment(MISCELLANEOUS)',
    'topological_domain': 'feature(TOPOLOGICAL DOMAIN)',
    'family': 'families',
    'interactor': 'interactor',
    'keywords': 'keywords',
}


def uniprot_data(field, organism = 9606, reviewed = True):
    """
    Retrieves a field from UniProt for all proteins of one organism, by
    default only the reviewed (SwissProt) proteins.
    For the available fields refer to the ``_uniprot_fields`` attribute of
    this module or the UniProt website.
    """

    rev = (
        ' AND reviewed: yes'
            if reviewed == True or reviewed == 'yes' else
        ' AND reviewed: no'
        if reviewed == False or reviewed == 'no' else
        ''
    )
    _field = _uniprot_fields[field] if field in _uniprot_fields else field
    url = urls.urls['uniprot_basic']['url']
    get = {
        'query': 'organism:%s%s' % (str(organism), rev),
        'format': 'tab',
        'columns': 'id,%s' % _field,
        'compress': 'yes',
    }

    c = curl.Curl(url, get = get, silent = False, large = True, compr = 'gz')
    _ = next(c.result)

    return dict(
        id_value
        for id_value in
        (
            line.strip('\n\r').split('\t')
            for line in c.result if line.strip('\n\r')
        )
        if id_value[1]
    )


def uniprot_preprocess(field, organism = 9606, reviewed = True):

    relabel = re.compile(r'[A-Z\s]+:\s')
    reisoform = re.compile(r'\[[-\w\s]+\]:?\s?')
    retermsep = re.compile(r'\s?[\.,]\s?')
    reref = re.compile(r'\{[-\w :\|,\.]*\}')

    result = collections.defaultdict(set)

    data = uniprot_data(
        field = field,
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, raw in iteritems(data):

        raw = raw.split('Note=')[0]
        raw = relabel.sub('', raw)
        raw = reref.sub('', raw)
        raw = reisoform.sub('', raw)
        raw = retermsep.split(raw)

        for item in raw:

            if item.startswith('Note'):

                continue

            item = item.split('{')[0]
            elements = tuple(
                it0
                for it0 in
                (
                    common.upper0(it.strip(' .;,'))
                    for it in item.split(';')
                )
                if it0
            )

            if elements:

                result[uniprot].add(elements)

    return result


def uniprot_locations(organism = 9606, reviewed = True):


    UniprotLocation = collections.namedtuple(
        'UniprotLocation',
        [
            'location',
            'features',
        ],
    )


    result = collections.defaultdict(set)

    data = uniprot_preprocess(
        field = 'subcellular_location',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, locations in iteritems(data):

        for location in locations:

            result[uniprot].add(
                UniprotLocation(
                    location = location[0],
                    features = location[1:] or None,
                )
            )

    return result


def uniprot_keywords(organism = 9606, reviewed = True):

    UniprotKeyword = collections.namedtuple(
        'UniprotKeyword',
        [
            'keyword',
        ],
    )


    result = collections.defaultdict(set)

    data = uniprot_data(
        field = 'keywords',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, keywords in iteritems(data):

        for keyword in keywords.split(';'):

            result[uniprot].add(
                UniprotKeyword(
                    keyword = keyword.strip(),
                )
            )

    return result


def uniprot_families(organism = 9606, reviewed = True):

    refamily = re.compile(r'(.+) (?:super)?family(?:, (.*) subfamily)?')


    UniprotFamily = collections.namedtuple(
        'UniprotFamily',
        [
            'family',
            'subfamily',
        ],
    )


    result = collections.defaultdict(set)

    data = uniprot_data(
        field = 'family',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, family in iteritems(data):

        if not family:

            continue

        family, subfamily = refamily.search(family).groups()

        result[uniprot].add(
            UniprotFamily(
                family = family,
                subfamily = subfamily,
            )
        )

    return result


def uniprot_topology(organism = 9606, reviewed = True):

    retopo = re.compile(r'TOPO_DOM (\d+)\.\.(\d+);\s+/note="(\w+)"')
    retm = re.compile(r'(TRANSMEM|INTRAMEM) (\d+)\.\.(\d+);')


    UniprotTopology = collections.namedtuple(
        'UniprotTopology',
        [
            'topology',
            'start',
            'end',
        ],
    )


    result = collections.defaultdict(set)

    transmem = uniprot_data(
        field = 'transmembrane',
        organism = organism,
        reviewed = reviewed,
    )

    intramem = uniprot_data(
        field = 'intramembrane',
        organism = organism,
        reviewed = reviewed,
    )

    signal = uniprot_data(
        field = 'signal_peptide',
        organism = organism,
        reviewed = reviewed,
    )

    data = uniprot_data(
        field = 'topological_domain',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, topo in iteritems(data):

        for topo_dom in retopo.findall(topo):

            start, end, topology = topo_dom
            start = int(start)
            end = int(end)

            result[uniprot].add(
                UniprotTopology(
                    topology = topology,
                    start = start,
                    end = end,
                )
            )

    for uniprot, tm in itertools.chain(
        iteritems(transmem),
        iteritems(intramem),
        iteritems(signal),
    ):

        for mem, start, end in retm.findall(tm):

            topology = (
                '%s%s' % (
                    mem.capitalize(),
                    'brane' if mem.endswith('MEM') else ''
                )
            )
            start = int(start)
            end = int(end)

            result[uniprot].add(
                UniprotTopology(
                    topology = topology,
                    start = start,
                    end = end,
                )
            )

    return result


def uniprot_tissues(organism = 9606, reviewed = True):

    reref = re.compile(r'\s?\{.*\}\s?')
    resep = re.compile(
        r',?(?:'
            r' in almost all |'
            r' but also in |'
            r' but also at |'
            r' within the |'
            r', in |'
            r' in |'
            r' but |'
            r', and |'
            r' and |'
            r' such as |'
            r' \(both |'
            r' as well as |'
            r' as |'
            r' or |'
            r' at the |'
            r' at |'
            r' including |'
            r' during |'
            r' especially |'
            r' to |'
            r' into |'
            r' = |'
            r' > |'
            r'; |'
            r', '
        r')(?=[^\d])'
    )
    relabel = re.compile(r'^TISSUE SPECIFICITY: ')
    repubmed = re.compile(r'\(?PubMed:?\d+\)?')
    respeci = re.compile(r'(\w+)[-\s]specific')
    rethe = re.compile(
        r'\s?(?:'
           r'[Tt]he |'
           r'[Ii]n |'
           r'[Ss]ome|'
           r'[Ii]n the|'
           r'[Ww]ithin the|'
           r'[Ww]ithin|'
           r'[Ii]nto|'
           r'[Ww]ith only|'
           r'[Ww]ith the|'
           r'[Ww]ith an|'
           r'[Ww]ith |'
           r'[Ii]s |'
           r'[Mm]any  |'
           r'[Aa] variety of '
           r'[Aa] |'
           r'[Ii]t |'
           r'[Tt]o |'
           r'[Oo]n |'
           r'[Oo]f |'
           r'[Tt]hose |'
           r'[Ff]rom |'
           r'[Aa]lso|'
           r'[Bb]y |'
           r'[Pp]articularly|'
           r'[Pp]articular|'
           r'[Pp]atients|'
           r'[Aa]n |'
           r'\'|'
           r':|'
           r'/'
        r')?(.*)'
    )
    reand = re.compile(r'(?: and| of| from| or| than)$')
    replevel = re.compile(r'\(at \w+ levels?\)')
    reiso = re.compile(r'[Ii]soform \w+')
    reindef = re.compile(
        r'\w'
        r'(?:'
           r'ifferent parts of |'
           r'ariety of tissues |'
           r' variety of tissues |'
           r' number of |'
           r'everal regions of '
        r')'
    )

    level_kw = (
        ('low', 'low'),
        ('weak', 'low'),
        ('lesser extent', 'low'),
        ('minimal level', 'low'),
        ('decrease', 'low'),
        ('moderate', 'low'),
        ('barely', 'low'),
        ('minor level', 'low'),
        ('reduced', 'low'),
        ('lesser', 'low'),
        ('down-regulated', 'low'),
        ('high', 'high'),
        ('elevated', 'high'),
        ('strong', 'high'),
        ('prominent', 'high'),
        ('greatest level', 'high'),
        ('concentrated', 'high'),
        ('predominant', 'high'),
        ('increase', 'high'),
        ('enrich', 'high'),
        ('abundant', 'high'),
        ('primarily', 'high'),
        ('induced', 'high'),
        ('up-regulated', 'high'),
        ('up regulated', 'high'),
        ('expression is restricted', 'high'),
        ('amplified', 'high'),
        ('basal l', 'basal'),
        ('not detected', 'none'),
        ('absent', 'none'),
        ('expressed', 'undefined'),
        ('detect', 'undefined'),
        ('found', 'undefined'),
        ('present', 'undefined'),
        ('expression', 'undefined'),
        ('localized', 'undefined'),
        ('produced', 'undefined'),
        ('confined', 'undefined'),
        ('transcribed', 'undefined'),
        ('xpressed', 'undefined'),
        ('synthesized', 'undefined'),
        ('secreted', 'undefined'),
        ('seen', 'undefined'),
        ('prevalent', 'undefined'),
        ('released', 'undefined'),
        ('appears', 'undefined'),
        ('varying levels', 'undefined'),
        ('various levels', 'undefined'),
        ('identified', 'undefined'),
        ('observed', 'undefined'),
        ('occurs', 'undefined'),
    )

    wide_kw = (
        ('widely', 'wide'),
        ('wide tissue distribution', 'wide'),
        ('wide range of tissues', 'wide'),
        ('wide range of adult tissues', 'wide'),
        ('wide range of cells', 'wide'),
        ('wide variety of normal adult tissues', 'wide'),
        ('widespread', 'wide'),
        ('ubiquitous', 'ubiquitous'),
        ('variety of tissues', 'wide'),
        ('many tissues', 'wide'),
        ('many organs', 'wide'),
        ('various organs', 'wide'),
        ('various tissues', 'wide'),
    )

    tissue_exclude = {
        'Adult',
        'All',
        'Apparently not',
        'Areas',
        'Are likely',
        'Both',
        'By contrast',
        'Normal cells',
        'Not only',
        'A',
        '[]: Localized',
        'Early',
        'Change from a quiescent',
        'Central',
        'Beta',
        'This layer',
        'With little',
        'Preferential occurrence',
        'Stage III',
        'Take up',
        'Hardly',
        'Only seen',
        'Prevalent',
        'Inner segment',
        'Memory',
        'Many fetal',
        'Tissues',
        '0 kb',
        '9 kb',
        'A 2',
        'A 3',
        'A 5',
        'A 6',
        '1-7',
        '1b-1',
        '2 is widely',
        '8 and 4',
        'Often amplified',
        'Other',
        'Others',
        'Those',
        'Tissues examined',
        'Tissues with',
        'Tissues (e)',
        'Probably shed',
        'Reports that',
        'Primitive',
        'Prolactin',
        'Overlap',
        'A smaller 0',
        'A smaller form',
        'A smaltissues',
        'Different levels',
        'Different amounts',
        'Disappears',
        'Digestion',
        'Very similar',
        'Vivo',
        'Contrary',
        'Contrast',
        'Not',
        'Not all',
        'Has it',
        'Has little',
        'All stages',
        'Soon',
        'Specific',
        'Stage',
        'Stage I',
        'Stage II',
        'Stages II',
        'Ends',
        'A minor degree',
        'A much smaller extent',
        'Lost',
        'Varies',
        'Various',
        'Mostly restricted',
        'Mostly',
        'Most probably',
        'Much more stable',
        'Naive',
        'Neither',
        'Nor',
        'None',
    }

    exclude_startswith = (
        'Were',
        'Where',
        'Which',
        'While',
        'When',
        'There',
        'Their',
        'Then',
        'These',
        'Level',
        'This',
        'Almost',
        'If',
        'Control',
        'Be ',
        'Although',
        'Than',
        'Addition',
    )

    exclude_in = (
        'kb transcript',
        'compared',
        'soform',
        'concentration of'
    )


    UniprotTissue = collections.namedtuple(
        'UniprotTissue',
        [
            'tissue',
            'level',
        ],
    )


    data = uniprot_data(
        'tissue_specificity',
        organism = organism,
        reviewed = reviewed,
    )

    result = collections.defaultdict(set)

    for uniprot, raw in iteritems(data):

        raw = relabel.sub('', raw)
        raw = reref.sub('', raw)
        raw = replevel.sub('', raw)
        raw = reiso.sub('', raw)
        raw = repubmed.sub('', raw)
        raw = reindef.sub('', raw)
        raw = raw.replace('adult and fetal', '')

        raw = raw.split('.')

        for phrase in raw:

            tokens = tuple(resep.split(phrase))
            level = None

            for token in tokens:

                level_token = False
                wide_token = False
                tissue = None

                token_lower = token.lower()

                for kw, lev in level_kw:

                    if kw in token_lower:

                        level = lev
                        level_token = True
                        break

                if level_token:

                    for kw, wide in wide_kw:

                        if kw in token_lower:

                            tissue = wide
                            wide_token = True
                            break

                if not level_token or wide_token:

                    if not wide_token:

                        specific = respeci.search(token)

                        tissue = (
                            specific.groups()[0].lower()
                                if specific else
                            token
                        )

                        if specific and not level:

                            level = 'high'

                    if tissue.strip():

                        if any(e in tissue for e in exclude_in):

                            continue

                        tissue = rethe.match(tissue).groups()[0]
                        tissue = rethe.match(tissue).groups()[0]
                        tissue = rethe.match(tissue).groups()[0]

                        if tissue.endswith('+'):

                            tissue = '%s cells' % tissue

                        tissue = tissue.strip(')(.,;- ')

                        if '(' in tissue and ')' not in tissue:

                            tissue = '%s)' % tissue

                        tissue = reand.sub('', tissue)
                        tissue = common.upper0(tissue)
                        tissue = tissue.replace('  ', ' ')

                        if any(
                            tissue.startswith(e)
                            for e in exclude_startswith
                        ):

                            continue

                        if tissue in tissue_exclude or len(tissue) < 3:

                            continue

                        result[uniprot].add(
                            UniprotTissue(
                                tissue = tissue,
                                level = level or 'undefined',
                            )
                        )

    return result


def uniprot_taxonomy():
    """
    Returns a dictionary with SwissProt IDs as keys and sets of various taxon
    names as values.
    """
    
    rename = re.compile(r'\(?(\w[\w\s\',/\.-]+\w)\)?')
    reac = re.compile(r'\s*\w+\s+\(([A-Z\d]+)\)\s*,')

    url = urls.urls['uniprot_basic']['speindex']
    c = curl.Curl(url, large = True, silent = False)

    result = collections.defaultdict(set)
    
    for line in c.result:
        
        if line[0] != ' ':
            
            names = set(rename.findall(line))
            
        else:
            
            for ac in reac.findall(line):
                
                result[ac].update(names)
    
    return result


Taxon = collections.namedtuple(
    'Taxon',
    [
        'ncbi_id',
        'latin',
        'english',
        'latin_synonym',
    ]
)
Taxon.__new__.__defaults__ = (None, None)


def uniprot_ncbi_taxids():
    
    url = urls.urls['uniprot_basic']['taxids']
    c = curl.Curl(url, large = True, silent = False, compr = 'gz')
    
    _ = next(c.result)
    
    result = {}
    
    for line in c.result:
        
        line = line.split('\t')
        
        if line[0].isdigit() and len(line) > 3:
            
            taxid = int(line[0])
            
            result[taxid] = Taxon(
                ncbi_id = taxid,
                latin = line[2],
                english = line[3],
            )
    
    return result


def uniprot_ncbi_taxids_2():
    
    reline = re.compile(
        r'(?:([A-Z\d]+)\s+)?' # code
        r'(?:([A-Z]))?\s+' # kingdom
        r'(?:(\d+): )?' # NCBI Taxonomy ID
        r'([A-Z])=' # name type
        r'([ \w\(\),/\.\'-]+)[\n\r\s]*' # the name
    )
    
    url = urls.urls['uniprot_basic']['speclist']
    c = curl.Curl(url, large = True, silent = False)
    
    result = {}
    entry = {}
    
    for line in c.result:
        
        m = reline.match(line)
        
        if m:
            
            _code, _kingdom, _taxid, _name_type, _name = m.groups()
            
            if _taxid:
                
                if entry and 'ncbi_id' in entry:
                    
                    result[entry['ncbi_id']] = Taxon(**entry)
                
                entry = {}
                entry['ncbi_id'] = int(_taxid)
            
            if _name_type == 'N':
                
                entry['latin'] = _name
                
            elif _name_type == 'C':
                
                entry['english'] = _name
                
            elif _name_type == 'S':
                
                entry['latin_synonym'] = _name
    
    if entry and 'ncbi_id' in entry:
        
        result[entry['ncbi_id']] = Taxon(**entry)
    
    return result
