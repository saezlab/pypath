def rolland_hi_ii_14():
    """
    Loads the HI-II-14 unbiased interactome from the large scale screening
    of from Rolland 2014.
    Returns list of interactions.
    """
    url = urls.urls['hiii14']['url']
    c = curl.Curl(url, silent = False, large = True)
    xlsname = c.fileobj.name
    c.fileobj.close()
    tbl = read_xls(xlsname, sheet = '2G')
    
    for row in tbl[1:]:
        
        yield [c.split('.')[0] for c in row]


def vidal_hi_iii_old(fname):
    """
    Loads the HI-III  unbiased interactome from preliminary data of
    the next large scale screening of Vidal Lab.

    The data is accessible here:
        http://interactome.dfci.harvard.edu/H_sapiens/dload_trk.php
    You need to register and accept the license terms.

    Returns list of interactions.
    """

    f = curl.FileOpener(fname)
    return \
        list(
            map(
                lambda l:
                    l.strip().split('\t'),
                f.result
            )
        )[1:]


def hi_iii():
    """
    Loads the unbiased human interactome version III (HI-III).
    This is an unpublished data and its use is limited.
    Please check the conditions and licensing terms carefully at
    http://interactome.baderlab.org.
    """
    
    HiiiiInteraction = collections.namedtuple(
        'HiiiiInteraction',
        [
            'id_a',
            'id_b',
            'isoform_a',
            'isoform_b',
            'screens',
            'score',
        ]
    )
    
    
    rescore = re.compile(r'author score: ([\d\.]+)')
    rescreens = re.compile(r'Found in screens ([\d,]+)')
    
    url = urls.urls['hid']['hi-iii']
    post_data = {
        'form[request_dataset]': '2',
        'form[request_file_format]': 'psi',
    }
    c = curl.Curl(url, silent = False, large = True, post = post_data)
    
    for row in c.result:
        
        if not row.strip():
            
            continue
        
        id_a, id_b, rest = row.split(' ', maxsplit = 2)
        id_a, isoform_a = id_a.split('-') if '-' in id_a else (id_a, 1)
        id_b, isoform_b = id_b.split('-') if '-' in id_b else (id_b, 1)
        
        sc = rescore.search(rest)
        score = float(sc.groups()[0]) if sc else None
        screens = tuple(
            int(i) for i in rescreens.search(rest).groups()[0].split(',')
        )
        
        yield HiiiiInteraction(
            id_a = id_a[10:],
            id_b = id_b[10:],
            isoform_a = int(isoform_a),
            isoform_b = int(isoform_b),
            screens = screens,
            score = score,
        )
