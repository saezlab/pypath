#!/usr/bin/env python

import os
import gzip
import requests

targets = (
    ('annotations', 'PROGENy'),
    ('annotations', 'MSigDB'),
    ('annotations', 'CytoSig'),
    ('annotations', 'PanglaoDB'),
    ('annotations', 'HPA'),
    ('interactions', 'collectri'),
    ('interactions', 'dorothea'),
)

for query, resource in targets:

    organisms = (9606,) if query == 'annotations' else (9606, 10090, 10116)
    resource_param = 'datasets' if query == 'interactions' else 'resources'

    for organism in organisms:

        outfile = os.path.expanduser(
            os.path.join(
                '~',
                'static',
                'resources',
                f'{query}_{resource}_{organism}.tsv.gz',
            )
        )

        dorothea = (
            ',dorothea_level&dorothea_levels=A,B,C,D'
                if resource == 'dorothea' else
            ''
        )
        query_param = (
            f'&organisms={organism}'
            '&genesymbols=yes'
            f'&fields=sources,references,curation_effort,evidences{dorothea}'
        ) if query == 'interactions' else ''

        with gzip.open(outfile, 'wb') as target:

            url = (
                f'https://omnipathdb.org/{query}?'
                f'{resource_param}={resource}{query_param}'
            )
            print(f'Saving `{url}` to `{outfile}`.')
            req = requests.get(url, stream = True)

            for chunk in req.iter_content(chunk_size=1024):
                target.write(chunk)
