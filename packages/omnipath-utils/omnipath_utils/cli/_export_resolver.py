"""``omnipath-utils export-resolver`` subcommand (Milestone I).

Emits canonical-target parquet projecting the symmetric ``id_mapping`` model for
an entity family, honouring a shared policy file. See
contracts/export-resolver-cli.md.
"""

from __future__ import annotations

import argparse
import logging
import os


def export_resolver_cmd(args: list[str]) -> None:
    parser = argparse.ArgumentParser(
        prog='omnipath-utils export-resolver',
        description='Export canonical-target resolver parquet for an entity family.',
    )
    parser.add_argument(
        '--entity-family',
        required=True,
        choices=['protein', 'chemical'],
        help='Entity family to project.',
    )
    parser.add_argument(
        '--policy',
        required=True,
        help='Path to the resolver policy YAML (entity_family → key-type rules).',
    )
    parser.add_argument(
        '--canonical-type',
        default=None,
        help="Canonical target id type (default: uniprot for protein, "
        "'Standard InChI Key' for chemical).",
    )
    parser.add_argument('--output', required=True, help='Output directory.')
    parser.add_argument('--db-url', default=None, help='Override OMNIPATH_UTILS_DB_URL.')
    parser.add_argument(
        '--organisms',
        nargs='+',
        type=int,
        default=None,
        help='Taxa for the protein family (default 9606); ignored for chemicals.',
    )
    parser.add_argument('--max-records', type=int, default=None)
    parser.add_argument(
        '--export-id',
        default=None,
        help='Stable id recorded in manifest.json (default: UTC timestamp).',
    )
    parser.add_argument('-v', '--verbose', action='store_true')

    opts = parser.parse_args(args)
    logging.basicConfig(level=logging.DEBUG if opts.verbose else logging.INFO)

    from omnipath_utils.db._resolver_export import (
        MANIFEST_FILENAME,
        export_resolver,
    )

    stats = export_resolver(
        family=opts.entity_family,
        policy_path=opts.policy,
        output_dir=opts.output,
        canonical_type=opts.canonical_type,
        db_url=opts.db_url,
        taxa=opts.organisms,
        max_records=opts.max_records,
        export_id=opts.export_id,
    )
    print(
        f'[export-resolver] family={stats.family} canonical={stats.canonical_type} '
        f'files={len(stats.files)} rows={stats.rows}'
    )
    for path in stats.files:
        print(f'  {path}')
    print(f'  manifest: {os.path.join(opts.output, MANIFEST_FILENAME)}')
