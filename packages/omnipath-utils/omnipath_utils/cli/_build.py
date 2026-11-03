"""Build command for omnipath-utils database."""

from __future__ import annotations

import logging
import argparse


def build_cmd(args: list[str]):
    """Execute the build command."""
    parser = argparse.ArgumentParser(
        description='Build omnipath-utils database'
    )
    parser.add_argument(
        '--db-url', default=None, help='PostgreSQL connection URL'
    )
    parser.add_argument(
        '--organisms',
        nargs='+',
        type=int,
        default=[9606],
        help='NCBI Taxonomy IDs to build for (default: 9606)',
    )
    parser.add_argument(
        '--ref-only',
        action='store_true',
        help='Only build reference tables (id_types, organisms, backends)',
    )
    parser.add_argument(
        '--ftp',
        action='store_true',
        help='Build from full UniProt FTP idmapping.dat.gz (all organisms, ~18GB)',
    )
    parser.add_argument(
        '--metabolites',
        action='store_true',
        help=(
            'Load ONLY the metabolite ID mappings (UniChem, RaMP, MetaNetX, '
            'BiGG, and structure->InChIKey incl. PubChem). Additive: does not '
            'touch the protein or FTP idmapping tables. Honours '
            '--pubchem-max-records (omit it for the full PubChem table).'
        ),
    )
    parser.add_argument(
        '--preset',
        choices=['minimal', 'standard', 'model', 'full'],
        default=None,
        help='Build preset',
    )
    parser.add_argument(
        '--scope',
        default=None,
        help=(
            'Build scope (007): a nested preset — only-human, core-model, '
            'extended-model, model-organisms, complete — or an explicit '
            'comma-separated list of NCBI taxonomy ids / organism names '
            '(e.g. "human,pig,chimpanzee" or "9606,9823,9598"). '
            'Overrides --organisms; drives build_all.'
        ),
    )
    parser.add_argument(
        '--parquet-dir',
        default=None,
        help='Directory for Parquet file exports',
    )
    parser.add_argument(
        '--list-presets',
        action='store_true',
        help='List available presets and exit',
    )
    parser.add_argument(
        '--max-records',
        type=int,
        default=None,
        help=(
            'Cap each mapping table at this many rows (for fast test builds). '
            'Applied at load time for inputs_v2 backends and as a write cap '
            'for all others. Default: no limit (full build).'
        ),
    )
    parser.add_argument(
        '--pubchem-max-records',
        type=int,
        default=None,
        help=(
            'Cap ONLY the PubChem -> InChIKey table at this many rows, while '
            'every other resource loads in full. Use this to run an otherwise '
            'complete build in reasonable time (PubChem is by far the largest '
            'chemical namespace). Overrides --max-records for PubChem only.'
        ),
    )
    parser.add_argument('-v', '--verbose', action='store_true')

    opts = parser.parse_args(args)

    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if opts.list_presets:
        from omnipath_utils.db._presets import PRESETS

        for name, config in PRESETS.items():
            print(f'{name}: {config["description"]}')
            print(
                f'  Organisms: {config["organisms"]}'
            )
            print(
                f'  Metabolites: {config["metabolite"]},'
                f' miRNA: {config["mirna"]},'
                f' Orthology: {config["orthology"]}'
            )
        return

    from omnipath_utils.db._build import DatabaseBuilder

    builder = DatabaseBuilder(
        db_url=opts.db_url,
        max_records=opts.max_records,
        pubchem_max_records=opts.pubchem_max_records,
    )

    if opts.preset:
        builder.build_preset(opts.preset, parquet_dir=opts.parquet_dir)
    elif opts.ftp:
        builder.build_reference_tables()
        builder.populate_from_ftp()
    elif opts.metabolites:
        builder.build_reference_tables()
        builder.populate_metabolites()
    elif opts.ref_only:
        builder.build_reference_tables()
    elif opts.scope:
        from omnipath_utils.db._presets import resolve_scope

        taxa = resolve_scope(opts.scope)
        builder.build_all(organisms=taxa, scope=opts.scope)
    else:
        builder.build_all(organisms=opts.organisms)

    # Always (re)apply the resolver projection views so an additive/incremental
    # load (e.g. --metabolites, --ftp) never leaves them stale or missing — they
    # are idempotent CREATE OR REPLACE VIEWs over the just-loaded tables, and
    # omnipath-build reads them via DuckDB ATTACH for the full build.
    builder.create_resolver_views()

    print('Build complete!')
