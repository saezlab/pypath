"""
Parse PTFI Discover data and emit Entity records.

PTFI Discover (Protein-Food Interaction Database) is a comprehensive resource
on food composition and metabolomics. This module creates Food entities with
compound members, including chemical formulas, statistical concentration data,
and compound classifications.

Data sources:
- https://ptfidiscover.markerlab.com/
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional

import requests

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.ptfi import _raw, MEMBER_DELIMITER
from pypath.share.downloads import DATA_DIR


config = ResourceConfig(
    id=ResourceCv.PTFI,
    name='PTFI Discover',
    url='https://ptfidiscover.markerlab.com/',
    license=LicenseCV.UNSPECIFIED,  # TODO: Verify license
    update_category=UpdateCategoryCV.REGULAR,
    description=(
        'PTFI Discover is a comprehensive database of food composition '
        'and metabolomics data, providing detailed information on chemical '
        'compounds found in foods with statistical measurements across '
        'multiple samples.'
    ),
)

# =============================================================================
# Download configuration
# =============================================================================

BASE_API_URL = "https://ptfidiscover.markerlab.com/food/{specimen_id}/detail"
DELAY_SECONDS = 1.0  # Delay between API requests

# List of all PTFI specimen IDs to download (241 total)
SPECIMEN_IDS = [
    "03311673", "03305289", "03310978", "00004515", "03309927", "03301716", "00003157", "03306726", "03309834", "03315874",
    "03301422", "03310193", "00003896", "00004138", "03301722", "03311121", "00004216", "00004415", "00004655", "00004660",
    "02021831", "03000023", "03301002", "03301449", "03301727", "03304854", "03304859", "03310120", "00003209", "00003618",
    "00004137", "00004257", "00004400", "00004428", "00004432", "00004433", "00004597", "00004624", "00004627", "00004653",
    "00004657", "00004676", "00004679", "00004680", "00004682", "00004687", "00004690", "00004699", "00004702", "00004708",
    "00004711", "02021818", "02021844", "03000252", "03000254", "03000261", "03301116", "03301205", "03301369", "03301545",
    "03301697", "03301713", "03301804", "03301828", "03302340", "03302734", "03303012", "03304999", "03305051", "03307110",
    "03307633", "03308022", "03309574", "03309928", "03310137", "03311315", "03311349", "03311757", "03316061", "00003255",
    "00003349", "00003350", "00003417", "00003515", "00003530", "00003532", "00003555", "00003717", "00003777", "00004108",
    "00004131", "00004151", "00004204", "00004206", "00004218", "00004236", "00004253", "00004258", "00004259", "00004365",
    "00004370", "00004379", "00004384", "00004389", "00004395", "00004397", "00004399", "00004403", "00004408", "00004409",
    "00004412", "00004419", "00004421", "00004425", "00004437", "00004440", "00004442", "00004584", "00004590", "00004593",
    "00004598", "00004601", "00004603", "00004606", "00004630", "00004632", "00004637", "00004638", "00004641", "00004645",
    "00004646", "00004650", "00004663", "00004664", "00004667", "00004670", "00004672", "00004675", "00004684", "00004685",
    "00004692", "00004694", "00004696", "00004706", "02000375", "02021161", "02021164", "02021784", "02021857", "02022078",
    "03000055", "03000069", "03000156", "03000244", "03000268", "03000273", "03301040", "03301192", "03301199", "03301223",
    "03301236", "03301312", "03301396", "03301424", "03301441", "03301446", "03301450", "03301556", "03301615", "03301702",
    "03301714", "03301717", "03301719", "03301723", "03301756", "03301777", "03301825", "03301830", "03302045", "03302068",
    "03302374", "03302515", "03302732", "03302805", "03302908", "03302918", "03302922", "03302992", "03303398", "03303502",
    "03303514", "03303578", "03303655", "03303659", "03303675", "03303700", "03303720", "03304495", "03304534", "03304590",
    "03304643", "03304696", "03304860", "03305086", "03305231", "03305454", "03305642", "03305769", "03305803", "03306104",
    "03306125", "03306516", "03306534", "03306597", "03307022", "03307161", "03308314", "03308356", "03308626", "03309091",
    "03309240", "03309584", "03309709", "03309932", "03309980", "03310147", "03310228", "03310357", "03310782", "03310851",
    "03310938", "03311276", "03311331", "03311336", "03311337", "03311407", "03311411", "03311457", "03311801", "03413321",
    "03601045",
]


def download_ptfi_data(
    output_dir: Path | str | None = None,
    specimen_ids: list[str] | None = None,
    force_refresh: bool = False,
) -> Path:
    """
    Download PTFI specimen data from the API.

    Args:
        output_dir: Directory to save JSON files (default: DATA_DIR/ptfi)
        specimen_ids: Optional list of specimen ID numbers (uses built-in list if None)
        force_refresh: If True, re-download existing files

    Returns:
        Path to directory containing downloaded JSON files
    """
    # Set up output directory
    if output_dir is None:
        output_dir = DATA_DIR / 'ptfi'
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Use built-in specimen IDs if none provided
    if specimen_ids is None:
        specimen_ids = SPECIMEN_IDS

    # Convert ID numbers to full specimen IDs
    full_specimen_ids = [f"specimen_FOODON_{id_num}" for id_num in specimen_ids]

    print(f"Found {len(full_specimen_ids)} specimen IDs to download")

    # Download each specimen
    with requests.Session() as session:
        session.headers.update({
            'User-Agent': 'Mozilla/5.0 (compatible; OmniPath/1.0; +https://omnipathdb.org/)'
        })

        successful = 0
        skipped = 0
        failed = 0

        for i, specimen_id in enumerate(full_specimen_ids, 1):
            output_file = output_dir / f"{specimen_id}.json"

            # Skip if already exists and not forcing refresh
            if output_file.exists() and not force_refresh:
                skipped += 1
                if i % 50 == 0:  # Print progress every 50 files
                    print(f"  [{i}/{len(full_specimen_ids)}] Processed (skipped existing files)")
                continue

            # Download from API
            url = BASE_API_URL.format(specimen_id=specimen_id)

            try:
                response = session.get(url, timeout=30)
                response.raise_for_status()
                data = response.json()

                # Save to file
                with open(output_file, 'w') as f:
                    json.dump(data, f)

                successful += 1

                if i % 10 == 0:  # Print progress every 10 downloads
                    print(f"  [{i}/{len(full_specimen_ids)}] Downloaded: {specimen_id}")

                # Be respectful with delays
                if i < len(full_specimen_ids):
                    time.sleep(DELAY_SECONDS)

            except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
                print(f"  ERROR downloading {specimen_id}: {e}")
                failed += 1

    print(f"\nDownload complete: {successful} downloaded, {skipped} skipped, {failed} failed")
    return output_dir


# Custom download handler
download_json = Download(
    url='https://ptfidiscover.markerlab.com/',  # Base URL for reference
    filename='ptfi_data',
    subfolder='ptfi',
    large=True,
)


# =============================================================================
# Field configuration
# =============================================================================

f = FieldConfig(delimiter=MEMBER_DELIMITER, preserve_indices=True)


# =============================================================================
# Food Schema with Compound Membership (declarative)
# =============================================================================

foods_schema = EntityBuilder(
    entity_type=EntityTypeCv.FOOD,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.PTFI, value=f('specimen_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.FOODON, value=f('foodon_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.SAMPLE_COUNT, value=f('num_samples')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.SMALL_MOLECULE,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.PTFI, value=f('member_ptfi_id')),
                CV(term=IdentifierNamespaceCv.NAME, value=f('member_name')),
                CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA, value=f('member_formula')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.COMPOUND_CLASS, value=f('member_category')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_VALUE, value=f('member_value')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MEAN, value=f('member_mean')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MEDIAN, value=f('member_median')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MIN, value=f('member_min')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MAX, value=f('member_max')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_UNIT, value=f('member_units')),
                CV(term=MoleculeAnnotationsCv.SAMPLE_COUNT, value=f('member_n_samples')),
            ),
        ),
    ),
)


# =============================================================================
# Resource definition
# =============================================================================

def _download_and_parse(opener=None, **kwargs):
    """
    Download PTFI data and parse it.

    This is a custom approach since PTFI requires downloading multiple files
    via API rather than a single archive.

    Args:
        opener: Unused (PTFI downloads via API, not a single file)
        **kwargs: Additional arguments passed to _raw parser
    """
    # Download the data (will skip if already exists)
    data_dir = download_ptfi_data()

    # Parse the downloaded data
    yield from _raw(data_dir=str(data_dir), **kwargs)


resource = Resource(
    config,
    foods=Dataset(
        download=download_json,  # Placeholder for metadata
        mapper=foods_schema,
        raw_parser=_download_and_parse,
    ),
)


__all__ = ['config', 'resource', 'download_json', 'download_ptfi_data']
