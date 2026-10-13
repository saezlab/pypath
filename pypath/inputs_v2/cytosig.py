"""
Cytokine response signatures of CytoSig.

CytoSig is a model of cytokine signaling activity. Its signature matrix gives,
for each of 4,881 genes and each of 43 cytokines, one value, the response of
the gene to the cytokine; it is the regressor matrix with which CytoSig infers
cytokine activities from expression profiles. The module imports the matrix of
the original signature (``signature.centroid``, 43 cytokines with high
confidence data) as two datasets:

* ``terms``: the 43 cytokines as entities with the names the matrix writes;
* ``signature``: a statement for each cell of the matrix, a gene associated
  with a cytokine, with the value. No cell is left out and none is thresholded,
  the zeros of the matrix included.

The cytokine labels are the matrix's own column names (``IL12``, ``Activin A``,
``NO``, ``IFN1``). Some name a family or a molecule that is not a protein, so
they are source terms and are not turned into protein entities.
"""

from __future__ import annotations

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.measurements import Measurement
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cytosig import iter_cytokines, iter_signature
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Relation
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.CYTOSIG,
    name='CytoSig',
    full='Cytokine Signaling Analyzer',
    url='https://cytosig.ccr.cancer.gov/',
    # LICENSE of https://github.com/data2intelligence/CytoSig: "Users can use CytoSig
    # under the GNU General Public License version 3:
    # https://www.gnu.org/licenses/gpl-3.0.en.html" (checked 2026-10-12). The legacy
    # registry records a "CytoSig" license, non-profit and share-alike.
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='34594031',
    topics=(Topic.SIGNALING,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL, KnowledgeOrigin.COMPUTATIONALLY_PREDICTED),
    summary='Response of 4,881 genes to 43 cytokines, the signature matrix of CytoSig.',
    description=(
        'CytoSig is a database of cytokine response signatures and a model that '
        'predicts cytokine signaling activity from transcriptomic profiles. The '
        'signature matrix, built from cytokine treatment experiments, gives the '
        'response of each gene to each of 43 cytokines.'
    ),
)

# ================================== DOWNLOAD ==================================

# The legacy URL is unchanged. The repository also has two larger signatures
# (``signature.centroid.expand``, 51 cytokines, and ``.beta``, which its README
# calls not fully validated); they are not imported.
download = Download(
    url=(
        'https://raw.githubusercontent.com/data2intelligence/CytoSig/'
        'master/CytoSig/signature.centroid'
    ),
    filename='cytosig_signature.centroid',
    subfolder='cytosig',
    large=True,
)

# =================================== SCHEMA ===================================

CYTOKINE = 'cytosig_cytokine'
# The README of CytoSig takes human gene symbols as the rows of the input
# profiles that the matrix is regressed on.
TAXON = 'NCBITaxon:9606'

# A cell of the matrix. The file gives no unit or definition: the README has the
# expression profiles that the signature is used with as log2(x+1) values,
# preferably differences between two conditions, so the value is a response on
# that scale. It is carried as the source gives it, from -2.46 to 5.27 in the
# file, 1,187 of the 209,883 cells being 0.
RESPONSE = 'cytosig:response'

f = FieldConfig()


def _response(row) -> Measurement | None:
    value = measurement(row.get('value'), source_field='signature.centroid')
    if value is None:
        raise ValueError(
            f'CytoSig cell is not a finite number: {row.get("gene")} / '
            f'{row.get("cytokine")}: {row.get("value")!r}'
        )
    return value


# The matrix says how much a gene responds to a cytokine, not that it is a
# target gene of one or that the cytokine regulates it.
signature_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=Gene,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('gene')),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
    ),
    predicate=slots.associated_with,
    object=EntityBuilder(
        entity_type=OntologyClass,
        identifiers=IdentifiersBuilder(CV(term=CYTOKINE, value=f('cytokine'))),
    ),
    annotations=AnnotationsBuilder(CV(term=RESPONSE, value=_response)),
)


def map_signature(row: dict) -> Relation | None:
    """A gene-cytokine statement; a cell that is not a number raises."""
    _response(row)
    return signature_schema(row)


cytokine_schema = EntityBuilder(
    entity_type=OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=CYTOKINE, value=f('cytokine')),
        CV(term=Namespace.NAME, value=f('cytokine')),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=cytokine_schema,
        raw_parser=iter_cytokines,
    ),
    signature=Dataset(
        download=download,
        mapper=map_signature,
        raw_parser=iter_signature,
    ),
)

# ================================= REFERENCE ==================================

# Activin A  BDNF  BMP2  ...  WNT3A      (43 cytokines; the header has no name for the first column)
# A4GALT  -0.1354...  0.0  -0.0185...    (gene symbol, then one value for each cytokine)
