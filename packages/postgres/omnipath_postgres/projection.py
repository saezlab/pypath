"""Stream unchanged published Parquets into the current main relational model.

This module has no PostgreSQL connection and performs no entity resolution.
The caller owns the DuckDB connection (prefer an on-disk database for a release),
PostgreSQL schema, transactions, partitions and COPY.  Only the small dimension
dictionaries cross the Python boundary; biological rows stay in DuckDB SQL.

Main's canonical relation is an endpoint/predicate triple.  Published qualified
statements retain their evidence links. Molecular references and occurrence
context are required PostgreSQL companions. Other published crosswalks and array
occurrences remain in pinned Parquets, with optional PostgreSQL audit copies.
Entity occurrences absent from the published contract are explicitly represented
as aggregate evidence, with status ``published`` rather than an invented match.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import Any

from omnipath_core.reactions import COMPARTMENT

from omnipath_postgres.parquet_queries import _literal, _scan, quantity_sql, validate_inputs
from omnipath_postgres.releases import PinnedRelease, ResourceSelection


PUBLISHED_STATUS_ID = 5
ENTITY_DATASET = "omnipath:published_entities"
CLAIM_DATASET = "omnipath:published_statement"
PUBLISHED_FALLBACK_NAMESPACE = "omnipath:unresolved_entity_key"
# The entity type of a reaction event; its compartments are the entity
# annotation omnipath_core.reactions.COMPARTMENT.
REACTION_ENTITY_TYPE = "molecular_activity"
ONTOLOGY_CV_NAMESPACE = "Cv Term Accession:OM:0204"
ONTOLOGY_NAME_NAMESPACE = "Name:OM:0202"
ONTOLOGY_SYNONYM_NAMESPACE = "Synonym:OM:0203"
STANDARD_INCHI_NAMESPACE = "Standard Inchi:MI:2010"
PUBLISHED_PROVENANCE_TABLES = (
    "parquet_entity",
    "parquet_statement",
    "parquet_evidence",
    "parquet_identifier_occurrence",
    "parquet_annotation_occurrence",
)

# Known outer CURIE prefixes only. Intrinsic identifier prefixes (CHEMBL,
# HMDB, SLM, C/R in KEGG IDs) remain part of the resulting identifier. Neither
# unknown namespaces nor arbitrary colon-containing strings are shortened.
BARE_IDENTIFIER_RULES = (
    ("chebi", "(?i)^chebi:([0-9]+)$"),
    ("pubchem", "(?i)^(?:pubchem|pubchem.compound|cid):([0-9]+)$"),
    ("chembl", "(?i)^chembl:(CHEMBL[0-9]+)$"),
    ("hmdb", "(?i)^hmdb:(HMDB[0-9]+)$"),
    ("kegg", "(?i)^(?:kegg|kegg.compound|cpd):(C[0-9]{5})$"),
    ("kegg_reaction", "(?i)^(?:kegg.reaction|rn):(R[0-9]{5})$"),
    ("lipidmaps", "(?i)^lipidmaps:(LM[A-Z0-9]+)$"),
    ("swisslipids", "(?i)^swisslipids:(SLM:[0-9]+)$"),
    ("entrez", "(?i)^(?:ncbigene|entrez):([0-9]+)$"),
    ("hgnc", "(?i)^hgnc:([0-9]+)$"),
    ("reactome", "(?i)^reactome:(R-[A-Z]{3}-[0-9]+(?:[.][0-9]+)?)$"),
)

# Explicit choices in inputs_v2. Namespace alone does not imply that a protein,
# pathway or activity typed as an ontology_class is an ontology term.
ONTOLOGY_SCOPES = (
    ("go", "go", "gene_ontology"),
    ("hpo", "hpo", "hpo"),
    ("mondo", "mondo", "mondo"),
    ("bto", "bto", "bto"),
    ("apo", "apo", "apo"),
    ("doid", "doid", "doid"),
    ("pato", "pato", "pato"),
    ("eco", "eco", "eco"),
    ("mp", "mp", "mp"),
    ("wbphenotype", "wbphenotype", "wbphenotype"),
    ("cl", "cl", "cl"),
    ("uberon", "uberon", "uberon"),
    ("oncotree", "oncotree", "oncotree"),
    ("chebi", "chebi", "chebi"),
    ("chemont", "chemont", "chemont"),
    ("psi_mi", "mi", "psi-mi"),
    ("omnipath_ontology", "om", "omnipath"),
    ("reactome", "reactome", "reactome_pathways"),
    ("kegg", "kegg_pathway", "kegg_pathways"),
    ("brenda", "ec", "enzyme_classification"),
    ("uniprot", "uniprot_keyword", "uniprot_keywords"),
    ("swisslipids", "swisslipids", "swisslipids"),
)


@dataclass(frozen=True)
class CopyQuery:
    table: str
    columns: tuple[str, ...]
    query: str


@dataclass(frozen=True)
class AlignedCopyPlan:
    queries: tuple[CopyQuery, ...]
    dimensions: Mapping[str, tuple[tuple[Any, ...], ...]]
    counts: Mapping[str, int]
    sources: tuple[tuple[int, str], ...]
    compatibility: Mapping[str, Any]


# These are additive provenance/quantity tables, not replacement main tables.
# DDL takes a quoted schema token supplied by the loader, not user input.
COMPANION_DDL = {
    "entity_reference_context": """resource text NOT NULL, version text NOT NULL,
        entity_key text NOT NULL, entity_id uuid NOT NULL,
        reference_entity_key text, gene_reference_keys jsonb,
        PRIMARY KEY(resource, version, entity_key)""",
    "statement_reference_context": """resource text NOT NULL, version text NOT NULL,
        relation_key text NOT NULL, subject_reference_entity_key text,
        object_reference_entity_key text,
        PRIMARY KEY(resource, version, relation_key)""",
    "molecular_evidence_context": """resource text NOT NULL, version text NOT NULL,
        owner_kind text NOT NULL CHECK(owner_kind IN ('entity','relation')),
        owner_key text NOT NULL, ordinal bigint NOT NULL CHECK(ordinal >= 0),
        relation_evidence_id uuid, source text, dataset text, row_id text,
        upstream_id text, occurrence_json jsonb NOT NULL,
        PRIMARY KEY(resource, version, owner_kind, owner_key, ordinal)""",
    "parquet_entity": """resource text NOT NULL, version text NOT NULL,
        entity_key text NOT NULL, entity_id uuid NOT NULL,
        entity_evidence_id uuid NOT NULL, label text, namespace text,
        identifier text, entity_type text, taxon text,
        has_hierarchy boolean, parent_count bigint, child_count bigint,
        identifiers_present boolean NOT NULL, annotations_present boolean NOT NULL,
        PRIMARY KEY(resource, version, entity_key)""",
    "parquet_statement": """resource text NOT NULL, version text NOT NULL,
        relation_key text NOT NULL, relation_id uuid, statement_kind text,
        subject_entity_id uuid NOT NULL, object_entity_id uuid NOT NULL,
        predicate_id bigint NOT NULL, subject_label text, subject_type text,
        object_label text, object_type text, taxon text, is_directed boolean,
        sign integer, category text, interaction_class text, evidence_count bigint,
        sources text[], evidence_present boolean NOT NULL,
        annotations_present boolean NOT NULL,
        PRIMARY KEY(resource, version, relation_key)""",
    "parquet_evidence": """resource text NOT NULL, version text NOT NULL,
        relation_key text NOT NULL, ordinal bigint NOT NULL,
        source_id bigint NOT NULL, relation_evidence_id uuid NOT NULL,
        source text, dataset text, original_row_id text, upstream_id text,
        annotations_present boolean NOT NULL, synthetic boolean NOT NULL,
        PRIMARY KEY(resource, version, relation_key, ordinal)""",
    "parquet_identifier_occurrence": """resource text NOT NULL, version text NOT NULL,
        entity_key text NOT NULL, ordinal bigint NOT NULL, identifier_id uuid,
        ns text, identifier text, is_canonical boolean, source text,
        PRIMARY KEY(resource, version, entity_key, ordinal)""",
    "parquet_annotation_occurrence": """resource text NOT NULL, version text NOT NULL,
        owner_kind text NOT NULL CHECK(owner_kind IN ('entity','relation','evidence')),
        owner_key text NOT NULL, evidence_ordinal bigint NOT NULL,
        ordinal bigint NOT NULL, annotation_key uuid NOT NULL,
        term text, untyped_value text, source text, dataset text, scope text,
        PRIMARY KEY(resource, version, owner_kind, owner_key, evidence_ordinal, ordinal)""",
    "annotation_quantity": """annotation_key uuid PRIMARY KEY,
        has_numeric_value double precision, has_unit text, has_unit_prefix text,
        has_binary_relation text, source_field text, comparator text, published_value text""",
}


def companion_ddl(
    schema: str,
    *,
    retain_published_provenance: bool = False,
) -> tuple[str, ...]:
    """Create required molecular/quantity context and requested audit copies."""
    if type(retain_published_provenance) is not bool:
        raise ValueError("retain_published_provenance must be a boolean")
    if not isinstance(schema, str) or not schema or "\x00" in schema:
        raise ValueError("schema must be nonempty text without NUL")
    quoted = '"' + schema.replace('"', '""') + '"'
    return tuple(
        f'CREATE TABLE {quoted}."{name}" ({definition})'
        for name, definition in COMPANION_DDL.items()
        if retain_published_provenance or name not in PUBLISHED_PROVENANCE_TABLES
    )


def _resources(
    release: PinnedRelease | Iterable[ResourceSelection],
) -> tuple[ResourceSelection, ...]:
    resources = release.resources if isinstance(release, PinnedRelease) else tuple(release)
    if not resources:
        raise ValueError("An aligned projection requires at least one pinned resource")
    if len({resource.source for resource in resources}) != len(resources):
        raise ValueError("A release must select exactly one version per resource")
    return tuple(sorted(resources, key=lambda resource: resource.source))


def _uuid(payload: str) -> str:
    return f"ap_uuid({payload})"


def _key(*expressions: str) -> str:
    # JSON arrays are unambiguous even when data includes separators or quotes.
    return "to_json(list_value(" + ",".join(f"{item}::VARCHAR" for item in expressions) + "))"


_ENTITY_COLUMNS = (
    "entity_key, entity_type, namespace, identifier, taxon, label, has_hierarchy, "
    "parent_count, child_count, reference_entity_key, gene_reference_keys, "
    "identifier_count, annotation_count"
)
_RELATION_COLUMNS = (
    "relation_key, statement_kind, subject_entity_key, subject_reference_entity_key, "
    "subject_label, subject_type, predicate, object_entity_key, object_reference_entity_key, "
    "object_label, object_type, taxon, is_directed, sign, category, interaction_class, "
    "sources, evidence_count, annotation_count"
)


def _create_inputs(con: Any, resources: tuple[ResourceSelection, ...]) -> None:
    """Stage each pinned resource's published tables, with resource and version.

    Child rows carry their owner's key (entity_key / relation_key) instead of the
    resource-local entity_id / relation_id, and an ``item`` struct of their fields.
    """
    con.execute("""CREATE OR REPLACE MACRO ap_uuid(value) AS (
        (substr(md5(value),1,8)||'-'||substr(md5(value),9,4)||'-'||
         substr(md5(value),13,4)||'-'||substr(md5(value),17,4)||'-'||
         substr(md5(value),21,12))::UUID)""")
    quantity = quantity_sql()

    def stage(table: str, select: Callable[[Callable[[str], str]], str]) -> None:
        parts = []
        for resource in resources:

            def scan(name: str, resource: ResourceSelection = resource) -> str:
                return _scan(resource.directory, f"{name}.parquet")

            parts.append(
                f"SELECT {_literal(resource.source)}::VARCHAR AS resource, "
                f'{_literal(resource.version)}::VARCHAR AS "version", q.* '
                f"FROM ({select(scan)}) q"
            )
        con.execute(f"CREATE OR REPLACE TABLE {table} AS " + " UNION ALL ".join(parts))

    stage("ap_entity_raw", lambda scan: f"SELECT {_ENTITY_COLUMNS} FROM {scan('entity')}")
    stage("ap_statement_raw", lambda scan: f"SELECT {_RELATION_COLUMNS} FROM {scan('relation')}")
    stage(
        "ap_identifier_raw",
        lambda scan: f"""SELECT e.entity_key,
            struct_pack(ns:=i.ns,id:=i.id,is_canonical:=i.is_canonical,source:=i.source) item,
            i.ordinal::BIGINT ordinal
            FROM {scan("entity_identifier")} i JOIN {scan("entity")} e USING(entity_id)""",
    )
    stage(
        "ap_evidence_raw",
        lambda scan: f"""SELECT r.relation_key,
            struct_pack(source:=v.source,dataset:=v.dataset,row_id:=v.row_id,
                upstream_id:=v.upstream_id,annotations:=v.annotations,
                subject_molecular_form:=v.subject_molecular_form,
                object_molecular_form:=v.object_molecular_form) item,
            v.ordinal::BIGINT ordinal
            FROM {scan("relation_evidence")} v JOIN {scan("relation")} r USING(relation_id)""",
    )
    stage(
        "ap_entity_evidence_raw",
        lambda scan: f"""SELECT e.entity_key,
            struct_pack(source:=v.source,dataset:=v.dataset,row_id:=v.row_id,
                upstream_id:=v.upstream_id,annotations:=v.annotations,
                molecular_form:=v.molecular_form) item,
            v.ordinal::BIGINT ordinal
            FROM {scan("entity_evidence")} v JOIN {scan("entity")} e USING(entity_id)""",
    )
    stage(
        "ap_annotation_raw",
        lambda scan: f"""SELECT 'entity'::VARCHAR owner_kind,e.entity_key owner_key,
            -1::BIGINT evidence_ordinal,a.ordinal::BIGINT ordinal,a.term,a.value,
            {quantity} quantity,a.source,a.dataset,NULL::VARCHAR AS "scope"
            FROM {scan("entity_annotation")} a JOIN {scan("entity")} e USING(entity_id)
            UNION ALL
            SELECT 'relation',r.relation_key,-1::BIGINT,a.ordinal::BIGINT,a.term,a.value,
            {quantity},a.source,a.dataset,a.scope
            FROM {scan("relation_annotation")} a JOIN {scan("relation")} r USING(relation_id)""",
    )
    con.execute("""INSERT INTO ap_annotation_raw
        SELECT resource,version,'evidence',relation_key,evidence_ordinal,ordinal,
            a.term,a.value,a.quantity,a.source,a.dataset,a.scope
        FROM (SELECT resource,version,relation_key,ordinal evidence_ordinal,
                     (generate_subscripts(item.annotations,1)-1)::BIGINT ordinal,
                     unnest(item.annotations) a
              FROM ap_evidence_raw)""")
    validate_inputs(con)


def _namespace_names() -> dict[str, str]:
    """Namespace spelling adapters only; no identifier translation/resolution."""
    from pypath.internals.cv_terms import IdentifierNamespaceCv, cv_term_label_accession

    members = {
        "uniprot": "UNIPROT",
        "uniprot_entry": "UNIPROT_ENTRY_NAME",
        "genesymbol": "GENE_NAME_PRIMARY",
        "genesymbol-syn": "GENE_NAME_SYNONYM",
        "entrez": "ENTREZ",
        "hgnc": "HGNC",
        "ensembl": "ENSEMBL",
        "ensg": "ENSEMBL",
        "ensp": "ENSEMBL",
        "enst": "ENSEMBL",
        "chebi": "CHEBI",
        "hmdb": "HMDB",
        "lipidmaps": "LIPIDMAPS",
        "swisslipids": "SWISSLIPIDS",
        "pubchem": "PUBCHEM_COMPOUND",
        "inchikey": "STANDARD_INCHI_KEY",
        "smiles": "SMILES",
        "chembl": "CHEMBL_COMPOUND",
        "kegg": "KEGG_COMPOUND",
        "kegg_reaction": "KEGG_REACTION",
        "cas": "CAS",
        "name": "NAME",
        "synonym": "SYNONYM",
        "cv_term": "CV_TERM_ACCESSION",
        "ramp": "RAMP_ID",
        "refmet": "REFMET",
        "mirbase_precursor": "MIRBASE_PRECURSOR",
        "mirbase_mature": "MIRBASE_MATURE",
        "drugbank": "DRUGBANK",
        "reactome": "REACTOME_STABLE_ID",
        "lipid_name": "LIPID_NAME",
    }
    return {
        name: cv_term_label_accession(getattr(IdentifierNamespaceCv, member))
        for name, member in members.items()
        if hasattr(IdentifierNamespaceCv, member)
    }


def _generic_identifier_included(namespace: str, value: str) -> str:
    """Main omits Standard InChI from generic identifiers, regardless of size.

    Raw published occurrences remain untouched. This is an exact namespace
    adapter for ingest.common.include_identifier, not a general length filter.
    """
    return f"""{namespace} IS NOT NULL AND {value} IS NOT NULL
        AND {namespace} NOT IN ('inchi',{_literal(STANDARD_INCHI_NAMESPACE)})"""


def _dimension(con: Any, table: str, required: Iterable[str], existing: Iterable[tuple]) -> tuple:
    rows = tuple(existing)
    mapping = {str(row[1]): int(row[0]) for row in rows}
    if len(mapping) != len(rows) or len({row[0] for row in rows}) != len(rows):
        raise ValueError(f"Duplicate seeded dimension name or ID in {table}")
    next_id = max(mapping.values(), default=0)
    for name in sorted(set(required) - mapping.keys()):
        next_id += 1
        mapping[name] = next_id
    con.execute(f"CREATE OR REPLACE TABLE ap_{table} (id BIGINT,name VARCHAR)")
    result = tuple(sorted((identifier, name) for name, identifier in mapping.items()))
    if result:
        con.executemany(f"INSERT INTO ap_{table} VALUES (?,?)", result)
    return result


def _ontology_cv_eligible(alias: str) -> str:
    # The former generic ontology mapper used CV_TERM_ACCESSION. Keep only
    # explicit provider/type namespaces whose current published spelling lost
    # that lookup gate; identifiers and canonical identities stay unchanged.
    return f"""(
        ({alias}.entity_type='ontology_class' AND lower({alias}.resource) IN
            ('go','hpo','mondo','chemont','psi_mi','omnipath_ontology','brenda','bto','cl','uberon','apo',
             'doid','pato','eco','mp','wbphenotype','oncotree'))
        OR ({alias}.entity_type='pathway' AND lower({alias}.resource)='kegg'
            AND {alias}.namespace IN ('kegg_pathway','kegg_pathway_category'))
        OR ({alias}.entity_type='ontology_class' AND lower({alias}.resource)='uniprot'
            AND {alias}.namespace='uniprot_keyword')
    ) AND {alias}.namespace NOT IN
        ('name','synonym','Name:OM:0202','Synonym:OM:0203',
         'omnipath:unresolved_entity_key')"""


def _ontology_name_eligible(alias: str) -> str:
    return f"""{alias}.entity_type IN
        ('ontology_class','pathway','chemical_entity','small_molecule')
        AND lower({alias}.resource) IN
            (SELECT resource FROM ap_ontology_source_scope)"""


def _prepare_dimensions(con: Any, existing: Mapping[str, Iterable[tuple]]) -> dict[str, tuple]:
    aliases = _namespace_names()
    con.execute("CREATE OR REPLACE TABLE ap_namespace_alias (raw VARCHAR,name VARCHAR)")
    if aliases:
        con.executemany("INSERT INTO ap_namespace_alias VALUES (?,?)", sorted(aliases.items()))
    con.execute(
        "CREATE OR REPLACE TABLE ap_ontology_scope (resource VARCHAR,namespace VARCHAR,ontology_id VARCHAR)"
    )
    con.executemany("INSERT INTO ap_ontology_scope VALUES (?,?,?)", ONTOLOGY_SCOPES)
    # The audited providers attach one explicit ontology root to each ontology
    # record. Resolution can replace a source namespace (e.g. ChEBI) with an
    # InChIKey, so that canonical namespace cannot select the record's family.
    # Ambiguous future declarations get no source-only fallback.
    con.execute("""CREATE OR REPLACE TABLE ap_ontology_source_scope AS
        SELECT resource,min(ontology_id) ontology_id FROM ap_ontology_scope
        GROUP BY resource HAVING min(ontology_id)=max(ontology_id)""")
    required_sql = {
        "data_source": """SELECT resource FROM ap_entity_raw UNION SELECT resource FROM
            ap_statement_raw UNION SELECT item.source FROM ap_evidence_raw
            WHERE item.source IS NOT NULL""",
        "vocab_identifier_type": f"""SELECT coalesce(a.name,e.namespace) FROM ap_entity_raw e
            LEFT JOIN ap_namespace_alias a ON a.raw=e.namespace
            UNION SELECT coalesce(a.name,i.item.ns) FROM ap_identifier_raw i
            LEFT JOIN ap_namespace_alias a ON a.raw=i.item.ns WHERE i.item.ns IS NOT NULL
            UNION SELECT {_literal(PUBLISHED_FALLBACK_NAMESPACE)}
            UNION SELECT {_literal(ONTOLOGY_CV_NAMESPACE)} FROM ap_entity_raw e
                  WHERE {_ontology_cv_eligible("e")}
            UNION SELECT {_literal(ONTOLOGY_NAME_NAMESPACE)} FROM ap_entity_raw e
                  WHERE {_ontology_name_eligible("e")}""",
        "vocab_entity_type": "SELECT entity_type FROM ap_entity_raw WHERE entity_type IS NOT NULL",
        "vocab_relation_predicate": "SELECT predicate FROM ap_statement_raw WHERE predicate IS NOT NULL",
        "vocab_relation_category": "SELECT coalesce(category,'omnipath:unspecified') FROM ap_statement_raw",
        "vocab_annotation_scope": "SELECT scope FROM ap_annotation_raw WHERE scope IS NOT NULL AND scope<>'evidence'",
    }
    result = {}
    seeds = {
        "vocab_entity_role": ((1, "parent"), (2, "member")),
        "vocab_annotation_scope": ((1, "relation"), (2, "subject"), (3, "object")),
        "vocab_resolution_status": (
            (1, "resolved"),
            (2, "unresolved"),
            (3, "ambiguous"),
            (4, "unsupported"),
            (PUBLISHED_STATUS_ID, "published"),
        ),
    }
    for table, query in required_sql.items():
        # Only vocabulary names cross into Python. Some source SELECTs contain
        # millions of biological rows; deduplicating in _dimension after
        # fetchall would first allocate those repeated names in Python memory.
        names_query = (
            f"SELECT DISTINCT dimension_name FROM ({query}) "
            "required(dimension_name) WHERE dimension_name IS NOT NULL"
        )
        required = [row[0] for row in con.execute(names_query).fetchall()]
        result[table] = _dimension(con, table, required, existing.get(table, seeds.get(table, ())))
    for table in ("vocab_entity_role", "vocab_resolution_status"):
        result[table] = _dimension(
            con, table, (row[1] for row in seeds[table]), existing.get(table, seeds[table])
        )
    status = dict((name, identifier) for identifier, name in result["vocab_resolution_status"])
    if status["published"] != PUBLISHED_STATUS_ID:
        raise ValueError("Published status must have the main-compatible explicit ID 5")
    datasets = tuple(existing.get("dataset", ()))
    dataset_map = {(int(row[1]), str(row[2])): int(row[0]) for row in datasets}
    source_ids = {name: identifier for identifier, name in result["data_source"]}
    required = con.execute(f"""SELECT resource,{_literal(ENTITY_DATASET)} FROM ap_entity_raw
        UNION SELECT coalesce(item.source,resource),coalesce(item.dataset,{_literal(CLAIM_DATASET)})
              FROM ap_evidence_raw
        UNION SELECT resource,{_literal(CLAIM_DATASET)} FROM ap_statement_raw""").fetchall()
    next_id = max(dataset_map.values(), default=0)
    for source, name in sorted(required):
        key = (source_ids[source], name)
        if key not in dataset_map:
            next_id += 1
            dataset_map[key] = next_id
    con.execute("CREATE OR REPLACE TABLE ap_dataset (id BIGINT,source_id BIGINT,name VARCHAR)")
    result["dataset"] = tuple(sorted((identifier, *key) for key, identifier in dataset_map.items()))
    if result["dataset"]:
        con.executemany("INSERT INTO ap_dataset VALUES (?,?,?)", result["dataset"])
    return result


def _checked(con: Any, query: str, message: str) -> None:
    count = con.execute(f"SELECT count(*) FROM ({query}) invalid").fetchone()[0]
    if count:
        # Deliberately aggregate only: diagnostics must not dump biological rows.
        raise ValueError(f"{message}: {count} conflicting rows/groups")


def _validate_canonical_inputs(con: Any) -> None:
    _checked(
        con,
        """SELECT entity_key FROM ap_entity_raw WHERE entity_type IS NULL
        OR namespace IS NULL OR identifier IS NULL""",
        "Main canonical entities require published type, namespace and identifier",
    )
    _checked(
        con,
        """SELECT entity_key FROM ap_entity_raw GROUP BY entity_key
        HAVING count(DISTINCT to_json(struct_pack(entity_type:=entity_type,
               namespace:=namespace,identifier:=identifier)))>1""",
        "A published entity identity has conflicting canonical fields across resources",
    )
    _checked(
        con,
        """SELECT resource,version,entity_key FROM ap_entity_raw
        GROUP BY resource,version,entity_key HAVING count(*)>1""",
        "Duplicate published entity identity within one resource",
    )
    _checked(
        con,
        """SELECT resource,version,relation_key FROM ap_statement_raw
        GROUP BY resource,version,relation_key HAVING count(*)>1""",
        "Duplicate published statement identity within one resource",
    )
    # The publisher hashes endpoints, predicate, kind and Biolink qualifiers.
    # Its scalar taxon is a per-resource consensus, not part of that key (see
    # writer._relations_sql). Keep that consensus and the exact in_taxon claims
    # scoped to resource/version; globally equating them rejects valid sources.
    # Retain the existing derived-field checks until a demonstrated publisher
    # version difference warrants a separately reviewed presentation policy.
    _checked(
        con,
        """SELECT relation_key FROM ap_statement_raw GROUP BY relation_key
        HAVING count(DISTINCT to_json(struct_pack(subject:=subject_entity_key,
            object:=object_entity_key,predicate:=predicate,kind:=statement_kind,
            directed:=is_directed,sign:=sign,category:=category,
            interaction_class:=interaction_class)))>1""",
        "A published statement identity has conflicting scientific fields",
    )
    _checked(
        con,
        """SELECT e.entity_key FROM ap_entity_raw e
        WHERE nullif(e.taxon,'') IS NOT NULL AND try_cast(e.taxon AS BIGINT) IS NULL""",
        "Published taxonomy cannot be represented in main bigint taxonomy_id",
    )
    _checked(
        con,
        """SELECT r.relation_key FROM ap_statement_raw r
        LEFT JOIN ap_entity_raw s ON s.resource=r.resource AND s.version=r.version
             AND s.entity_key=r.subject_entity_key
        LEFT JOIN ap_entity_raw o ON o.resource=r.resource AND o.version=r.version
             AND o.entity_key=r.object_entity_key
        WHERE s.entity_key IS NULL OR o.entity_key IS NULL OR r.predicate IS NULL""",
        "Published statement has missing canonical endpoint or predicate",
    )


def _create_entity_projection(con: Any) -> None:
    # Published identity is authoritative. A null/known taxon disagreement does
    # not mint a second entity; two distinct known taxa make canonical taxonomy
    # unknown. Every exact published occurrence remains in evidence/crosswalk.
    con.execute("""CREATE OR REPLACE TABLE ap_entity_canonical_input AS SELECT
        entity_key,entity_type,namespace,identifier,
        CASE WHEN min(try_cast(taxon AS BIGINT))=max(try_cast(taxon AS BIGINT))
             THEN min(try_cast(taxon AS BIGINT)) ELSE NULL END taxonomy_id
        FROM ap_entity_raw GROUP BY entity_key,entity_type,namespace,identifier""")
    # A reaction is published once per set of compartments it takes place in:
    # the same MetaNetX or Rhea reaction in the cytosol and in the lysosome are
    # two events, two entities that share the accession. They are told apart by
    # the compartments each one states.
    con.execute(f"""CREATE OR REPLACE TABLE ap_entity_compartments AS
        SELECT owner_key entity_key,
            array_to_string(list(DISTINCT "value" ORDER BY "value"),',') compartments
        FROM ap_annotation_raw
        WHERE owner_kind='entity' AND term={_literal(COMPARTMENT)}
        GROUP BY owner_key""")
    con.execute(f"""CREATE OR REPLACE TABLE ap_entity_natural_conflict AS
        SELECT e.entity_type,e.taxonomy_id,coalesce(a.name,e.namespace) namespace_name,e.identifier,
            bool_and(e.entity_type={_literal(REACTION_ENTITY_TYPE)})
            AND count(DISTINCT coalesce(p.compartments,''))=count(*) compartment_variants
        FROM ap_entity_canonical_input e LEFT JOIN ap_namespace_alias a ON a.raw=e.namespace
        LEFT JOIN ap_entity_compartments p ON p.entity_key=e.entity_key
        GROUP BY e.entity_type,e.taxonomy_id,coalesce(a.name,e.namespace),e.identifier
        HAVING count(*)>1""")
    # Entities sharing a natural key take their published key as canonical
    # identifier (below) when the key is a name, or when they are compartment
    # variants of one reaction, each with compartments of its own. Any other
    # shared accession is two identities the publisher should have merged.
    _checked(
        con,
        """SELECT e.entity_key FROM ap_entity_canonical_input e
        JOIN ap_entity_natural_conflict c ON c.entity_type=e.entity_type
            AND c.taxonomy_id IS NOT DISTINCT FROM e.taxonomy_id AND c.identifier=e.identifier
        LEFT JOIN ap_namespace_alias a ON a.raw=e.namespace
        WHERE c.namespace_name=coalesce(a.name,e.namespace)
          AND e.namespace NOT IN ('name','synonym')
          AND NOT c.compartment_variants""",
        "Distinct published entities violate main canonical natural-key uniqueness",
    )
    con.execute(f"""CREATE OR REPLACE TABLE ap_entity AS SELECT
        {_uuid(_key(_literal("published-entity"), "e.entity_key"))} entity_id,
        et.id entity_type_id,e.taxonomy_id,
        CASE WHEN c.identifier IS NOT NULL THEN fallback.id ELSE it.id END canonical_identifier_type_id,
        CASE WHEN c.identifier IS NOT NULL THEN e.entity_key ELSE e.identifier END canonical_identifier,
        {PUBLISHED_STATUS_ID}::SMALLINT resolution_status_id,
        'published_parquet'::VARCHAR resolution_mechanism
        FROM ap_entity_canonical_input e
        JOIN ap_vocab_entity_type et ON et.name=e.entity_type
        LEFT JOIN ap_namespace_alias a ON a.raw=e.namespace
        JOIN ap_vocab_identifier_type it ON it.name=coalesce(a.name,e.namespace)
        JOIN ap_vocab_identifier_type fallback ON fallback.name={_literal(PUBLISHED_FALLBACK_NAMESPACE)}
        LEFT JOIN ap_entity_natural_conflict c ON c.entity_type=e.entity_type
            AND c.taxonomy_id IS NOT DISTINCT FROM e.taxonomy_id AND c.identifier=e.identifier
            AND c.namespace_name=coalesce(a.name,e.namespace)""")
    _checked(
        con,
        """SELECT entity_type_id,taxonomy_id,canonical_identifier_type_id,
        canonical_identifier FROM ap_entity GROUP BY ALL HAVING count(*)>1""",
        "Distinct published entities violate main canonical natural-key uniqueness",
    )
    con.execute(f"""CREATE OR REPLACE VIEW ap_entity_occurrence AS SELECT e.*,
        {_uuid(_key(_literal("published-entity"), "e.entity_key"))} entity_id,
        {_uuid(_key(_literal("aggregate-entity-evidence"), "e.resource", "e.version", "e.entity_key"))}
            entity_evidence_id,
        ds.id source_id,d.id dataset_id,et.id entity_type_id,
        hash({_key("e.resource", "e.entity_key")})::UBIGINT % 9223372036854775807 row_id
        FROM ap_entity_raw e JOIN ap_data_source ds ON ds.name=e.resource
        JOIN ap_dataset d ON d.source_id=ds.id AND d.name={_literal(ENTITY_DATASET)}
        JOIN ap_vocab_entity_type et ON et.name=e.entity_type""")
    # Dictionary deduplication mirrors main. Original null structs/occurrences
    # remain in the narrow companion rather than acquiring invented namespaces.
    con.execute("""CREATE OR REPLACE TABLE ap_identifier_input AS
        SELECT namespace ns,identifier FROM ap_entity_raw
        UNION SELECT item.ns,item.id FROM ap_identifier_raw
        UNION SELECT it.name,e.canonical_identifier FROM ap_entity e
              JOIN ap_vocab_identifier_type it ON it.id=e.canonical_identifier_type_id""")
    con.execute(
        "CREATE OR REPLACE TABLE ap_bare_identifier_rule (namespace VARCHAR,pattern VARCHAR)"
    )
    con.executemany("INSERT INTO ap_bare_identifier_rule VALUES (?,?)", BARE_IDENTIFIER_RULES)
    con.execute(f"""CREATE OR REPLACE TABLE ap_identifier_alias AS SELECT DISTINCT
        {_uuid(_key(_literal("identifier"), "coalesce(a.name,i.ns)", "i.identifier"))} original_identifier_id,
        {_uuid(_key(_literal("identifier"), "coalesce(a.name,i.ns)", "regexp_extract(i.identifier,rule.pattern,1)"))} alias_identifier_id,
        coalesce(a.name,i.ns) namespace_name,
        regexp_extract(i.identifier,rule.pattern,1) alias_value
        FROM ap_identifier_input i JOIN ap_bare_identifier_rule rule ON rule.namespace=i.ns
        LEFT JOIN ap_namespace_alias a ON a.raw=i.ns
        WHERE regexp_full_match(i.identifier,rule.pattern)""")
    con.execute(f"""CREATE OR REPLACE TABLE ap_identifier AS SELECT DISTINCT
        {_uuid(_key(_literal("identifier"), "coalesce(a.name,i.ns)", "i.identifier"))} identifier_id,
        it.id identifier_type_id,i.identifier AS "value",
        NULL::VARCHAR value_normalized
        FROM (SELECT ns,identifier FROM ap_identifier_input
              UNION SELECT namespace_name,alias_value FROM ap_identifier_alias) i
        LEFT JOIN ap_namespace_alias a ON a.raw=i.ns
        JOIN ap_vocab_identifier_type it ON it.name=coalesce(a.name,i.ns)
        WHERE {_generic_identifier_included("coalesce(a.name,i.ns)", "i.identifier")}""")
    con.execute(f"""CREATE OR REPLACE VIEW ap_identifier_occurrence AS SELECT i.*,
        CASE WHEN {_generic_identifier_included("coalesce(a.name,i.item.ns)", "i.item.id")} THEN
            {_uuid(_key(_literal("identifier"), "coalesce(a.name,i.item.ns)", "i.item.id"))}
        END identifier_id
        FROM ap_identifier_raw i LEFT JOIN ap_namespace_alias a ON a.raw=i.item.ns""")


def _create_statement_projection(con: Any) -> None:
    con.execute(f"""CREATE OR REPLACE VIEW ap_statement AS SELECT r.*,
        s.entity_id subject_entity_id,o.entity_id object_entity_id,rp.id predicate_id,
        rc.id relation_category_id,
        coalesce(source_scope.ontology_id,os.ontology_id,
                 'omnipath:published-scope:'||r.resource||':'||s.namespace)
            ontology_id,
        (source_scope.ontology_id IS NOT NULL OR os.ontology_id IS NOT NULL)
            ontology_scope_recovered,
        CASE WHEN r.statement_kind='relation' THEN
            {_uuid(_key(_literal("canonical-relation"), "s.entity_id", "rp.name", "o.entity_id"))}
        END relation_id
        FROM ap_statement_raw r
        JOIN ap_entity_occurrence s ON s.resource=r.resource AND s.version=r.version
             AND s.entity_key=r.subject_entity_key
        JOIN ap_entity_occurrence o ON o.resource=r.resource AND o.version=r.version
             AND o.entity_key=r.object_entity_key
        JOIN ap_vocab_relation_predicate rp ON rp.name=r.predicate
        JOIN ap_vocab_relation_category rc ON rc.name=coalesce(r.category,'omnipath:unspecified')
        LEFT JOIN ap_ontology_scope os ON os.resource=lower(r.resource) AND os.namespace=s.namespace
        AND (s.namespace<>'reactome' OR s.entity_type='pathway')
        LEFT JOIN ap_ontology_source_scope source_scope
          ON source_scope.resource=lower(r.resource)
         AND r.statement_kind IN ('ontology','ontology_axiom')""")
    con.execute("""CREATE OR REPLACE TABLE ap_relation AS SELECT relation_id,
        subject_entity_id,predicate_id,object_entity_id,
        CASE WHEN count(DISTINCT relation_category_id)=1 THEN min(relation_category_id)
             ELSE NULL END relation_category_id
        FROM ap_statement WHERE relation_id IS NOT NULL
        GROUP BY relation_id,subject_entity_id,predicate_id,object_entity_id""")
    con.execute(f"""CREATE OR REPLACE TABLE ap_evidence AS
        SELECT e.resource,e.version,e.relation_key,e.ordinal,
            {
        _uuid(
            _key(
                _literal("published-relation-evidence"),
                "e.resource",
                "e.version",
                "e.relation_key",
                "e.ordinal",
            )
        )
    } relation_evidence_id,
            ds.id source_id,d.id dataset_id,
            coalesce(try_cast(e.item.row_id AS BIGINT),
                     (hash({_key("e.resource", "e.item.source", "e.item.dataset", "e.item.row_id")})
                      % 9223372036854775807)::BIGINT) row_id,
            e.item.source AS "source",e.item.dataset AS dataset,e.item.row_id original_row_id,
            e.item.upstream_id upstream_id,e.item.annotations IS NOT NULL annotations_present,
            false synthetic
        FROM ap_evidence_raw e
        JOIN ap_data_source ds ON ds.name=coalesce(e.item.source,e.resource)
        JOIN ap_dataset d ON d.source_id=ds.id
             AND d.name=coalesce(e.item.dataset,{_literal(CLAIM_DATASET)})
        UNION ALL
        SELECT r.resource,r.version,r.relation_key,-1::BIGINT,
            {
        _uuid(
            _key(
                _literal("published-statement-evidence"),
                "r.resource",
                "r.version",
                "r.relation_key",
            )
        )
    }, ds.id,d.id,
            (hash({_key("r.resource", "r.relation_key")}) % 9223372036854775807)::BIGINT,
            NULL::VARCHAR,NULL::VARCHAR,NULL::VARCHAR,NULL::VARCHAR,false,true
        FROM ap_statement r
        JOIN ap_data_source ds ON ds.name=r.resource
        JOIN ap_dataset d ON d.source_id=ds.id AND d.name={_literal(CLAIM_DATASET)}
        ANTI JOIN ap_evidence_raw e USING(resource,version,relation_key)""")
    _checked(
        con,
        """SELECT source_id,dataset_id,row_id FROM ap_evidence
        WHERE NOT synthetic AND original_row_id IS NOT NULL
        GROUP BY source_id,dataset_id,row_id
        HAVING count(DISTINCT coalesce(try_cast(original_row_id AS BIGINT)::VARCHAR,
                                      original_row_id))>1""",
        "Original source row identities collide in the main bigint row representation",
    )


def _create_ontology_lookup_aliases(con: Any) -> None:
    """Add narrow ontology reader aliases without inventing source occurrences.

    Ownership comes from actual known-source ontology edges. Only primitive
    positive owners are staged; global name/value exclusions first restrict the
    published identifier scan to their entity keys. Neither raw arrays nor
    evidence/annotation payloads enter this stage.
    """
    con.execute("""CREATE OR REPLACE TABLE ap_ontology_lookup_owner AS
        SELECT e.resource,e.version,e.entity_key,e.entity_type,e.namespace,
               e.identifier,e.label
        FROM ap_entity_raw e SEMI JOIN (
            SELECT resource,version,subject_entity_key AS entity_key
            FROM ap_statement_raw WHERE statement_kind IN ('ontology','ontology_axiom')
            UNION
            SELECT resource,version,object_entity_key AS entity_key
            FROM ap_statement_raw WHERE statement_kind IN ('ontology','ontology_axiom')
        ) edges USING(resource,version,entity_key)
        WHERE lower(e.resource) IN (SELECT resource FROM ap_ontology_source_scope)""")
    cv_namespace = _literal(ONTOLOGY_CV_NAMESPACE)
    name_namespace = _literal(ONTOLOGY_NAME_NAMESPACE)
    con.execute(f"""CREATE OR REPLACE TABLE ap_ontology_lookup_alias AS
        SELECT resource,version,entity_key,{cv_namespace} AS namespace_name,
               identifier AS alias_value,
               {_uuid(_key(_literal("identifier"), cv_namespace, "identifier"))} alias_identifier_id
        FROM ap_ontology_lookup_owner e
        WHERE {_ontology_cv_eligible("e")} AND identifier<>''""")
    # A scalar published label is useful only if it is a genuine missing name.
    # Respect every source's original names and synonyms for the same identity.
    con.execute(f"""CREATE OR REPLACE TABLE ap_ontology_name_candidate AS
        SELECT e.resource,e.version,e.entity_key,e.label
        FROM ap_ontology_lookup_owner e
        WHERE {_ontology_name_eligible("e")}
          AND e.label IS NOT NULL AND e.label<>''
          AND e.label<>e.identifier AND e.label<>e.entity_key""")
    con.execute(f"""CREATE OR REPLACE TABLE ap_ontology_original_named AS
        SELECT DISTINCT i.entity_key
        FROM ap_identifier_raw i SEMI JOIN ap_ontology_name_candidate c
             ON c.entity_key=i.entity_key
        LEFT JOIN ap_namespace_alias ns ON ns.raw=i.item.ns
        WHERE i.item.id IS NOT NULL AND i.item.id<>''
          AND coalesce(ns.name,i.item.ns) IN
              ({name_namespace},{_literal(ONTOLOGY_SYNONYM_NAMESPACE)})""")
    con.execute("""DELETE FROM ap_ontology_name_candidate c
        USING ap_ontology_original_named n WHERE n.entity_key=c.entity_key""")
    # An accession presented as a scalar label must not become a Name alias.
    # This check is global across source occurrences, not just the owning edge.
    con.execute("""CREATE OR REPLACE TABLE ap_ontology_original_label_value AS
        SELECT DISTINCT i.entity_key,i.item.id AS value
        FROM ap_identifier_raw i SEMI JOIN ap_ontology_name_candidate c
             ON c.entity_key=i.entity_key
        WHERE i.item.id IS NOT NULL""")
    con.execute(f"""INSERT INTO ap_ontology_lookup_alias
        SELECT c.resource,c.version,c.entity_key,{name_namespace},c.label,
               {_uuid(_key(_literal("identifier"), name_namespace, "c.label"))}
        FROM ap_ontology_name_candidate c
        ANTI JOIN ap_ontology_original_label_value original
             ON original.entity_key=c.entity_key AND original.value=c.label""")
    # Deduplicate the small alias set, then insert only missing UUIDs. Do not
    # regroup the release-wide identifier dictionary to add these lookup rows.
    con.execute("""INSERT INTO ap_identifier
        SELECT aliases.alias_identifier_id,it.id,aliases.alias_value,NULL::VARCHAR
        FROM (SELECT DISTINCT alias_identifier_id,namespace_name,alias_value
              FROM ap_ontology_lookup_alias) aliases
        JOIN ap_vocab_identifier_type it ON it.name=aliases.namespace_name
        ANTI JOIN ap_identifier original
             ON original.identifier_id=aliases.alias_identifier_id""")
    for table in (
        "ap_ontology_lookup_owner",
        "ap_ontology_name_candidate",
        "ap_ontology_original_named",
        "ap_ontology_original_label_value",
    ):
        con.execute(f"DROP TABLE {table}")


def _create_annotation_projection(con: Any) -> None:
    con.execute(f"""CREATE OR REPLACE TABLE ap_annotation_occurrence AS SELECT a.*,
        {_uuid("to_json(struct_pack(kind:='annotation',term:=term,value:=value,quantity:=quantity))")}
            annotation_key
        FROM ap_annotation_raw a""")
    con.execute("""CREATE OR REPLACE TABLE ap_annotation AS SELECT DISTINCT
        annotation_key,term,
        CASE WHEN quantity.has_numeric_value IS NOT NULL THEN
             quantity.has_numeric_value::VARCHAR ELSE value END AS "value",
        quantity.has_unit unit
        FROM ap_annotation_occurrence WHERE term IS NOT NULL""")
    con.execute("""CREATE OR REPLACE TABLE ap_annotation_quantity AS SELECT DISTINCT
        annotation_key,quantity.has_numeric_value has_numeric_value,
        quantity.has_unit has_unit,quantity.has_unit_prefix has_unit_prefix,
        quantity.has_binary_relation has_binary_relation,quantity.source_field source_field,
        quantity.comparator comparator,value published_value
        FROM ap_annotation_occurrence WHERE quantity IS NOT NULL""")
    _checked(
        con,
        "SELECT annotation_key FROM ap_annotation GROUP BY annotation_key HAVING count(*)>1",
        "Annotation dictionary UUID collision",
    )
    _checked(
        con,
        "SELECT annotation_key FROM ap_annotation_quantity GROUP BY annotation_key HAVING count(*)>1",
        "Quantity dictionary UUID collision",
    )


def _pg_array(expression: str) -> str:
    # The COPY cell is PostgreSQL array input, including quoting/NULL elements.
    return f"""CASE WHEN {expression} IS NULL THEN NULL ELSE '{{' || coalesce(
        array_to_string(list_transform({expression}, x -> CASE WHEN x IS NULL THEN 'NULL'
            ELSE chr(34)||replace(replace(x,chr(92),chr(92)||chr(92)),
                 chr(34),chr(92)||chr(34))||chr(34) END), ','),'') || '}}' END"""


def _prepare_relation_evidence_annotation_copy(
    connection: Any,
    *,
    shard_count: int = 128,
    on_progress: Callable[[Mapping[str, Any]], None] | None = None,
) -> None:
    """Materialize exact annotation links without one release-wide anti-join.

    Persist only primitive link inputs, clustered by a complete owner's shard.
    Values, quantities and arrays already have their own tables and never enter
    these joins. A published evidence UUID belongs to one resource/version/owner,
    so four-column deduplication within an owner-complete shard is also global.
    """
    if isinstance(shard_count, bool) or not isinstance(shard_count, int) or shard_count < 1:
        raise ValueError("shard_count must be a positive integer")

    def report(stage: str, state: str, **details: Any) -> None:
        if on_progress is not None:
            on_progress({"stage": stage, "state": state, "shard_count": shard_count, **details})

    report("annotation_input", "start")
    connection.execute(f"""CREATE OR REPLACE TABLE ap_annotation_link_input AS
        SELECT resource,version,owner_key,evidence_ordinal,annotation_key,source,dataset,
            CASE WHEN scope='evidence' THEN 'relation'
                 ELSE coalesce(scope,'relation') END normalized_scope,
            owner_kind='evidence' is_evidence,term IS NOT NULL is_typed,
            (hash(resource,version,owner_key)%{shard_count})::INTEGER owner_shard
        FROM ap_annotation_occurrence
        WHERE owner_kind IN ('relation','evidence')
        ORDER BY owner_shard""")
    report("annotation_input", "done")
    report("evidence_input", "start")
    connection.execute(f"""CREATE OR REPLACE TABLE ap_evidence_link_input AS
        SELECT resource,version,relation_key owner_key,ordinal,source_id,
            relation_evidence_id,source,dataset,synthetic,
            (hash(resource,version,relation_key)%{shard_count})::INTEGER owner_shard
        FROM ap_evidence ORDER BY owner_shard""")
    report("evidence_input", "done")
    connection.execute("""CREATE OR REPLACE TABLE ap_copy_relation_evidence_annotation (
        source_id BIGINT,relation_evidence_id UUID,annotation_key UUID,
        annotation_scope_id SMALLINT)""")
    for shard in range(shard_count):
        report("shard", "start", shard=shard)
        connection.execute(f"""INSERT INTO ap_copy_relation_evidence_annotation
            WITH annotations AS MATERIALIZED (
                SELECT resource,version,owner_key,evidence_ordinal,annotation_key,
                    source,dataset,normalized_scope,is_evidence,is_typed
                FROM ap_annotation_link_input WHERE owner_shard={shard}
            ), evidence AS MATERIALIZED (
                SELECT resource,version,owner_key,ordinal,source_id,
                    relation_evidence_id,source,dataset,synthetic
                FROM ap_evidence_link_input WHERE owner_shard={shard}
            ), observed AS MATERIALIZED (
                SELECT DISTINCT resource,version,owner_key,annotation_key,
                    source,dataset,normalized_scope
                FROM annotations WHERE is_evidence
            ), true_statement_annotation AS (
                SELECT a.resource,a.version,a.owner_key,a.annotation_key,
                    a.source,a.dataset,a.normalized_scope
                FROM annotations a WHERE NOT a.is_evidence AND a.is_typed
                  AND NOT EXISTS (
                    SELECT 1 FROM observed o
                    WHERE o.resource=a.resource AND o.version=a.version
                      AND o.owner_key=a.owner_key AND o.annotation_key=a.annotation_key
                      AND o.source IS NOT DISTINCT FROM a.source
                      AND o.dataset IS NOT DISTINCT FROM a.dataset
                      AND o.normalized_scope=a.normalized_scope
                  )
            )
            SELECT e.source_id,e.relation_evidence_id,a.annotation_key,
                sc.id::SMALLINT annotation_scope_id
            FROM annotations a JOIN evidence e
              ON a.resource=e.resource AND a.version=e.version AND a.owner_key=e.owner_key
             AND a.evidence_ordinal=e.ordinal
            JOIN ap_vocab_annotation_scope sc ON sc.name=a.normalized_scope
            WHERE a.is_evidence AND a.is_typed
            UNION
            SELECT e.source_id,e.relation_evidence_id,a.annotation_key,
                sc.id::SMALLINT annotation_scope_id
            FROM true_statement_annotation a JOIN evidence e
              ON a.resource=e.resource AND a.version=e.version AND a.owner_key=e.owner_key
             AND (e.synthetic OR
                  ((a.source IS NULL OR a.source=e.source)
                   AND (a.dataset IS NULL OR a.dataset=e.dataset)))
            JOIN ap_vocab_annotation_scope sc ON sc.name=a.normalized_scope""")
        report("shard", "done", shard=shard)
    connection.execute("DROP TABLE ap_annotation_link_input")
    connection.execute("DROP TABLE ap_evidence_link_input")


def _queries(*, retain_published_provenance: bool = False) -> tuple[CopyQuery, ...]:
    result: list[CopyQuery] = []

    def add(table: str, columns: str, query: str) -> None:
        result.append(CopyQuery(table, tuple(columns.split()), query))

    add(
        "entity",
        "entity_id entity_type_id taxonomy_id canonical_identifier_type_id "
        "canonical_identifier resolution_status_id resolution_mechanism",
        "SELECT * FROM ap_entity",
    )
    add(
        "identifier_evidence",
        "identifier_id identifier_type_id value value_normalized",
        "SELECT * FROM ap_identifier",
    )
    add("annotation", "annotation_key term value unit", "SELECT * FROM ap_annotation")
    add(
        "entity_evidence",
        "source_id entity_evidence_id dataset_id row_id "
        "parent_entity_evidence_id entity_role_id entity_type_id taxonomy_id",
        """
        SELECT source_id,entity_evidence_id,dataset_id,row_id::BIGINT,NULL::UUID,
            1::SMALLINT,entity_type_id,try_cast(taxon AS BIGINT)
        FROM ap_entity_occurrence""",
    )
    add(
        "entity_evidence_resolution",
        "source_id entity_evidence_id status_id entity_id reason_id molecular_type_id",
        f"""SELECT source_id,entity_evidence_id,
            {PUBLISHED_STATUS_ID}::SMALLINT,entity_id,NULL::SMALLINT,NULL::SMALLINT
        FROM ap_entity_occurrence""",
    )
    add(
        "entity_identifier",
        "source_id entity_id identifier_id",
        f"""WITH published_links AS (
        SELECT DISTINCT e.source_id,e.entity_id,i.identifier_id
        FROM ap_identifier_occurrence i JOIN ap_entity_occurrence e
        USING(resource,version,entity_key) WHERE i.identifier_id IS NOT NULL
        UNION SELECT DISTINCT e.source_id,e.entity_id,
            {_uuid(_key(_literal("identifier"), "coalesce(a.name,e.namespace)", "e.identifier"))}
        FROM ap_entity_occurrence e LEFT JOIN ap_namespace_alias a ON a.raw=e.namespace
        WHERE {_generic_identifier_included("coalesce(a.name,e.namespace)", "e.identifier")}
        UNION SELECT DISTINCT occurrence.source_id,e.entity_id,
            {_uuid(_key(_literal("identifier"), "it.name", "e.canonical_identifier"))}
        FROM ap_entity_occurrence occurrence JOIN ap_entity e ON e.entity_id=occurrence.entity_id
        JOIN ap_vocab_identifier_type it ON it.id=e.canonical_identifier_type_id
        WHERE {_generic_identifier_included("it.name", "e.canonical_identifier")})
        SELECT * FROM published_links
        UNION SELECT links.source_id,links.entity_id,a.alias_identifier_id
        FROM published_links links JOIN ap_identifier_alias a
             ON a.original_identifier_id=links.identifier_id
        UNION SELECT e.source_id,e.entity_id,a.alias_identifier_id
        FROM ap_ontology_lookup_alias a JOIN ap_entity_occurrence e
             USING(resource,version,entity_key)""",
    )
    # Computed lookup aliases do not claim source evidence and consequently do
    # not acquire main's authority/reference role from an occurrence join.
    add(
        "entity_evidence_identifier",
        "source_id entity_evidence_id identifier_id",
        """
        SELECT DISTINCT e.source_id,e.entity_evidence_id,i.identifier_id
        FROM ap_identifier_occurrence i JOIN ap_entity_occurrence e
        USING(resource,version,entity_key) WHERE i.identifier_id IS NOT NULL""",
    )
    add(
        "entity_evidence_annotation",
        "source_id entity_evidence_id annotation_key",
        """
        SELECT DISTINCT e.source_id,e.entity_evidence_id,a.annotation_key
        FROM ap_annotation_occurrence a JOIN ap_entity_occurrence e
        ON a.resource=e.resource AND a.version=e.version AND a.owner_key=e.entity_key
        WHERE a.owner_kind='entity' AND a.term IS NOT NULL""",
    )
    add(
        "relation",
        "relation_id subject_entity_id predicate_id object_entity_id relation_category_id",
        "SELECT * FROM ap_relation",
    )
    add(
        "relation_evidence",
        "source_id relation_evidence_id dataset_id row_id "
        "subject_entity_evidence_id subject_entity_id predicate_id object_entity_evidence_id "
        "object_entity_id relation_category_id",
        """SELECT e.source_id,e.relation_evidence_id,
        e.dataset_id,e.row_id,NULL::UUID,r.subject_entity_id,r.predicate_id,NULL::UUID,
        r.object_entity_id,r.relation_category_id
        FROM ap_evidence e JOIN ap_statement r USING(resource,version,relation_key)""",
    )
    add(
        "relation_evidence_relation",
        "source_id relation_id relation_evidence_id",
        """
        SELECT e.source_id,r.relation_id,e.relation_evidence_id FROM ap_evidence e
        JOIN ap_statement r USING(resource,version,relation_key) WHERE r.relation_id IS NOT NULL""",
    )
    add(
        "relation_evidence_annotation",
        "source_id relation_evidence_id annotation_key annotation_scope_id",
        """WITH true_statement_annotation AS (
            SELECT a.* FROM ap_annotation_occurrence a
            WHERE a.owner_kind='relation' AND a.term IS NOT NULL AND NOT EXISTS (
                SELECT 1 FROM ap_annotation_occurrence observed
                WHERE observed.owner_kind='evidence' AND observed.resource=a.resource
                  AND observed.version=a.version AND observed.owner_key=a.owner_key
                  AND observed.annotation_key=a.annotation_key
                  AND observed.source IS NOT DISTINCT FROM a.source
                  AND observed.dataset IS NOT DISTINCT FROM a.dataset
                  AND (CASE WHEN observed.scope='evidence' THEN 'relation'
                            ELSE coalesce(observed.scope,'relation') END)
                    = (CASE WHEN a.scope='evidence' THEN 'relation'
                            ELSE coalesce(a.scope,'relation') END)
            )
        ), own_annotation AS (
            SELECT * FROM ap_annotation_occurrence WHERE owner_kind='evidence' AND term IS NOT NULL
            UNION ALL SELECT * FROM true_statement_annotation
        )
        SELECT DISTINCT e.source_id,e.relation_evidence_id,
        a.annotation_key,sc.id::SMALLINT
        FROM own_annotation a JOIN ap_evidence e
        ON a.resource=e.resource AND a.version=e.version AND a.owner_key=e.relation_key
        AND (a.owner_kind='relation' OR (a.owner_kind='evidence' AND a.evidence_ordinal=e.ordinal))
        AND (a.owner_kind='evidence' OR e.synthetic OR
             ((a.source IS NULL OR a.source=e.source) AND (a.dataset IS NULL OR a.dataset=e.dataset)))
        JOIN ap_vocab_annotation_scope sc ON sc.name=CASE WHEN a.scope='evidence' THEN 'relation'
                                                       ELSE coalesce(a.scope,'relation') END
        WHERE a.term IS NOT NULL""",
    )
    add(
        "entity_ontology_relation",
        "source_id subject_entity_id predicate_id object_entity_id ontology_id",
        """SELECT DISTINCT ds.id,r.subject_entity_id,r.predicate_id,r.object_entity_id,r.ontology_id
        FROM ap_statement r JOIN ap_data_source ds ON ds.name=r.resource
        WHERE r.statement_kind IN ('ontology','ontology_axiom')""",
    )
    # Current main's fresh EvidenceProjector never populates ontology_terms_raw.
    # Preserve its base table/COPY contract; the main shared derivation produces
    # serving entity_ontology_term from ontology edges, identifiers and annotations.
    add(
        "ontology_terms",
        "source_id term_entity_id term_id ontology_prefix label definition "
        "ontology_id synonyms synonyms_text sources",
        """SELECT
        NULL::BIGINT AS source_id,NULL::UUID AS term_entity_id,NULL::VARCHAR AS term_id,
        NULL::VARCHAR AS ontology_prefix,NULL::VARCHAR AS label,NULL::VARCHAR AS definition,
        NULL::VARCHAR AS ontology_id,NULL::VARCHAR AS synonyms,NULL::VARCHAR AS synonyms_text,
        NULL::VARCHAR AS sources WHERE false""",
    )
    # Reference keys are published identities, including virtual gene anchors.
    # No resolution, target-row requirement, or state interpretation belongs here.
    add(
        "entity_reference_context",
        "resource version entity_key entity_id reference_entity_key gene_reference_keys",
        """SELECT resource,version,entity_key,entity_id,reference_entity_key,
        CASE WHEN gene_reference_keys IS NULL THEN NULL ELSE to_json(gene_reference_keys)::VARCHAR END
        FROM ap_entity_occurrence""",
    )
    add(
        "statement_reference_context",
        "resource version relation_key subject_reference_entity_key object_reference_entity_key",
        """SELECT resource,version,relation_key,subject_reference_entity_key,
        object_reference_entity_key FROM ap_statement_raw""",
    )
    add(
        "molecular_evidence_context",
        "resource version owner_kind owner_key ordinal relation_evidence_id "
        "source dataset row_id upstream_id occurrence_json",
        """SELECT r.resource,r.version,'relation',r.relation_key,r.ordinal,
        e.relation_evidence_id,r.item.source,r.item.dataset,r.item.row_id,r.item.upstream_id,
        to_json(r.item)::VARCHAR FROM ap_evidence_raw r
        JOIN ap_evidence e USING(resource,version,relation_key,ordinal)
        UNION ALL
        SELECT resource,version,'entity',entity_key,ordinal,NULL::UUID,
        item.source,item.dataset,item.row_id,item.upstream_id,to_json(item)::VARCHAR
        FROM ap_entity_evidence_raw""",
    )
    if retain_published_provenance:
        add(
            "parquet_entity",
            "resource version entity_key entity_id entity_evidence_id label namespace "
            "identifier entity_type taxon has_hierarchy parent_count child_count identifiers_present "
            "annotations_present",
            """SELECT resource,version,entity_key,entity_id,entity_evidence_id,
            label,namespace,identifier,entity_type,taxon,has_hierarchy,parent_count,child_count,
            identifier_count>0,annotation_count>0 FROM ap_entity_occurrence""",
        )
        add(
            "parquet_statement",
            "resource version relation_key relation_id statement_kind "
            "subject_entity_id object_entity_id predicate_id subject_label subject_type object_label "
            "object_type taxon is_directed sign category interaction_class evidence_count sources "
            "evidence_present annotations_present",
            f"""SELECT resource,version,relation_key,relation_id,
            statement_kind,subject_entity_id,object_entity_id,predicate_id,subject_label,subject_type,
            object_label,object_type,taxon,is_directed,sign,category,interaction_class,evidence_count,
            {_pg_array("sources")},EXISTS (SELECT 1 FROM ap_evidence_raw e WHERE
                e.resource=s.resource AND e.version=s.version AND e.relation_key=s.relation_key),
            annotation_count>0 FROM ap_statement s""",
        )
        add(
            "parquet_evidence",
            "resource version relation_key ordinal source_id relation_evidence_id "
            "source dataset original_row_id upstream_id annotations_present synthetic",
            """SELECT resource,version,relation_key,ordinal,source_id,relation_evidence_id,
            source,dataset,original_row_id,upstream_id,annotations_present,synthetic FROM ap_evidence""",
        )
        add(
            "parquet_identifier_occurrence",
            "resource version entity_key ordinal identifier_id ns identifier is_canonical source",
            """SELECT resource,version,entity_key,ordinal,
            identifier_id,item.ns,item.id,item.is_canonical,item.source FROM ap_identifier_occurrence""",
        )
        add(
            "parquet_annotation_occurrence",
            "resource version owner_kind owner_key evidence_ordinal "
            "ordinal annotation_key term untyped_value source dataset scope",
            """
            SELECT resource,version,owner_kind,owner_key,evidence_ordinal,ordinal,annotation_key,
            term,CASE WHEN term IS NULL THEN value END,source,dataset,scope FROM ap_annotation_occurrence""",
        )
    add(
        "annotation_quantity",
        "annotation_key has_numeric_value has_unit has_unit_prefix "
        "has_binary_relation source_field comparator published_value",
        "SELECT * FROM ap_annotation_quantity",
    )
    return tuple(result)


def prepare_aligned_release(
    connection: Any,
    release: PinnedRelease | Iterable[ResourceSelection],
    *,
    dimension_rows: Mapping[str, Iterable[tuple]] | None = None,
    progress: Callable[[Mapping[str, Any]], None] | None = None,
    retain_published_provenance: bool = False,
) -> AlignedCopyPlan:
    """Prepare streaming COPY SELECTs, keeping main's strict identity contracts.

    ``dimension_rows`` contains existing ``(id,name)`` rows by main vocabulary
    table; dataset rows are ``(dataset_id,source_id,name)``. Returned dimensions
    include those seed rows; the caller inserts only new rows and advances its
    sequences. Data sources in the result tell it which source partitions to
    create. ``counts`` are exact prepared COPY row counts, not estimates.

    By default, exact published crosswalks and repeated array occurrences stay
    in the pinned Parquets. ``retain_published_provenance=True`` additionally
    prepares the five PostgreSQL audit tables; normalized scientific tables,
    mandatory molecular context and quantity metadata are identical in both modes.

    SQL work is cancellable through DuckDB. Nothing mutates a PostgreSQL schema
    or an input artifact. Full raw payload Parquets are never scanned here.
    """
    if type(retain_published_provenance) is not bool:
        raise ValueError("retain_published_provenance must be a boolean")

    def report(phase: str, **details: Any) -> None:
        if progress is not None:
            progress({"phase": phase, **details})

    resources = _resources(release)
    report("validate_and_stage", resources=len(resources))
    _create_inputs(connection, resources)
    report("validate_canonical_inputs")
    _validate_canonical_inputs(connection)
    report("dimensions")
    dimensions = _prepare_dimensions(connection, dimension_rows or {})
    report("entities")
    _create_entity_projection(connection)
    report("statements_and_evidence")
    _create_statement_projection(connection)
    report("ontology_lookup_aliases")
    _create_ontology_lookup_aliases(connection)
    report("annotation_dictionary")
    _create_annotation_projection(connection)
    # Materialize each final SELECT once. DISTINCT link reductions and shared
    # joins must not execute independently for counting and PostgreSQL COPY.
    # Caller-owned on-disk DuckDB keeps this bounded by its memory limit.
    queries = []
    counts = {}
    for item in _queries(retain_published_provenance=retain_published_provenance):
        report("prepare_copy", table=item.table)
        names = ",".join('"' + column + '"' for column in item.columns)
        staged = f"ap_copy_{item.table}"
        if item.table == "relation_evidence_annotation":
            _prepare_relation_evidence_annotation_copy(
                connection,
                on_progress=lambda details: report("prepare_annotation_links", **details),
            )
        else:
            connection.execute(f"CREATE OR REPLACE TABLE {staged} ({names}) AS {item.query}")
        counts[item.table] = int(connection.execute(f"SELECT count(*) FROM {staged}").fetchone()[0])
        queries.append(CopyQuery(item.table, item.columns, f"SELECT {names} FROM {staged}"))
        report("copy_prepared", table=item.table, rows=counts[item.table])
    return AlignedCopyPlan(
        tuple(queries),
        dimensions,
        counts,
        dimensions["data_source"],
        {
            "resolution": "published identities; no original matched flag inferred",
            "entity_evidence": "one explicit aggregate per published resource/entity",
            "original_entity_occurrences_available": False,
            "original_resolution_diagnostics_available": False,
            "relation_identity": (
                "main triple; published statement crosswalk retained in PostgreSQL"
                if retain_published_provenance
                else "main triple; published statement identity retained in pinned Parquets"
            ),
            "molecular_context_location": "postgresql_and_pinned_parquet",
            "retain_published_provenance": retain_published_provenance,
            "published_provenance_location": (
                "postgresql_and_pinned_parquet" if retain_published_provenance else "pinned_parquet"
            ),
            "published_input_counts": {
                name: int(connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0])
                for name, table in (
                    ("entity_rows", "ap_entity_raw"),
                    ("statement_rows", "ap_statement_raw"),
                    ("evidence_occurrences", "ap_evidence_raw"),
                    ("identifier_occurrences", "ap_identifier_raw"),
                    ("annotation_occurrences", "ap_annotation_raw"),
                )
            },
            "raw_payloads_loaded": False,
            "record_json_loaded": False,
            "generic_identifier_policy": "main_standard_inchi_excluded",
            "excluded_standard_inchi_identifier_occurrences": int(
                connection.execute(
                    "SELECT count(*) FROM ap_identifier_raw WHERE item.ns IN ('inchi',?) "
                    "AND item.id IS NOT NULL",
                    [STANDARD_INCHI_NAMESPACE],
                ).fetchone()[0]
            ),
            "excluded_standard_inchi_dictionary_values": int(
                connection.execute(
                    "SELECT count(*) FROM ap_identifier_input WHERE ns IN ('inchi',?) "
                    "AND identifier IS NOT NULL",
                    [STANDARD_INCHI_NAMESPACE],
                ).fetchone()[0]
            ),
            "safe_bare_identifier_alias_pairs": int(
                connection.execute("SELECT count(*) FROM ap_identifier_alias").fetchone()[0]
            ),
            "ontology_cv_lookup_alias_pairs": int(
                connection.execute(
                    "SELECT count(*) FROM ap_ontology_lookup_alias WHERE namespace_name=?",
                    [ONTOLOGY_CV_NAMESPACE],
                ).fetchone()[0]
            ),
            "ontology_published_name_lookup_alias_pairs": int(
                connection.execute(
                    "SELECT count(*) FROM ap_ontology_lookup_alias WHERE namespace_name=?",
                    [ONTOLOGY_NAME_NAMESPACE],
                ).fetchone()[0]
            ),
            "unknown_ontology_scope_statements": int(
                connection.execute(
                    "SELECT count(*) FROM ap_statement WHERE statement_kind IN ('ontology','ontology_axiom') "
                    "AND NOT ontology_scope_recovered"
                ).fetchone()[0]
            ),
            "source_scoped_name_fallback_entities": int(
                connection.execute(
                    "SELECT count(*) FROM ap_entity e JOIN ap_vocab_identifier_type it "
                    "ON it.id=e.canonical_identifier_type_id WHERE it.name=?",
                    [PUBLISHED_FALLBACK_NAMESPACE],
                ).fetchone()[0]
            ),
            "shared_identity_taxonomy_conflicts": int(
                connection.execute(
                    "SELECT count(*) FROM (SELECT entity_key FROM ap_entity_raw GROUP BY entity_key "
                    "HAVING min(nullif(taxon,''))<>max(nullif(taxon,''))) q"
                ).fetchone()[0]
            ),
            "shared_identity_taxonomy_null_and_known": int(
                connection.execute(
                    "SELECT count(*) FROM (SELECT entity_key FROM ap_entity_raw GROUP BY entity_key "
                    "HAVING min(nullif(taxon,''))=max(nullif(taxon,'')) AND bool_or(nullif(taxon,'') IS NULL)) q"
                ).fetchone()[0]
            ),
            "published_resolution_status_id": PUBLISHED_STATUS_ID,
        },
    )
