"""Small release-specific presentation cache for the explorer landing page."""

import logging
import hashlib
import json
import os
import tempfile
from pathlib import Path


logger = logging.getLogger(__name__)

# One example per kind of entity, in landing-page order, by identifier: labels are not
# unique ('lactate' also names a wax ester). Examples missing from a release are skipped.
EXAMPLES = (
    "entrez:7157",  # TP53 (gene, with its products)
    "CHEBI:15422",  # ATP
    "R-HSA-70171",  # Glycolysis (Reactome pathway)
    "MONDO:0005148",  # type 2 diabetes mellitus
    "entrez:3630",  # insulin
    "CHEBI:17234",  # glucose
    "LMGP01010005",  # POPC (a phosphatidylcholine)
    "GO:0006915",  # apoptotic process
    "MIMAT0000076",  # hsa-miR-21-5p
    "entrez:1956",  # EGFR
    "CHEBI:16113",  # cholesterol
    "HP:0000822",  # hypertension (phenotype)
    "CHEBI:15365",  # aspirin
    "RHEA:17825",  # hexokinase: D-glucose + ATP = D-glucose 6-phosphate + ADP + H(+)
    "CHEBI:27732",  # caffeine
    # PD-1/PD-L1 complex: complexes have no public accession; this id hashes its members
    "complex:35d325c2447b9e812bae10301e71ec21ad17524f91f7c0f63ce2cca7f357da1a",
)


# Bumped when the shape of a cached example page changes (group cards), so a release's
# prepared page is rebuilt rather than served stale.
CACHE_FORMAT = "examples-v3"


def _example(queries, token):
    """The entity an example identifier names, with its group key: the query
    example_entity."""
    text = queries.run_query("example_entity", {"token": token})
    return json.loads(text) if text else None


def get_examples(queries):
    """Landing-page examples, as entities and as the groups the grouped view shows."""
    from omnipath_api.queries.services.groups import search_groups

    infos = queries._selected_resource_infos()
    identity = [
        [info["key"] for info in infos],
        queries._inventory_fingerprint,
        EXAMPLES,
        CACHE_FORMAT,
        queries._release_scope.get(),
        queries._served().identity,
        [m.get("references") for m in queries.releases.list()],
    ]
    key = hashlib.sha256(json.dumps(identity, sort_keys=True, default=str).encode()).hexdigest()

    def load():
        target = queries.data_root / ".presentation" / CACHE_FORMAT / (key + ".json")
        try:
            return json.loads(target.read_text())
        except FileNotFoundError:
            pass
        except (OSError, ValueError):
            logger.warning("Cannot load example cache %s", target, exc_info=True)
        found = [found for token in (EXAMPLES if infos else ()) if (found := _example(queries, token))]
        found = list({item["entity"]["entityPk"]: item for item in found}.values())
        groups = []
        for group_key in dict.fromkeys(item["group_key"] for item in found):
            page = search_groups(queries, strategy="auto", group_key=group_key, member_limit=1)
            groups.extend(page["groups"][:1])
        result = dict(
            entities=[item["entity"] for item in found],
            groups=groups,
            nextCursor=None,
            total=len(found),
            kind="examples",
        )
        if queries.read_only:
            return result
        # A read-only data mount can still serve/cache examples in memory.
        staging = None
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(mode="w", dir=target.parent, delete=False) as handle:
                staging = Path(handle.name)
                json.dump(result, handle, default=str)
            os.replace(staging, target)
        except OSError:
            logger.warning("Cannot persist optional example cache %s", target, exc_info=True)
        finally:
            if staging is not None:
                staging.unlink(missing_ok=True)
        return result

    return queries._example_cache.get(key, load)
