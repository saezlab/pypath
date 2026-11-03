from pypath.inputs_v2 import chembl
from biolink_model.datamodel.model import Protein
from omnipath_build.silver import SilverExtractor
from omnipath_resolver.canonical.policy import get_policy
from omnipath_resolver.canonical.match import votes_for


def test_structured_peptide_is_chemical_but_target_stays_protein():
    row = {
        "chembl_id": "CHEMBL1",
        "molecule_type": "Protein",
        "standard_inchi_key": "CFOZOOXNCIABHM-FKBYEOEOSA-N",
    }
    mapped = chembl.molecules_schema(row)
    ex = SilverExtractor("chembl", "molecules")
    ex.process_record(mapped, row, "molecules:0", 0)
    obs = next(iter(ex.entities.values()))
    assert get_policy(obs.entity_type).library == "chemical"
    assert any(v.ns == "inchikey" for v in votes_for(obs, get_policy(obs.entity_type))[0])
    assert chembl.TARGET_TYPE_MAP["SINGLE PROTEIN"] == Protein
    assert chembl._molecule_entity_type({"molecule_type": "Protein"}) == Protein
