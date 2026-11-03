"""Tiny hub Parquet files with one case per identity rule (see test_identity_build.py)."""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA = pa.schema(
    [(n, pa.string()) for n in ("source_type", "source_id", "hub_id", "taxonomy_id", "backend")]
)


def key(letter: str) -> str:
    """A syntactically valid full InChIKey."""
    return f"{letter * 14}-{letter * 10}-N"


class Hub:
    def __init__(self, name, taxon="0"):
        self.name, self.taxon, self.rows = name, taxon, []

    def record(self, hub_id, *pairs, taxon=None):
        """One record: its identity row plus (source_type, value) pairs."""
        for source_type, value in ((self.name, hub_id), *pairs):
            self.rows.append((source_type, value, hub_id, taxon or self.taxon, self.name))
        return self

    def write(self, directory: Path):
        columns = list(zip(*self.rows)) if self.rows else [[]] * 5
        table = pa.table(dict(zip(SCHEMA.names, columns)), schema=SCHEMA)
        pq.write_table(table, directory / f"{self.name}.parquet")


def chemical_hubs():
    chebi, pubchem, hmdb = Hub("chebi"), Hub("pubchem"), Hub("hmdb")
    chembl, bigg, metanetx = Hub("chembl"), Hub("bigg"), Hub("metanetx")
    lipidmaps, swisslipids, ramp = Hub("lipidmaps"), Hub("swisslipids"), Hub("ramp")

    # Rule 1: InChIKey anchors; two records share one structure; secondary ChEBI id; cas; kegg.
    chebi.record(
        "CHEBI:1",
        ("inchikey", key("A")),
        ("name", "Alpha"),
        ("synonym", "alpha synonym"),
        ("cas", "50-00-0"),
        ("kegg", "C00001"),
        ("chebi", "CHEBI:99"),
        ("smiles", "CC"),
        ("formula", "C2H6"),
    )
    chebi.record("CHEBI:3", ("inchikey", key("D")), ("name", "Delta"), ("metanetx", "MNXM3"))
    pubchem.record("100", ("inchikey", key("A")), ("smiles", "CC"))
    # An excluded "empty" InChIKey is not an anchor.
    pubchem.record("101", ("inchikey", "MOSFIJXAXDLOML-UHFFFAOYSA-N"))
    hmdb.record(
        "HMDB4", ("inchikey", key("D")), ("name", "Delta from HMDB"), ("drugbank", "DB00001")
    )

    # Rule 2: two InChIKeys quarantine the record; its keys still get their own entities.
    hmdb.record("HMDB2", ("inchikey", key("B")), ("inchikey", key("C")), ("name", "Twin"))
    # Two lipid names at the most specific level quarantine a record without InChIKey.
    lipidmaps.record("LMQ", ("name", "PE 16:0/18:1"), ("name", "PE 18:1/16:0"))

    # Rule 1 (lipid name): records without InChIKey are anchored by their most specific Goslin name.
    swisslipids.record("SLM:1", ("lipid_shorthand", "PC 16:0/18:1"), ("name", "PC 34:1"))
    lipidmaps.record("LM1", ("name", "PC 16:0/18:1"))

    # Rule 3: one-step attach, both directions; chains are not followed; disagreement blocks.
    chembl.record("CHEMBL1", ("chebi", "CHEBI:1"))  # outgoing edge to an anchored record
    metanetx.record("MNXM3")  # incoming edge: CHEBI:3 names it
    metanetx.record("MNXM4", ("metanetx", "MNXM3"))  # only reaches an anchorless record: no chain
    bigg.record("b1", ("chebi", "CHEBI:1"), ("chebi", "CHEBI:3"))  # two different anchors
    bigg.record("b2", ("hmdb", "HMDB2"))  # reaches a quarantined record
    bigg.record("b3", ("bigg", "b1"))  # a bigg claim for b1's id: fallback against b1's native row

    # Rule 4: group anchorless records; one record per hub, else every record stays alone.
    metanetx.record("MNXM10", ("kegg", "C00009"))
    bigg.record("g10", ("metanetx", "MNXM10"))
    metanetx.record("MNXM20", ("ramp", "RAMP_C_1"))
    metanetx.record("MNXM21", ("ramp", "RAMP_C_1"))
    ramp.record("RAMP_C_1")

    # Rule 5: a full-structure Goslin name resolves to an InChIKey only if the keys agree.
    chebi.record("CHEBI:5", ("inchikey", key("E")), ("name", "FA 18:0;5Me"))
    lipidmaps.record("LM5", ("inchikey", key("E")), ("name", "FA 18:0;5Me"))
    swisslipids.record("SLM:5", ("lipid_shorthand", "FA 18:0;5Me"))
    chebi.record("CHEBI:6", ("inchikey", key("F")), ("name", "FA 20:0;5Me"))
    lipidmaps.record("LM6", ("inchikey", key("G")), ("name", "FA 20:0;5Me"))
    # A species-level name never resolves to a structure.
    chebi.record("CHEBI:7", ("inchikey", key("H")), ("name", "PC 34:1"))
    return [chebi, pubchem, hmdb, chembl, bigg, metanetx, lipidmaps, swisslipids, ramp]


def protein_hubs():
    uniprot, entrez, ramp_gene = Hub("uniprot"), Hub("entrez"), Hub("ramp_gene")
    human = dict(taxon="9606")

    entrez.record(
        "7157",
        ("ensg", "ENSG00000141510.17"),
        ("enst", "ENST00000269305.9"),
        ("refseq", "NM_000546.6"),
        ("hgnc", "HGNC:11998"),
        ("genesymbol", "TP53"),
        **human,
    )
    entrez.record("7158", ("ensg", "ENSG0000007158"), **human)
    entrez.record("7159", **human)
    entrez.record("9999", **dict(taxon="10090"))
    uniprot.record(
        "P04637",
        ("uniprot", "P04637"),
        ("uniprot_entry", "P53_HUMAN"),
        ("genesymbol", "TP53"),
        ("genesymbol-syn", "P53"),
        ("ensg", "ENSG00000141510.17"),
        ("hgnc", "HGNC:11998"),
        ("ensp", "ENSP00000269305.4"),
        ("refseq_protein", "NP_000537.3"),
        ("genbank", "AAA12345.1"),
        ("entrez", "7157"),
        ("uniprot-sec", "Q15086"),
        ("uniprot-sec", "Q16535"),
        **human,
    )
    uniprot.record("P04637-2", ("uniprot", "P04637-2"), ("enst", "ENST00000999999"), **human)
    uniprot.record(
        "P99999",
        ("uniprot", "P99999"),
        ("uniprot_entry", "CYC_HUMAN"),
        ("genesymbol", "CYCS"),
        ("genesymbol-syn", "TP53"),
        ("hgnc", "HGNC:19986"),
        ("entrez", "7158"),
        ("entrez", "7159"),
        ("uniprot-sec", "Q16535"),
        **human,
    )
    uniprot.record("A0A000", ("uniprot", "A0A000"), ("uniprot_entry", "A0A000_HUMAN"), **human)
    # Known and different organisms: no gene-product link; an unknown gene: link without record.
    uniprot.record(
        "Q00002", ("uniprot", "Q00002"), ("entrez", "9999"), ("genesymbol", "MOUSEY"), **human
    )
    uniprot.record(
        "Q00003", ("uniprot", "Q00003"), ("entrez", "5555"), ("hgnc", "HGNC:5555"), **human
    )
    # Several genes: the one whose symbol is the protein's gene name keeps it (C4A, not
    # C4B); a protein naming two of its three genes keeps all three.
    for gene, symbol in (("720", "C4A"), ("721", "C4B"), ("1082", "CGB3"), ("93659", "CGB5"), ("94115", "CGB8")):
        entrez.record(gene, ("genesymbol", symbol), **human)
    uniprot.record(
        "P0C0L4", ("uniprot", "P0C0L4"), ("genesymbol", "C4A"), ("entrez", "720"), ("entrez", "721"), **human
    )
    uniprot.record(
        "P0DN86",
        ("uniprot", "P0DN86"),
        ("genesymbol", "CGB3"),
        ("genesymbol", "CGB5"),
        *(("entrez", gene) for gene in ("1082", "93659", "94115")),
        **human,
    )
    # More than ten genes behind one symbol: every gene stays a candidate.
    for n in range(12):
        entrez.record(str(100 + n), **human)
        uniprot.record(
            f"M{n:05d}",
            ("uniprot", f"M{n:05d}"),
            ("uniprot_entry", f"MANY{n}_HUMAN"),
            ("genesymbol", "MANYSYM"),
            ("entrez", str(100 + n)),
            **human,
        )

    ramp_gene.record("RAMP_G_1", ("uniprot", "P04637"))
    ramp_gene.record("RAMP_G_2", ("entrez", "7157"))
    ramp_gene.record("RAMP_G_3", ("name", "lonely"))
    ramp_gene.record("RAMP_G_4", ("uniprot", "P99999"))
    ramp_gene.record("RAMP_G_5", ("uniprot", "P04637"), ("entrez", "7157"))
    return [uniprot, entrez, ramp_gene]


EQUATION = "pentanamide + H2O = pentanoate + NH4(+)"


def reaction_hubs():
    rhea, mnx = Hub("rhea"), Hub("metanetx_reaction")
    # Master reactions are the anchors; directional ids are claims of the master's record.
    rhea.record(
        "10000",
        ("rhea", "10001"),
        ("rhea", "10002"),
        ("rhea", "10003"),
        ("kegg_reaction", "R00001"),
        ("reactome", "R-HSA-1.3"),  # Rhea writes Reactome ids with a version
        ("ec", "3.5.1.50"),  # an attribute, never an identity
        ("name", EQUATION),
    )
    rhea.record("20000", ("kegg_reaction", "R00002"), ("metacyc_reaction", "RXN-1"))
    # MetaNetX reactions: attach to the one Rhea master they name.
    mnx.record(
        "MNXR1", ("rhea", "10000"), ("bigg_reaction", "AMIDASE"), ("kegg_reaction", "R00001")
    )
    # No Rhea reaction: a model reaction keeps MetaNetX's reconciled identity.
    mnx.record("MNXR2", ("bigg_reaction", "MODELONLY"), ("vmh_reaction", "MODELONLY"))
    # Same KEGG id as Rhea 20000 but no Rhea link: Rhea's own cross-reference wins.
    mnx.record("MNXR3", ("kegg_reaction", "R00002"))
    # Two Rhea masters: not merged, evidence for both.
    mnx.record("MNXR4", ("rhea", "10000"), ("rhea", "20000"), ("bigg_reaction", "TWOWAY"))
    return [rhea, mnx]


def write_hubs(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    for hub in chemical_hubs() + protein_hubs() + reaction_hubs():
        hub.write(directory)
    return directory
