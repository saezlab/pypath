from collections.abc import Generator

import re
import pathlib as pl
import collections

import pandas as pd
import pyparsing

from pypath_common import _misc as _common
from pypath.share import curl
from pypath.utils import mapping, taxonomy

REORGANISM = re.compile(
    r'#(\d+)# '
    r'((?:no activity in )?)'
    r'([-\w\s\.\[\]]+[^\s#\{<\(]).*'
)
REEC = re.compile(r'EC ([\d\.]+)')
REID = re.compile(r'\{([\w\.]+); source: (\w+)\}')
REISOFORM = re.compile('isoform ([\d\w/ ]+), cf\. EC ([\d\.]+)')
REEFFECT = re.compile(
    r'#([\d,]+)# '     # proteins (by numeric reference)
    r'([^\(][-\w\+\s,%\.]+) ' # compound name
    r' ?\(?((?:.*)?)\)? ' # within parentheses concentration & time
    r'<([\d,]+)>'  # literature references (numeric)
)
REKIKM = re.compile(
    r'#([\d,]+)# '  # proteins (by numeric reference)
    r'([-\d\.]+) '          # concentration
    r'\{([-\w\+\s,%\.\(\)/\'\[\]]+)\} '  # compound name
    r' ?\(?((?:.*)?)(?:>\))? ' # within parentheses concentration & time
    r'<([\d,]+)>'  # literature references (numeric)
)
REREFE = re.compile(
    r'<([\d]+)>'  # literature references (numeric)
    r'.*'
    r'Pubmed:(\d+)'
)


ALLOSTERIC_ROLES = {
    'IN': 'inhibitor',
    'AC': 'activator',
    'CF': 'cofactor',
}

RECORDS_ENABLED = {'ID', 'PR', 'AC', 'IN', 'CF', 'KI', 'KM', 'RF'}

# Example of an effects line:
#7,13,20,101,118,127,153,178# Cu2+ (#20# no effect <75>; #101# 1 mM, 99% loss of activity <168>; #127# 1 mM, 89% of initial activity <218>; #153# 1 mM, no% residual activity <243>; #7# 100 mM, 76% of initial activity <169>; #178# over 90% inhibition at 0.25 mM <307>; #118# 1 mM, 44.2% of initial activity <308>) <45,75,168,169,218,243,307,308>

AllostericRegulation = collections.namedtuple(
    'AllostericRegulation',
    [
        'action',
        'compound',
        'organism',
        'protein',
        'id_type',
        'wrong_ec',
        'pubmeds',
        'reaction_constants',
    ]
)

ReactionConstant = collections.namedtuple(
    'ReactionConstant',
    [
        'type',
        'value',
        'conditions',
        'pubmeds',
    ]
)


def main(
        output_dir: str = 'brenda_output',
        limit: int = 20,
    ) -> Generator[tuple[str, str], None, None]:

    #output of directory
    output_dir = pl.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outfile0 = output_dir / 'brenda.xlsx'
    outfile1 = output_dir / 'brenda-translated.xlsx'
    # Download data from BRENDA
    url = "https://www.brenda-enzymes.org/download.php"
    c = curl.Curl(
        url,
        post = {'dlfile': 'dl-textfile', 'accept-license': '1'},
        large = True,
        compr = 'tgz',
    )
    infile = _common.first(c.result.values())

    # Prepare for the loop
    yielded_records = 0
    current_record = []

    # Clean data, extract activator and inhibitor for each enzyme
    for ln in infile:

        if (
            not (ln := ln.strip(b'\r\n')) or
            ln.startswith(b'*') or
            b'\t' not in ln
        ):

            continue

        try:
            label, line = ln.decode('utf-8').split('\t', maxsplit = 1)

        except:
            print(ln.decode('utf-8').split('\t'))
            raise

        if label:

            if current_record:

                yielded_records += 1
                yield current_label, "".join(current_record)

                if yielded_records == limit:

                    break

            current_record = []
            current_label = label

        current_record.append(line)


def allosteric_regulation(
        organisms: str | int | list[str | int] = 'mouse',
        limit: int = 20,
    ) -> Generator[tuple[str, str], None, None]:

    organisms = _common.to_list(organisms)
    organisms = {taxonomy.ensure_latin_name(o) for o in organisms}
    record = None

    i = 0
    for ln in main(limit = None):

        label, data = ln

        if label not in RECORDS_ENABLED:

            continue

        if label == 'ID':

            if record:
                i += 1
                if i == limit:
                    break
                yield from allosteric_regulation_records(record)

            ec = 'ec:' + data.strip(',')

            record = {
                'ec': ec,
                'proteins': {},
                'actions': [],
                'km_ki': [],
                'references': {}
            }

        elif label == 'PR':

            pr_match = REORGANISM.match(data).groups()
            org_id, negation, organism = REORGANISM.match(data).groups()

            if organism not in organisms or negation:
                continue

            ecs = REEC.findall(data)
            ids = REID.findall(data)
            isoforms = REISOFORM.findall(data)
            record['proteins'][org_id] = (organism, ids, ecs, isoforms)

        elif label in {'IN', 'AC', 'CF'}:

            role = ALLOSTERIC_ROLES[label]
            effects = REEFFECT.findall(data)
            record['actions'].append((role, effects))


        elif label in {'KI', 'KM'}:

            values = REKIKM.findall(data)
            record['km_ki'].append((label, _common.first(values)))

        elif label == 'RF':

            if ref_match := REREFE.match(data):

                ref, pub = ref_match.groups()
                record['references'][ref] = pub

    yield from allosteric_regulation_records(record)


def allosteric_regulation_records(
        stage0: dict,
    ) -> Generator[AllostericRegulation, None, None]:

    def collect_refs(ref_idx: str) -> list[str]:

        return ';'.join(
            stage0['references'][i]
            for i in ref_idx.split(',')
            if i in stage0['references']
        )

    for actions in stage0['actions']:

        reactions_constants = collections.defaultdict(list)

        for k, (k_proteins, k_value, compound, k_details, refs) in stage0['km_ki']:

            for k_prot_idx in k_proteins.split(','):

                try:
                    k_value = float(k_value)
                except ValueError: pass
                k_pubmeds = collect_refs(refs)
                reactions_constants[(k_prot_idx, compound)].append(
                    ReactionConstant(k, k_value, k_details, k_pubmeds)
                )

        for action in actions[1]:

            protein_indices = (
                set(action[0].split(',')) &
                set(stage0['proteins'].keys())
            )

            for protein_idx in protein_indices:

                protein_data = stage0['proteins'][protein_idx]
                protein_ids = None
                wrong_ec = None
                id_type = None

                if uniprots := protein_data[1]:

                    protein_ids = [u[0] for u in uniprots]
                    id_type = 'uniprot'

                elif genesymbols := protein_data[3]:

                    wrong_ec = protein_data[2]
                    protein_ids = [g[0] for g in genesymbols]
                    id_type = 'genesymbol'

                pubmed = collect_refs(action[3])
                compound = action[1]

                yield AllostericRegulation(
                    action = actions[0],
                    compound = compound,
                    organism = protein_data[0],
                    protein = protein_ids,
                    id_type = id_type,
                    wrong_ec = wrong_ec,
                    pubmeds = pubmed,
                    reaction_constants = reactions_constants[
                        (protein_idx, compound)
                    ]
                )


def rest():
    # Clean EC with ()
    df_new = df_new[~df_new['ENTRY'].str.contains(r'\(.*\)', na=False)]
    df_new = df_new[1:].reset_index(drop=True)

    df_new.to_excel(outfile0)

    # ID conversion
    tmp = []


    def ID_conversion(df, column_name, label, organism):
        for _, row in df.iterrows():
            entry = row['ENTRY']
            uniprot_id = row['Uniprot']
            act_value = row[column_name]

            if pd.isna(act_value):
                continue

            metabolites = act_value.split(';')
            for name in metabolites:
                name = name.strip()
                if name.lower() == "more":
                    continue
                if name:
                    metabolite_id = mapping.map_name(name, 'iupac', 'pubchem')
                    metabolite_id = ';'.join(str(x) for x in metabolite_id if x is not None)
                    tmp.append([entry, uniprot_id, name, metabolite_id, label, organism])


        ID_conversion(df_new, "Act", "Activator", Name_list[0])
        ID_conversion(df_new, "Inh", "Inhibitor", Name_list[0])
        final_df = pd.DataFrame(tmp, columns=[
        "ENTRY", "Uniprot", "metabolite name", "Pubchem ID", "Act/Inh", "Organism"])
        final_df.to_excel(outfile1, index=False)
        return final_df

