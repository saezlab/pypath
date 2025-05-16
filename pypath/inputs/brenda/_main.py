from collections.abc import Generator

import re
import pathlib as pl
import collections

import pandas as pd

from pypath_common import _misc as _common
from pypath.share import curl
from pypath.utils import mapping, taxonomy

REORGANISM = re.compile(r'#(\d+)# ((?:no activity in )?)([-\w\s\.]+[^\s#\{]).*')
REEC = re.compile(r'EC ([\d\.]+)')
REID = re.compile(r'\{([\w\.]+); source: (\w+)\}')
REEFFECT = re.compile(
    r'((?:\()?)\s?#'
    r'([\d,]+)# '
    r'([^\(][\w\+\s,%\.]+)'
    r'((?: \([^\(\)]*\))?) '
    r'<([\d,]+)>'
)

AllostericRegulation = collections.namedtuple(
    'AllostericRegulation',
    [
        'ec',
        'uniprot',
        'regulator',
        'effect',
    ]
)


def main(
        organisms: str | int | list[str | int] = 'mouse',
        output_dir: str = 'brenda_output',
    ) -> Generator[tuple[str], None, None]:

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
    infile = c.result['brenda_2024_1.txt']

    # Prepare for the loop
    organisms = _common.to_list(organisms)
    organisms = {taxonomy.ensure_latin_name(o) for o in organisms}
    colums4now = ['EC','Uniprot','Act', 'Inh']
    df_new = pd.DataFrame(columns = colums4now)
    j = -1

    # Clean data, extract activator and inhibitor for each enzyme
    for ln in infile:

        if not (ln := ln.strip()):

            continue

        label, line = ln.split('\t')

        tmp_step = 1
        indent_of_the_line = ''

        while indent_of_the_line == '':

            if i == length - 1:
                break
            orig_nextline_splitted = re.split('\t', data[i + tmp_step])
            if len(body_now) < 2:
                body_now = re.sub('\n', '', body_now)
            else:
                if str.isdigit(body_now[-2]):
                    body_now = re.sub('\n', ',', body_now)
                else:
                    body_now = re.sub('\n', '', body_now)
            if not orig_nextline_splitted[0] == '':
                break
            body_now = body_now + ' ' + '/t'.join(orig_nextline_splitted[1:])
            tmp_step = tmp_step + 1


        if label == 'ID':

            record_organisms = {}
            j = j + 1
            df_new.loc[j, 'ENTRY'] = 'ec:' + line.strip(',')
            num_inh_thisentry = 0
            num_act_thisentry = 0
            id4species_here = []

        elif label == 'PR':

            org_id, negation, organism = REORGANISM.match(line).groups()[0]
            if organism not in organisms or negation:
                continue

            ecs = REEC.findall(line)
            ids = REID.findall(line)
            record_organisms[org_id] = (organism, ids, ecs)

            df_new.loc[j, 'Uniprot'] = ids
            df_new.loc[j, 'EC'] = ecs

        elif label in {'IN', 'AC'}:

            effect = 'Inh' if label == 'IN' else 'Act'
            effects = REEFFECT.findall(line)

            splitted_further_inh = re.split(r'# | <| [(]', line)
            prs_now = re.split(',', re.sub('#', '', splitted_further_inh[0]))
            IsThisINisOfthisspecies = False
            for k in range(0, len(prs_now)):
                for kk in range(0, len(id4species_here)):
                    if prs_now[k] == id4species_here[kk]:
                        IsThisINisOfthisspecies = True
            if IsThisINisOfthisspecies:
                if num_inh_thisentry == 0:
                    df_new.loc[j, effect] = splitted_further_inh[1]
                else:
                    df_new.loc[j, effect] = df_new.loc[j, effect] + ';' + splitted_further_inh[1]
                num_inh_thisentry = num_inh_thisentry + 1


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


# Example of an effects line:
#7,13,20,101,118,127,153,178# Cu2+ (#20# no effect <75>; #101# 1 mM, 99% loss of activity <168>; #127# 1 mM, 89% of initial activity <218>; #153# 1 mM, no% residual activity <243>; #7# 100 mM, 76% of initial activity <169>; #178# over 90% inhibition at 0.25 mM <307>; #118# 1 mM, 44.2% of initial activity <308>) <45,75,168,169,218,243,307,308>
