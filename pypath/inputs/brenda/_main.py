#
from pypath.share import curl
from pypath.utils import mapping
import pandas as pd
import re
import os
import numpy as np
import pathlib as pl


def main(output_dir = 'brenda_output'):
    #output of directory
    output_dir = pl.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outfile0 = output_dir / 'brenda.xlsx'
    outfile1 = output_dir / 'brenda-translated.xlsx'
    # Download data from BRENDA
    url = "https://www.brenda-enzymes.org/download.php"
    c = curl.Curl(url, post = {'dlfile': 'dl-textfile', 'accept-license': "1"}, large = True,compr = 'tgz')
    brenda_gen = c.result['brenda_2024_1.txt']
    data = [line.decode("utf-8") for line in brenda_gen]

    # Prepare for the loop
    Name_list = ['Mus musculus']
    length = len(data)
    colums4now = ['ENTRY','Uniprot','Act', 'Inh']
    df_new = pd.DataFrame(columns = colums4now)

    # Define function
    def contains_element(A, B):
        for element in B:
            if element in A:
                return True
        return False


    # Clean data, extract activator and inhibitor for each enzyme
    j = -1
    for i in range(0, length):
        body_now = data[i]
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
        splitted = re.split('\t', body_now)
        if splitted[0] == 'ID':
            j = j + 1
            df_new.loc[j, 'ENTRY'] = 'ec:' + splitted[1].strip(',')
            num_inh_thisentry = 0
            num_act_thisentry = 0
            id4species_here = []
        elif splitted[0] == 'PR':
            spname = splitted[1].split('#')[2].strip().split('  ')[0]
            if contains_element(spname, Name_list):
                if not 'no activity in' in spname:
                    splitted_further_pr = re.split(' ', splitted[1])
                    id4species_here.append(re.sub('#', '', splitted_further_pr[0]))
                    match = re.search(r'EC ([\d\.]+)', splitted[1])
                    uni_ID = None
                    if match:
                        content = match.group(1)
                        uni_ID = content.split(';')[0]
                    df_new.loc[j, 'Uniprot'] = uni_ID
        elif splitted[0] == 'IN':
            splitted_further_inh = re.split(r'# | <| [(]', splitted[1])
            prs_now = re.split(',', re.sub('#', '', splitted_further_inh[0]))
            IsThisINisOfthisspecies = False
            for k in range(0, len(prs_now)):
                for kk in range(0, len(id4species_here)):
                    if prs_now[k] == id4species_here[kk]:
                        IsThisINisOfthisspecies = True
            if IsThisINisOfthisspecies:
                if num_inh_thisentry == 0:
                    df_new.loc[j, 'Inh'] = splitted_further_inh[1]
                else:
                    df_new.loc[j, 'Inh'] = df_new.loc[j, 'Inh'] + ';' + splitted_further_inh[1]
                num_inh_thisentry = num_inh_thisentry + 1

        elif splitted[0] == 'AC':
            splitted_further_act = re.split(r'# | <| [(]', splitted[1])
            prs_now = re.split(',', re.sub('#', '', splitted_further_act[0]))
            IsThisACisOfthisspecies = False
            for k in range(0, len(prs_now)):
                for kk in range(0, len(id4species_here)):
                    if prs_now[k] == id4species_here[kk]:
                        IsThisACisOfthisspecies = True
            if IsThisACisOfthisspecies:
                if num_act_thisentry == 0:
                    df_new.loc[j, 'Act'] = splitted_further_act[1]
                else:
                    df_new.loc[j, 'Act'] = df_new.loc[j, 'Act'] + ';' + splitted_further_act[1]
                num_act_thisentry = num_act_thisentry + 1

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

