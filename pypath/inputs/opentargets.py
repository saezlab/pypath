from pypath.share import curl
from pypath.resources.urls import urls

import json
import re

def opentargets_general(base_url, primary_k):
    file_url = base_url + "/%s"

    c = curl.Curl(
        base_url,
        silent=False,
        large=False,
        encoding="utf-8",
        default_mode="r",
    )

    html = c.result

    regex = r'"(part.*\.json)"'
    matcher = re.compile(regex)
    filenames = matcher.findall(html)

    result = {}

    for filename in filenames:

        c = curl.Curl(
            file_url % filename,
            silent=False,
            large=False,
            encoding="utf-8",
            default_mode="r",
        )

        lines = c.result.split("\n")

        for line in lines:

            if line == "":
                continue

            try:

                tmp_dict = json.loads(line)

            except json.JSONDecodeError:

                print(
                    f'Following line could not be read properly:\n{line}'
                )

                continue

            key = tmp_dict[primary_k]
            del tmp_dict[primary_k]

            try:
                result[key].append(tmp_dict)
            except KeyError:
                result[key] = [tmp_dict]
    
    return result


def overall_indirect_score():
    base_url = urls['opentargets']['url'] % 'associationByOverallIndirect'
    return opentargets_general(base_url, 'diseaseId')


def overall_direct_score():
    base_url = urls['opentargets']['url'] % 'associationByOverallDirect'
    return opentargets_general(base_url, 'diseaseId')


def significant_adverse_drug_reactions():
    base_url = urls['opentargets']['url'] % 'fda/significantAdverseDrugReactions'
    return opentargets_general(base_url, 'chembl_id')


def baseline_expression():
    base_url = urls['opentargets']['url'] % 'baselineExpression'
    return opentargets_general(base_url, 'id')