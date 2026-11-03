from pprint import pprint

import pypath.resources.urls as urls
from pypath.inputs import new_stitch


def test_actions(mode_type = "activation"):
    
    test_list = []
    url = urls.urls['stitch']['actions'] % 9606

    test_actions = new_stitch.tables(url, max_lines = 5000)

    for test_action in test_actions:
        if test_action['mode'] == mode_type:
            test_list.append(test_action)

    return test_list

def prop_actions(mode_type = "activation"):

    action_list = []
    actions = new_stitch.actions(max_lines = 5000)

    for action in actions:
        if action.mode == mode_type:
            action_list.append(action)

    return action_list

def run_test():

    mode_list = ["expression", 
                 "binding", 
                 "pred_bind", 
                 "reaction", 
                 "activation", 
                 "inhibition", 
                 "catalysis"]
    
    for mode in mode_list:

        print(f"""
################
mode = "{mode}"
################
""")
    
        test_list = test_actions(mode)
        prop_list = prop_actions(mode)

        print(len(test_list))
        print(len(prop_list))

        if len(test_list) / 2 == len(prop_list):
            print("success")
        else:
            print("fail")

        pprint(test_list[0:10])
        print("###############")
        pprint(prop_list[0:5])