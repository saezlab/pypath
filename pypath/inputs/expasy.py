
import re
from pypath.share import curl
from pypath.resources.urls import urls

def enzyme_classes():

    url = urls['expasy']['url'] % 'enzclass.txt'

    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding="utf-8",
        default_mode="r",
    )

    regex = r"^([\d ]+\.[\d\- ]+\.[\d\- ]+\. ?[\d\-]+) +(.+)\.$"
    matcher = re.compile(regex)
    tree = Tree(delimiter='.')

    for line in c.result:

        try:
            result = matcher.findall(line)[0]
        except IndexError:
            continue
        
        keys, value = result
        tree.parse_add(keys, value)
    
    return tree

def enzymes():

    url = urls['expasy']['url'] % 'enzyme.dat'

    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding="utf-8",
        default_mode="r",
    )

    id_regex = r'^ID +([\d\w]+\.[\d\w]+\.[\d\w]+\.[\d\w]+)$'
    id_matcher = re.compile(id_regex)

    de_regex = r'^DE +(.+)\.?$'
    de_matcher = re.compile(de_regex)

    # Not implemented
    an_regex = r''
    an_matcher = re.compile(an_regex)

    dr_regex = r'([^ ,]+), {0,}[^;]+ {0,};'
    dr_matcher = re.compile(dr_regex)

    enzyme = {
        'name': [],
        'class_description': None,
        'annotations': []
    }

    ec_tree = process_classes()

    result = {}

    for line in c.result:
        
        if line.startswith('CC'):
            continue

        elif line.startswith('//'):

            try:
                enzyme['class_description'] = ec_tree.describe(enzyme_id)
                result[enzyme_id] = enzyme

                enzyme = {
                    'name': [],
                    'class_description': None,
                    'annotations': []
                }
            except NameError:
                continue

        elif line.startswith('ID'):
            try:
                enzyme_id = id_matcher.findall(line)[0]
            except IndexError:
                print(line)
                print(id_matcher.findall(line))

        elif line.startswith('DE'):
            enzyme_desc = de_matcher.findall(line)
            enzyme['name'].extend(enzyme_desc)

        elif line.startswith('DR'):
            enzyme_xref = dr_matcher.findall(line)
            enzyme['annotations'].extend(enzyme_xref)
    
    return result



class Tree:

    def __init__(self, delimiter=None):
        self.root = Node('root')
        self.delimiter = delimiter

    def parse_add(self, keys, value):
        keys = keys.split(self.delimiter)
        keys = [int(key) for key in keys if key.strip() != '-']

        tmp_node = self.root

        for key in keys:
            tmp_node.add_child(key)
            tmp_node = tmp_node.get_child(key)
        
        tmp_node.set_data(value)
    
    def describe(self, keys):
        description = []

        keys = keys.split(self.delimiter)

        tmp_node = self.root

        for key in keys:

            if key.strip() == '-':
                break

            try:
                tmp_node = tmp_node.get_child(int(key))
            except ValueError:
                break
        
            try:
                description.append(tmp_node.get_data())
            except AttributeError:
                continue
        
        return ', '.join(description)

    def __str__(self):
        string = str(self.root)
        string = string.strip('\n')

        return string

    def __len__(self):
        pass

class Node:

    def __init__(self, key, data=None):

        self.key = key
        self.data = data
        self.children = {}
    
    def add_child(self, child_key, data=None):
        
        if child_key not in self.children:
            self.children[child_key] = Node(child_key, data)
            return True
        
        return False
    
    def remove_child(self, child_key):

        try:
            del self.children[child_key]
            return True
        except KeyError:
            return False
    
    def get_child(self, child_key):

        try:
            return self.children[child_key]
        except KeyError:
            return None

    def set_data(self, data):
        self.data = data

    def get_data(self):
        return self.data
    
    def __str__(self):
        
        string = ''
        prefix = ''

        if self.key != 'root':
            string += f'{self.key}: {self.data}\n'
            prefix = '\t'

        for child_key in self.children:

            child_string = str(self.children[child_key])
            lines = child_string.split('\n')
            lines = [line for line in lines if line]

            for line in lines:
                string += f'{prefix}{line}\n'

        return string

    def __len__(self):
        return len(self.children)
