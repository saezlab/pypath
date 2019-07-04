import os
import copy

from pypath import main
from pypath import data_formats
from pypath import annot
from pypath import intercell
from pypath import complex
from pypath import settings


class OmniPath(object):
    
    
    def __init__(
        self,
        output_dir = None,
        network_pickle = None,
        annotation_pickle = None,
        intercell_pickle = None,
        complex_pickle = None,
    ):
        
        pass
    
    
    def main(self):
        
        self.load()
    
    
    def load(self):
        
        self.load_complex()
        self.load_network()
        self.load_annotations()
        self.load_intercell()
    
    
    def load_complex(self):
        
        complex.get_db(
            pickle_file = self.ensure_path_exists(self.complex_pickle)
        )
    
    
    def load_network(self):
        
        self.network = main.PyPath()
        
        if os.path.exists(self.network_pickle):
            
            self.network.init_network(pfile = self.network_pickle)
            
        else:
            
            network_input = copy.deepcopy(data_formats.omnipath)
            network_input.update(data_formats.ligand_receptor)
            network_input.update(data_formats.ptm_misc)
            self.network.load_omnipath(
                omnipath = network_input,
                remove_htp = False,
            )
    
    
    def load_annotations(self):
        
        self.annot = annot.get_db(
            pickle_file = self.ensure_path_exists(self.annotation_pickle)
        )
    
    
    def load_intercell(self):
        
        self.intercell = (
            intercell.get_db(
                pickle_file = self.ensure_path_exists(self.intercell_pickle)
            )
        )
    
    
    @staticmethod
    def ensure_path_exists(path):
        
        return path if os.path.exists(path) else None
