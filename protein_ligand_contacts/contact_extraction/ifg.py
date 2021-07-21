from functional_groups import SMARTS_ATOM_MAP
from typing import List, Tuple

import numpy as np

#typedef
xyzcoord = Tuple[float, float, float]

class iFG():
    '''
    iFG class to store information about ligand atoms that make up an iFG, mapping to amino acid residues,
    as well as the protein residues that interact with the iFG
    '''
    def __init__(self, name:str, cfg_entry:dict):
        
        self.canonical_name = name.split('-')[0]
        self.name = name
        
        atom_list = []
        coord_list = []
        
        for atm in cfg_entry:
            atom_list.append(atm)
            coord_list.append(cfg_entry[atm])
        
        self.lig_atoms = atom_list
        self.heavy_atom_coords = np.vstack(coord_list)
        
        self.aa_atoms = None
        
        self.contacts = dict()

        self.ala_ref = np.array([[ 0.11800566, -2.45811057,  0. ], \
                                 [-0.77976233, -1.30825949,  0. ], \
                                 [ 0.        ,  0.        ,  0. ]])

    def __repr__(self):

        return self.name

    def get_aa_atom_names(self, SMARTS_ATOM_MAP):
        
        lig_aa_mapping = {}
        
        for restype in SMARTS_ATOM_MAP[self.canonical_name]:
            
            atoms = SMARTS_ATOM_MAP[self.canonical_name][restype]
            
            lig_aa_mapping[restype] = dict(zip(self.lig_atoms, atoms))

        self.aa_atoms = lig_aa_mapping

    def add_interacting_residue(self, interaction_type:str, resid:str, name2:str, name1:str, idx_coord_map:dict):
        '''
        if interaction_type == 'cc':
        
            if resid not in self.ccs:
                self.ccs.update({resid: {'atoms': idx_coord_map[resid]['names'], 'coords': idx_coord_map[resid]['coords'], 'interaction': [(name2 ,name1)]}})
            else:
                self.ccs[resid]['interaction'].append((name2, name1))
        
        elif interaction_type == 'hb':
            
            if resid not in self.hbs:
                self.hbs.update({resid: {'atoms': idx_coord_map[resid]['names'], 'coords': idx_coord_map[resid]['coords'], 'interaction': [(name2 ,name1)]}})
            else:
                self.hbs[resid]['interaction'].append((name2, name1))
        '''
        if resid not in self.contacts:
            
            self.contacts.update({resid: {'residue_atoms': idx_coord_map[resid]['names'], \
                                          'residue_coords': np.array(idx_coord_map[resid]['coords']),
                                          'interactions': {'hb': set(), 'cc': set()}}})

        self.contacts[resid]['interactions'][interaction_type].add((name2,name1))

class VDM():

    def __init__(
        self,
        name: str,
        atom_names: List[str],
        atom_coords: List[xyzcoord],
    ):
        self.name = name
        self.canonical_name = name.split('_')[0]
        self.atoms = atom_names
        self.coords = np.array(atom_coords)

        self.aa_atoms = None

        self.contacts = dict()

    
    def get_aa_atoms(self):
        return self.aa_atoms

    def set_aa_atoms(self):
        global SMARTS_ATOM_MAP
        lig_aa_mapping = dict()
        for restype in SMARTS_ATOM_MAP[self.canonical_name]:
            
            atoms = SMARTS_ATOM_MAP[self.canonical_name][restype]
            
            lig_aa_mapping[restype] = dict(zip(self.atoms, atoms))

        self.aa_atoms = lig_aa_mapping

    def __repr__(self):
        return self.name

    def update_contacts(self, new_contact: dict):
        self.contacts.update(new_contact)