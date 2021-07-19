SMARTS_ATOM_MAP = {
    'carbonyl': {
        'ASN': ['CG', 'OD1'],
        'GLN': ['CD', 'OE1'],
        'ASP': ['CG', 'OD1'],
        'GLU': ['CD', 'OE1'],
        'GLY': ['C', 'O'],
        'ALA': ['C', 'O'],
    },
    'carboxylate': {
        'ASP': ['CB', 'CG', 'OD1', 'OD2'],
        'GLU': ['CG', 'CD', 'OE1', 'OE2'],
    },
    'isopropyl': {
        'LEU': ['CB', 'CG', 'CD1', 'CD2']
    },
    'carboxamide': {
        'ASN': ['CB', 'CG', 'OD1', 'ND2'],
        'GLN': ['CG', 'CD', 'OE1', 'NE2']
    },
    'NH' : {
        'ASN': ['CG', 'ND2'],
        'GLN': ['CD', 'NE2'],
        'GLY': ['CA', 'N'],
        'ALA': ['CA', 'N']
    },
    'OH' : {
        'ASP': ['CG', 'OD2'],
        'GLU': ['CD', 'OE2'],
        'TYR': ['CZ', 'OH'],
        'SER': ['CB', 'OG'],
        'THR': ['CB', 'OG1']
    },
    'thiol' :{
        'CYS': ['CB', 'SG']
    },
    'thioether' : {
        'MET': ['CG', 'SD', 'CE']
    },
    'guanidine': {
        'ARG': ['NH2', 'CZ', 'NE', 'NH1']
    },
    'phenol': {
        # vv canonical, alt: ['OH', 'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2']
        'TYR': ['OH', 'CZ', 'CE2', 'CD2', 'CG', 'CD1', 'CE1'],
    },
    'posN': {
        'LYS': ['CE', 'NZ']
    },
    'indole': {
        'TRP': ['CG', 'CD1', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2']
    },
    'imidazole': {
        'HIS': ['CG', 'CD2', 'NE2', 'CE1', 'ND1']
    },
    'phenyl': {
        'PHE': ['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']
    },
    'isopropanol': {
        'THR': ['CG2', 'CB', 'CA', 'OG1']
    },
    'ethanol': {
        'SER': ['CA', 'CB', 'OG']
    },
    'secondary_amine': {
        'HIS': ['CE1', 'NE2', 'CD2']
    },
    'val_side': {
        'VAL': ['CG1','CB', 'CG2']
    },
    'ile_side': {
        'ILE': ['CG2', 'CB', 'CG1','CD1']
    },
    'pro_ring': {
        'PRO': ['N', 'CD', 'CG', 'CB', 'CA']
    }
}

PROTON_ATOM_MAP = {'ALA': {'HA': 'CA', 'H': 'N', 'HB1': 'CB', 'HB2': 'CB', 'HB3': 'CB'},
 'ARG': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HG2': 'CG',
  'HG3': 'CG',
  'HD2': 'CD',
  'HD3': 'CD',
  'HE': 'NE',
  'HH11': 'NH1',
  'HH12': 'NH1',
  'HH21': 'NH2',
  'HH22': 'NH2'},
 'ASN': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HD21': 'ND2',
  'HD22': 'ND2'},
 'ASP': {'HA': 'CA', 'H': 'N', 'HB2': 'CB', 'HB3': 'CB'},
 'CYS': {'HA': 'CA', 'H': 'N', 'HB2': 'CB', 'HB3': 'CB', 'HG': 'SG'},
 'GLN': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HG2': 'CG',
  'HG3': 'CG',
  'HE21': 'NE2',
  'HE22': 'NE2'},
 'GLU': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HG2': 'CG',
  'HG3': 'CG'},
 'GLY': {'HA2': 'CA', 'HA3': 'CA', 'H': 'N'},
 'HIS': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HD1': 'ND1',
  'HE1': 'CE1',
  'HD2': 'CD2'},
 'ILE': {'HA': 'CA',
  'H': 'N',
  'HB': 'CB',
  'HG22': 'CG2',
  'HG21': 'CG2',
  'HG23': 'CG2',
  'HG12': 'CG1',
  'HG13': 'CG1',
  'HD11': 'CD1',
  'HD12': 'CD1',
  'HD13': 'CD1'},
 'LEU': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HG': 'CG',
  'HD11': 'CD1',
  'HD12': 'CD1',
  'HD13': 'CD1',
  'HD21': 'CD2',
  'HD22': 'CD2',
  'HD23': 'CD2'},
 'LYS': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HG2': 'CG',
  'HG3': 'CG',
  'HD2': 'CD',
  'HD3': 'CD',
  'HE2': 'CE',
  'HE3': 'CE',
  'HZ1': 'NZ',
  'HZ2': 'NZ',
  'HZ3': 'NZ'},
 'MET': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HG2': 'CG',
  'HG3': 'CG',
  'HE1': 'CE',
  'HE2': 'CE',
  'HE3': 'CE'},
 'PHE': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HD1': 'CD1',
  'HD2': 'CD2',
  'HE1': 'CE1',
  'HE2': 'CE2',
  'HZ': 'CZ'},
 'PRO': {'HA': 'CA',
  'HB2': 'CB',
  'HB3': 'CB',
  'HG2': 'CG',
  'HG3': 'CG',
  'HD2': 'CD',
  'HD3': 'CD'},
 'SER': {'HA': 'CA', 'H': 'N', 'HB2': 'CB', 'HB3': 'CB', 'HG': 'OG'},
 'THR': {'HA': 'CA',
  'H': 'N',
  'HB': 'CB',
  'HG1': 'OG1',
  'HG21': 'CG2',
  'HG22': 'CG2',
  'HG23': 'CG2'},
 'TRP': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HD1': 'CD1',
  'HE1': 'NE1',
  'HZ2': 'CZ2',
  'HH2': 'CH2',
  'HZ3': 'CZ3',
  'HE3': 'CE3'},
 'TYR': {'HA': 'CA',
  'H': 'N',
  'HB2': 'CB',
  'HB3': 'CB',
  'HD1': 'CD1',
  'HD2': 'CD2',
  'HE1': 'CE1',
  'HE2': 'CE2',
  'HH': 'OH'},
 'VAL': {'HA': 'CA',
  'H': 'N',
  'HB': 'CB',
  'HG11': 'CG1',
  'HG12': 'CG1',
  'HG13': 'CG1',
  'HG21': 'CG2',
  'HG22': 'CG2',
  'HG23': 'CG2'}}
