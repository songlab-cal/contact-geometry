import numpy as np
import pickle as pkl

from functional_groups import SMARTS_ATOM_MAP

DATA_FOLDER         = 'raw_contacts'
VDMS_FOLDER         = 'contacts'
IFG_REPS_FOLDER     = 'ifg_reps'

AMINO_ACIDS = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU',
        'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT', 'H', 'H2', 'HA2', 'HXT'}
SYMMETRIC_IFGS = {'phenol': [[2,6], [3,5]], 'val_side': [[0,2]], 'guanidine': [[0,3]], \
        'secondary_amine': [[0,2]], 'isopropyl': [[0,2,3]], 'phenyl': [[0,1,2,3,4,5]], \
        'carboxylate': [[2,3]], 'thioether': [[0,2]], 'pro_ring': [[1,4], [2,3]]}

VDM_TYPES = [(aa, ifg, h_bonding, sidechain_bonding) for aa in AMINO_ACIDS for ifg in \
        SMARTS_ATOM_MAP.keys() for h_bonding in (True, False) for sidechain_bonding in \
        (True, False)]

RMSD_CUTOFF = 0.5
MAX_NUM_CONTACTS = float('inf')
MIN_DATA_POINTS = 30
ALPHA = 0.001

with open('data/{0}/ideal_ala_CO.pkl'.format(DATA_FOLDER), 'rb') as pkl_file:
    ala = pkl.load(pkl_file)

IDEAL_ALA_REF = np.array(ala.loc[[0,1,2], ['c_x', 'c_y', 'c_z']])

def is_H(atom):
    return (atom[1] if atom[0].isdigit() else atom[0])=='H'

def get_label(vdm_type):
    aa, ifg, h_bonding, sidechain_bonding = vdm_type
    return '{0}-{1}-h_bonding_{2}-sidechain_bonding_{3}'.format(aa, ifg, \
                h_bonding, sidechain_bonding)

def kabsch_algorithm(mob, targ):
    
    mob_com = np.mean(mob, 0)
    targ_com = np.mean(targ, 0)
    
    mob_cent = mob - mob_com
    targ_cent = targ - targ_com
    
    cov = mob_cent.T @ targ_cent
    
    V, S, W = np.linalg.svd(cov)
    
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    
    if d:
        V[:,-1] *= -1
        
    R = V @ W
    
    return R, mob_com, targ_com

# Align N, CA, C atoms to those of IDEAL_ALA_REF
def get_ref_alignment(backbone):
    R, m, t = kabsch_algorithm(backbone, IDEAL_ALA_REF)
    return lambda pts : ((pts - m) @ R) + t
