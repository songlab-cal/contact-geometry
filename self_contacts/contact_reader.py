import pickle as pkl
import numpy as np
import math

from functional_groups import SMARTS_ATOM_MAP
from utils import AMINO_ACIDS, BACKBONE_ATOMS, DATA_FOLDER, VDMS_FOLDER, \
        IFG_REPS_FOLDER, VDM_TYPES, get_label

def compute_vdms():

    with open('data/{0}/aa_seqs.pkl'.format(DATA_FOLDER), 'rb') as pkl_file:
        aa_seqs = pkl.load(pkl_file)

    with open('data/{0}/ifg_Hs.pkl'.format(DATA_FOLDER), 'rb') as pkl_file:
        ifg_Hs = pkl.load(pkl_file)

    with open('data/{0}/all_contacts.pkl'.format(DATA_FOLDER), 'rb') as pkl_file:
        contacts = pkl.load(pkl_file)

    vdms = {vdm_type: [[], []] for vdm_type in VDM_TYPES}

    contained_ifgs = {aa: set() for aa in AMINO_ACIDS}
    for ifg in SMARTS_ATOM_MAP.keys():
        for aa in SMARTS_ATOM_MAP[ifg].keys():
            contained_ifgs[aa].add(ifg)
    
    remaining_ifgs = set(SMARTS_ATOM_MAP.keys())

    def update_vdms(aa1_data, aa2_data, aa1_fxl_atoms, aa2_fxl_atoms, bonds):

        aa1            = aa1_data[0]
        aa2, aa2_atoms = aa2_data[0], aa2_data[1]

        for ifg in contained_ifgs[aa2]:

            if not any([atom in aa2_fxl_atoms for atom in SMARTS_ATOM_MAP[ifg][aa2]]):
                continue

            ifg_i = [aa2_atoms.index(atom) for atom in SMARTS_ATOM_MAP[ifg][aa2]]
            ifg_pts = [aa2_data[2][i] for i in ifg_i]

            h_bonding = any([atom_pair[1] in SMARTS_ATOM_MAP[ifg][aa2] or \
                    atom_pair[1] in ifg_Hs[(ifg, aa2)] for atom_pair in bonds['hb']])

            sidechain_bonding = not all([atom in BACKBONE_ATOMS for atom in aa1_fxl_atoms])

            label = (aa1, ifg, h_bonding, sidechain_bonding)

            # Separate backbone/ifg and sidechain atoms to perform downstream 
            # rmsd computations on backbone/ifg only
            vdms[label][0].append(np.array(list(aa1_data[2][:3])+ifg_pts).astype('float32'))
            vdms[label][1].append(aa1_data[2][3:len(aa_seqs[aa1])].astype('float32'))

            if not ifg in remaining_ifgs:
                continue

            with open('data/{0}/{1}.pkl'.format(IFG_REPS_FOLDER, ifg), 'wb') as pkl_file:
                pkl.dump(ifg_pts, pkl_file)

            remaining_ifgs.remove(ifg)

    for aa_pair in contacts.keys():

        aa1, aa2 = aa_pair

        if aa1 not in AMINO_ACIDS or aa2 not in AMINO_ACIDS:
            continue

        for contact in contacts[aa_pair]:

            contact[0][1][:] = [x.strip() for x in contact[0][1]]
            contact[1][1][:] = [x.strip() for x in contact[1][1]]

            # Contact atoms should either match AA_SEQS[AA_NAME] exactly up to the
            # second-to-last non-H atom. We ignore the following atoms, which are
            # likely 'OXT', followed by any hydrogens.
            if len(contact[0][1]) < len(aa_seqs[aa1]) or \
                    len(contact[1][1]) < len(aa_seqs[aa2]) or \
                    not all([contact[0][1][i]==aa_seqs[aa1][i]
                        for i in range(len(aa_seqs[aa1]))]) or \
                    not all([contact[1][1][i]==aa_seqs[aa2][i]
                        for i in range(len(aa_seqs[aa2]))]):
                continue
            
            aa1_fxl_atoms, aa2_fxl_atoms = set(), set()

            for typed_contacts in contact[2].values():
                for pair in typed_contacts:
                    aa1_fxl_atoms.add(pair[0])
                    aa2_fxl_atoms.add(pair[1])
            
            update_vdms(contact[0], contact[1], aa1_fxl_atoms, aa2_fxl_atoms, contact[2])
            update_vdms(contact[1], contact[0], aa2_fxl_atoms, aa1_fxl_atoms, \
                    {bond_type: {(x[1], x[0]) for x in contact[2][bond_type]} for \
                    bond_type in contact[2].keys()})

    for vdm_type in VDM_TYPES:

        vdms[vdm_type][0] = np.array(vdms[vdm_type][0])
        vdms[vdm_type][1] = np.array(vdms[vdm_type][1])

        label = get_label(vdm_type)

        folder = VDMS_FOLDER

        with open('data/{0}/{1}.pkl'.format(folder, label), 'wb') as pkl_file:
            pkl.dump(vdms[vdm_type], pkl_file)
    
    return vdms

if __name__ == "__main__":
    compute_vdms()
