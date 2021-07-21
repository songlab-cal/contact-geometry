import argparse
import glob
import json
from pathlib import Path
from os import PathLike
import pickle
from typing import Tuple, List, Dict, Union

import numpy as np
import pandas as pd

from schrodinger import structure
from schrodinger.structutils import analyze
from ifg import VDM

xyxcoord = Tuple[float, float, float]

VDM_SMARTS = ['[C,c]=O', 'CC(=O)[O-]','CC(=O)[OH]', 'CC(C)[C;!R]', '[C,c]C([N;H2,H1])=O', '[C,c][N;H2,H1]', '[C,c][O;H1]', '[CX4;H2][S;H1,H0]', 'CSC', 'N=C(N)N', '[O;H1]c1ccccc1', '[C,c][N+;H3]', 'c1c[nH]c2ccccc12', 'c1c[n;H0,H1]cn1', 'c1ccccc1', '[C,c][C;!R]([C;!R])[O;H1]', '[C,c][C;H2;!R][O;H1]', 'c[n;H0,H1]c', '[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)','[CHX4]([CH3X4])[CH3X4]', '[CHX4]([CH3X4])[CH2X4][CH3X4]']
VDM_NAMES = ['carbonyl', 'carboxylate', 'carboxylate', 'isopropyl', 'carboxamide', 'NH', 'OH', 'thiol', 'thioether', 'guanidine', 'phenol', 'posN', 'indole', 'imidazole', 'phenyl', 'isopropanol', 'ethanol', 'secondary_amine', 'pro_ring', 'val_side', 'ile_side']

def _all_ifgs(ligand) -> Tuple[list, list]:
    '''
    Gets all ifgs and their names for a particular ligand

    Args:
    ligand - ligand structure object (only one needed)

    Returns:
    list of ifg atoms, list of ifg names
    '''

    substruc_matches = list()
    substruc_names = list()

    for i, vdm in enumerate(VDM_SMARTS):

        substructs = analyze.evaluate_smarts_canvas(ligand, smarts = vdm, start_index = 1, uniqueFilter=True)

        if len(substructs) == 1 and substructs[0]:
            substruc_matches.append([ligand.atom[a].pdbname.strip() for a in substructs[0]])
            substruc_names.append(VDM_NAMES[i])
        else:
            for match in substructs:
                substruc_matches.append([ligand.atom[a].pdbname.strip() for a in match])
                substruc_names.append(VDM_NAMES[i])
    
    return substruc_matches, substruc_names

def _unique_ifgs(all_matches: List[list], all_names: List[str]) -> Dict[str, list]:
    '''
    Returns unique iFGs from a list of all matches. Any iFG that is fully contained in another vdM
    in the list is removed.

    Args:
    all_matches - all iFGs in ligand
    all_names - all atom names in ligand

    Returns:
    vdm_dict - dictionary indexed by iFG names returning the corresponding atoms that make up that iFG
    '''

    final_matches = []
    final_names = []

    for i,q in enumerate(all_matches):
        
        #check whether iFG is contained in any other vdM iFG
        if not sum([set(q).issubset(set(t)) for j,t in enumerate(all_matches) if j != i]):
            final_matches.append(q)
            final_names.append(all_names[i])
        

    vdm_dict = dict()

    for i in range(len(final_matches)):
        vdm_dict[f"{final_names[i]}_{i}"] = final_matches[i]

    return vdm_dict

def parse_structure_file(filename: Union[str, PathLike]) -> Tuple[dict, dict]:

    '''
    Parses a Maestro file to get information about atoms in the structure, and extracts the iFGs from the ligand.

    Args:
    filename - input filename (.mae). Ideally this file has been preprocessed by the protein preparation wizard

    Returns:
    atomdict - Dictionary of all atoms and their coordinates in the structure, indexed by residue in the following format:
                CH:RESNUM(RESNAME)
    vdm - dictionary of each unique iFG and the ligand atoms that make it up
    '''

    infile = Path(filename)
    res_info_file = infile.with_suffix('.json')

    with structure.StructureReader(filename) as st:
        cplx = next(st)

    a_info = [(a.index, a.chain, a.pdbres, a.resnum, a.pdbname, a.xyz, a.inscode) for a in cplx.atom]

    atomdict = dict()

    for at in a_info:
        idx = at[0] 
        if at[6] != ' ': 
            resid = f"{at[1]}:{at[3]}{at[6]}({at[2].strip()})" 
        else: 
            resid = f"{at[1]}:{at[3]}({at[2].strip()})" 
        at_name = at[4].strip() 
        coords = at[5] 
        if resid not in atomdict: 
            atomdict[resid] = {at_name: {'index':idx, 'coords': coords}} 
        else: 
            atomdict[resid].update({at_name: {'index': idx, 'coords': coords}})
    #save processed data
    with open(res_info_file, 'w+') as fout:
        json.dump(atomdict, fout)
    #ligand chain object will be Chain L following our preprocessing
    ligchain = cplx.chain['L']
    ligst = ligchain.extractStructure()

    m,n = _all_ifgs(ligst)
    vdm = _unique_ifgs(m,n)

    return atomdict, vdm

def parse_poseviewer_file(poseviewer_file: Union[str, PathLike]) -> Dict[str, List[Tuple[str, str]]]:
    '''
    Parses the output of SCRODINGER's poseviwer_interactions.py

    Args:
    poseviewer_file - csv file output from poseviewer_interactions.py, requires the "-csv" flag when running

    Returns:
    Dictionary of interactions broken down by ligand. Each entry is a ligand atom, and the values are lists of tuples of 
    receptor atom, and type of interaction as identified by the poseviewer_interactions script
    '''

    lig_interactions = dict()

    df = pd.read_csv(poseviewer_file, sep=',')
    for lig in np.array(df.Ligand.unique()):

        ligdf = df[df.Ligand == lig]

        lig_interactions[f"Ligand_{lig}"] = {}

        for i in range(len(ligdf)):
            la = ligdf['LigAtom'].iloc[i]
            latom = la[la.find('(')+1: la.find(')')].strip()

            ra = ligdf['RecAtom'].iloc[i]
            ratom = ra[ra.find('(')+1 : ra.find(')')].strip()

            rr = ligdf['RecResidue'].iloc[i].strip()

            interaction_type = ligdf.Type.iloc[i]
            if interaction_type.startswith('HDonor') or interaction_type.startswith('HAccep') or interaction_type == 'Ar-Hbond':
                int_type = 'hbond'
            elif interaction_type != 'Ugly' or interaction_type != 'Bad':
                int_type = 'vdw'
            else:
                #ignore 'Ugly' and 'Bad' contacts. These were also ignored in protein-protein dataset, so should be kosher
                continue

            if latom not in lig_interactions[f"Ligand_{lig}"]:
                lig_interactions[f"Ligand_{lig}"][latom] = dict()
            
            if rr not in lig_interactions[f"Ligand_{lig}"][latom]:
                lig_interactions[f"Ligand_{lig}"][latom][rr] = list()
            lig_interactions[f"Ligand_{lig}"][latom][rr].append((ratom, int_type))
    
    return lig_interactions

def gen_vdms(vdm_dict: Dict[str, str], atom_dict: dict) -> dict:
    '''
    Initializes vdM objects from the iFGs and atomic coordinates extracted by parse_structure_file

    Args:
    vdm_dict - dictionary of iFGs and their corresponding atoms
    atom_dict - dictionary of atomic coordinates and indices, indexed first by residue

    Returns:
    Dictionary with a single Ligand (for this demonstration) that yields a list of iFGs
    '''

    vdm_dictionary = dict()

    for entry in atom_dict:
        if entry.startswith('L:'): #ligand is chain L, entries in the atom dict are CH:RESID(NAME)
            vdm_dictionary['Ligand_1'] = []
            lig_info = atom_dict[entry]

            for ifg in vdm_dict:

                #atoms = [lig_info[a] for a in vdm_dict[ifg]]
                atoms = vdm_dict[ifg]
                coords = [lig_info[a]['coords'] for a in vdm_dict[ifg]]
                #instantiate the VDM object with name and coordinates
                vdm = VDM(ifg, atoms, coords)

                vdm_dictionary['Ligand_1'].append(vdm)
    
    return vdm_dictionary

def get_residue_info(res:str, atom_dict:dict) -> Tuple[list, list]:
    '''
    Returns atoms and coordinates for a given residue as two lists

    Args:
    res - residue name in format CH:RESNUM(RESNAME) as established in parse_structure_file
    atom_dict - dictionary of atomic coordinates and indices, indexed first by residue

    Returns:
    list of atoms, list of coordinates for that residue
    '''

    atoms = []
    coords = []

    for atm in atom_dict[res]:
        atoms.append(atm)
        coords.append(atom_dict[res][atm]['coords'])
    
    return atoms, coords

def parse_contact(contact: List[Tuple[str, str]], atom_dict: Dict[str, dict]) -> dict:
    '''
    Parses a set of contacts for a given ligand atom

    Args:
    contact - list of contacts to a given atom (of the form (atom, contact type))
    atom_dict - dictionary of atomic coordinates and indices, indexed first by residue

    Returns:
    Dictionary of contacts for the ligand atom
    '''

    contact_dict = dict()

    for ct in contact:
        la = ct[0]
        ctx = ct[1]

        for entry in ctx:
            if entry not in atom_dict:
                continue
            if entry not in contact_dict:
                atoms, coords = get_residue_info(entry, atom_dict)
                contact_dict[entry] = {'residue_atoms': atoms, 'residue_coords': coords, 'interactions': {'hb': set(), 'cc': set()}}
            
            hbs = set((la, item[0]) for item in ctx[entry] if item[1] == 'hbond')
            ccs = set((la, item[0]) for item in ctx[entry] if item[1] == 'vdw')

            ccs.difference_update(hbs)

            contact_dict[entry]['interactions']['hb'].update(hbs)
            contact_dict[entry]['interactions']['cc'].update(ccs)

    return contact_dict

def populate_vdms(vdm_dict, atom_dict, interaction_dict) -> None:

    '''
    Populate contact information for the vdMs identified from the ligand

    Args:
    vdm_dict - dictionary of vdM objects from gen_vdMs
    atom_dict - dictionary of atomic coordinates and indices, indexed first by residue
    interaction_dict - dictionary of interactions indexed by ligand atom

    Returns:
    None - updates the contacts for each vdM
    '''

    for lig in vdm_dict:

        lig_interactions = interaction_dict[lig]

        for vdm in vdm_dict[lig]:

            contacts = [(a, lig_interactions[a]) for a in vdm.atoms if a in lig_interactions]

            r = parse_contact(contacts, atom_dict)

            vdm.update_contacts(r)

def generate_all_vdms(base_dir: str) -> None:

    '''
    Generates vdMs using a pair of preprocessed maestro and poseviewer interaction files
    vdms are saved as pdb_name + '_vdms.pkl). ie if the pdb file was originally 4DKL, the processed maestro file is 4DKL_proc_out.mae
    and the vdm file will be 4DKL_vdms.pkl

    Args:
    Base directory with subdirectories for each PDB file. 

    Returns:
    None - generates and saves vdMs
    '''

    for dir in glob.glob(base_dir + '*'):

        d = Path(dir)
        pdb_id = d.parts[-1]

        print(f"Working on {pdb_id}")

        maefile = d / (pdb_id + '_proc_out.mae')
        pv_file = d / (pdb_id +'_proc_out_pv_interactions.csv')

        atomd, vdms = parse_structure_file(maefile)
        interactiond = parse_poseviewer_file(pv_file)

        vdmd = gen_vdms(vdms, atomd)

        populate_vdms(vdmd, atomd, interactiond)

        vdm_file = d / (pdb_id +'_vdms.pkl')

        with open(vdm_file, 'wb') as fout:
            pickle.dump(vdmd, fout)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process Maestro Files To Extract VDMS")

    parser.add_argument("-i", "--input_dir", type=str, help="Input Directory")

    args = parser.parse_args()

    generate_all_vdms(args.input_dir)