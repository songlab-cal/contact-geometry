import argparse
import glob
import json
from pathlib import Path
from os import PathLike
import pickle
from typing import Tuple, List, Dict, Union
import sys

import numpy as np
import pandas as pd

from schrodinger import structure
from schrodinger.structutils import analyze, rings
from ifg import VDM
from functional_groups import ALPHABET

xyxcoord = Tuple[float, float, float]

VDM_SMARTS = ['[C,c]=O', 'CC(=O)[O-]','CC(=O)[OH]', 'CC(C)[C;!R]', '[C,c]C(=O)[N;H2,H1]', '[C,c][N;H2,H1]', '[C,c][O;H1]', '[CX4;H2][S;H1,H0]', 'CSC', 'N=C(N)N', '[O;H1]c1ccccc1', '[C,c][N+;H3]', 'c1c[nH]c2ccccc12', 'c1c[n;H0,H1]cn1', 'c1ccccc1', '[C,c][C;!R]([C;!R])[O;H1]', '[C,c][C;H2;!R][O;H1]', 'c[n;H0,H1]c', '[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)','[CHX4]([CH3X4])[CH3X4]', '[CHX4]([CH3X4])[CH2X4][CH3X4]']
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

def _all_ifgs_restrained(ligand):

    VDM_SMARTS_2 = ['[C,c]=O', '[C,c][CX3](=[OX1])[OH0-,OH]', 'CC(C)[C;!R]', '[C,c][CX3](=[OX1])[NX3H2]','[C,c][N;H2,H1]', '[C,c][O;H1]', '[CX4;H2][S;H1,H0]', 'CSC', 'N=C(N)N', '[O;H1]c1ccccc1', '[C,c][N+;H3]', 'c1c[nH]c2ccccc12', 'c1c[n;H0,H1]cn1', 'c1ccccc1', '[C,c][C;!R]([C;!R])[O;H1]', '[C,c][C;H2;!R][O;H1]', 'c[n;H0,H1]c', '[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)','[CHX4]([CH3X4])[CH3X4]', '[CHX4]([CH3X4])[CH2X4][CH3X4]'] #'[NX3][CX3](=[OX1])[#6]'
    VDM_NAMES_2 = ['carbonyl', 'carboxylate', 'isopropyl', 'carboxamide', 'NH', 'OH', 'thiol', 'thioether', 'guanidine', 'phenol', 'posN', 'indole', 'imidazole', 'phenyl', 'isopropanol', 'ethanol', 'secondary_amine', 'pro_ring', 'val_side', 'ile_side'] #'amide'

    ring_vdms = ['phenol', 'phenyl', 'indole', 'imidazole']
    substruc_matches = list()
    substruc_names = list()

    #lists of ring systems. Two rings are merged into a ring system if they share a bond in common
    ring_sys = rings.find_ring_systems(ligand)

    for i, vdm in enumerate(VDM_SMARTS_2):
        name = VDM_NAMES_2[i]
        substructs = analyze.evaluate_smarts_canvas(ligand, smarts = vdm, start_index = 1, uniqueFilter=True)

        if substructs:
            for match in substructs:
                #avoid including ring vdms if they are part of a ring system
                if name in ring_vdms and any([set(match) < set(r) for r in ring_sys]):
                    continue
                
                substruc_matches.append([ligand.atom[a].pdbname.strip() for a in match])
                substruc_names.append(name)
    
    return substruc_matches, substruc_names

def _all_ifgs_restrained_new_amide(ligand):
    VDM_SMARTS_2 = ['[C,c]=O', '[C,c][CX3](=[OX1])[OH0-,OH]', 'CC(C)[C;!R]', '[C,c][CX3](=[OX1])[NX3H2]','[NX3][CX3](=[OX1])[#6]', '[C,c][O;H1]', '[CX4;H2][S;H1,H0]', 'CSC', 'N=C(N)N', '[O;H1]c1ccccc1', '[C,c][N+;H3]', 'c1c[nH]c2ccccc12', 'c1c[n;H0,H1]cn1', 'c1ccccc1', '[C,c][C;!R]([C;!R])[O;H1]', '[C,c][C;H2;!R][O;H1]', 'c[n;H0,H1]c', '[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)','[CHX4]([CH3X4])[CH3X4]', '[CHX4]([CH3X4])[CH2X4][CH3X4]'] #'[NX3][CX3](=[OX1])[#6]'
    VDM_NAMES_2 = ['carbonyl', 'carboxylate', 'isopropyl', 'carboxamide', 'NH', 'OH', 'thiol', 'thioether', 'guanidine', 'phenol', 'posN', 'indole', 'imidazole', 'phenyl', 'isopropanol', 'ethanol', 'secondary_amine', 'pro_ring', 'val_side', 'ile_side'] #'amide'

    ring_vdms = ['phenol', 'phenyl', 'indole', 'imidazole']
    substruc_matches = list()
    substruc_names = list()

    #if amide, atoms 1 and 0
    #lists of ring systems. Two rings are merged into a ring system if they share a bond in common
    ring_sys = rings.find_ring_systems(ligand)

    for i, vdm in enumerate(VDM_SMARTS_2):
        name = VDM_NAMES_2[i]
        substructs = analyze.evaluate_smarts_canvas(ligand, smarts = vdm, start_index = 1, uniqueFilter=True)

        if substructs:
            for match in substructs:
                #avoid including ring vdms if they are part of a ring system
                if name in ring_vdms and any([set(match) < set(r) for r in ring_sys]):
                    continue
                
                if name == 'NH':
                    print("Found an NH")
                    substruc_matches.append([ligand.atom[a].pdbname.strip() for a in [match[1], match[0]]])
                else:
                    substruc_matches.append([ligand.atom[a].pdbname.strip() for a in match])
                substruc_names.append(name)
    
    return substruc_matches, substruc_names


def _unique_ifgs(all_matches: List[list], all_names: List[str]) -> Tuple[list, list]:
    '''
    Returns unique vdMs from a list of all matches. Any vdM that is fully contained in another vdM
    in the list is removed.

    Args:
    all_matches - all vdms in ligand
    all_names - all names

    Returns:
    final_matches - all unique vdM instances
    final_names - corresponding names
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

def parse_structure_file(filename: Union[str, PathLike], restrictive: bool = False) -> Tuple[dict, dict]:

    infile = Path(filename)
    res_info_file = infile.with_suffix('.json')

    with structure.StructureReader(filename) as st:
        cplx = next(st)

    a_info = [(a.index, a.chain, a.pdbres, a.resnum, a.pdbname, a.xyz, a.inscode) for a in cplx.atom]

    atomdict = dict()

    for at in a_info:
        if at[2].strip() in ALPHABET:
            pdbres = ALPHABET[at[2].strip()]
        else:
            pdbres = at[2].strip()
        idx = at[0] 
        if at[6] != ' ': 
            resid = f"{at[1]}:{at[3]}{at[6]}({pdbres})" 
        else: 
            resid = f"{at[1]}:{at[3]}({pdbres})" 
        at_name = at[4].strip() 
        coords = at[5] 
        if resid not in atomdict: 
            atomdict[resid] = {at_name: {'index':idx, 'coords': coords}} 
        else: 
            atomdict[resid].update({at_name: {'index': idx, 'coords': coords}})

    with open(res_info_file, 'w+') as fout:
        json.dump(atomdict, fout)

    ligchain = cplx.chain['L']
    ligst = ligchain.extractStructure()

    if not restrictive:
        m,n = _all_ifgs(ligst)
    else:
        m, n = _all_ifgs_restrained_new_amide(ligst)
    vdm = _unique_ifgs(m,n)

    #gets mapping from ligand hydrogen to heavy atom for better parsing probe
    #each list only has one entry, so take 0 element
    protonmap = {atom.pdbname.strip():[a.pdbname.strip() for a in atom.bonded_atoms][0] for atom in ligst.atom if atom.element == 'H'}

    return atomdict, vdm, protonmap

def parse_poseviewer_file(poseviewer_file: Union[str, PathLike]) -> dict:

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
            
            rec_res = ligdf['RecResidue'].iloc[i].strip()
            resid = rec_res[rec_res.find('(') +1 : rec_res.find(')')].strip()
            other_info = rec_res[0:rec_res.find('(')]

            if resid in ALPHABET:
                rr = f"{other_info}({ALPHABET[resid]})"
            else:
                rr = rec_res

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

def parse_probe_line(line: str) -> tuple:
    '''
    Parses a probe line to extract salient information about 'hb' or 'cc' contacts to the ligand
    (chain L)
    
    Args:
    line (str) - line of output from probe

    Returns:
    tuple of (interaction type, chain1, resnum1, resname1, atom name1, chain2, resnum2, resname2, atom name2) where 1 and 2
    are the interacting pair
    '''

    spl = line.split(':')[1:]
    INTERACTION = spl[1]
    CHAIN1 = spl[2][:2].strip()
    RESNUM1 = int(spl[2][2:6])
    RESNAME1 = spl[2][6:10].strip()
    NAME1 = spl[2][10:15].strip()
    CHAIN2 = spl[3][:2].strip()
    RESNUM2 = int(spl[3][2:6])
    RESNAME2 = spl[3][6:10].strip()
    NAME2 = spl[3][10:15].strip()
    #deal with waters
    if RESNAME1 == 'HOH':
        NAME1 = 'O'
    if RESNAME2 == 'HOH':
        NAME2 = 'O'

    return (INTERACTION, CHAIN1, RESNUM1, RESNAME1, NAME1, CHAIN2, RESNUM2, RESNAME2, NAME2)

def parse_probe_file(probe_file: Union[str, Path], Hmap: Dict[str, str]) -> dict:
    '''
    Parse the output from probe into a dictionary of interactions keyed by the ligand atoms

    Args:
    probe_file - the probe file name
    Hmap - mapping between ligand Hydrogens and their heavy atoms
    Returns:
    dictionary of interactions between the ligand and protein.
    '''
    probe_dict = {'Ligand_1' : {}}

    with open(probe_file, 'r+') as fin:
        for line in fin:
            interact, ch1, resnum1, resname1, name1, ch2, resnum2, resname2, name2 = parse_probe_line(line)
            if ch1 != 'L' or (ch1 == 'L' and ch2 == 'L'):
                continue
            else:
                if name1.startswith('H') and name1 in Hmap:
                    #proton, assign it to its bonded heavy atom
                    l_atom = Hmap[name1]
                else:
                    l_atom = name1

                rr = f"{ch2}:{resnum2}({resname2})"
                ratom = name2

                if interact == 'hb' or interact == 'cc':

                    int_type = 'hbond' if interact == 'hb' else 'vdw'

                    if l_atom not in probe_dict["Ligand_1"]:
                        probe_dict["Ligand_1"][l_atom] = dict()

                    if rr not in probe_dict["Ligand_1"][l_atom]:
                        probe_dict["Ligand_1"][l_atom][rr] = list()
                    
                    probe_dict["Ligand_1"][l_atom][rr].append((ratom, int_type))

    return probe_dict

def gen_vdms(vdm_dict: Dict[str, str], atom_dict: dict) -> dict:

    vdm_dictionary = dict()

    for entry in atom_dict:
        if entry.startswith('L:'): #ligand is chain L, entries in the atom dict are CH:RESID(NAME)
            vdm_dictionary['Ligand_1'] = []
            lig_info = atom_dict[entry]

            for ifg in vdm_dict:

                #atoms = [lig_info[a] for a in vdm_dict[ifg]]
                
                atoms = vdm_dict[ifg]
                coords = [lig_info[a]['coords'] for a in atoms]
                vdm = VDM(ifg, atoms, coords)

                vdm_dictionary['Ligand_1'].append(vdm)
    
    return vdm_dictionary

def get_residue_info(res, atom_dict):

    atoms = []
    coords = []

    for atm in atom_dict[res]:
        atoms.append(atm)
        coords.append(atom_dict[res][atm]['coords'])
    
    return atoms, coords

def parse_contact(contact, atom_dict):

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

    for lig in vdm_dict:

        lig_interactions = interaction_dict[lig]

        for vdm in vdm_dict[lig]:

            contacts = [(a, lig_interactions[a]) for a in vdm.atoms if a in lig_interactions]

            r = parse_contact(contacts, atom_dict)

            vdm.update_contacts(r)

def generate_all_vdms(base_dir: str, restrictive: bool = False, read_probe: bool =False):

    all_dirs = [d for d in glob.glob(base_dir + '*') if Path(d).is_dir()]
    num_dirs = len(all_dirs)

    for i,dir in enumerate(all_dirs):

        #enumerate directory
        d = Path(dir)
        pdb_id = d.parts[-1]

        if pdb_id.isupper():
            pref = pdb_id +'_proc_out'
        else:
            pref = pdb_id + '_protein_ligand_out'

        maefile = d / (pref + '.mae')

        if not maefile.exists():
            continue

        sys.stdout.write(f"\r{((i+1)/num_dirs) * 100:.1f} % ({i+1}/{num_dirs}) | {pdb_id}")
        sys.stdout.flush()

        #print(f"Working on {pdb_id}")
        pv_file = d / (pdb_id +'_out_pv_interactions.csv')
        probe_file = maefile.with_suffix('.probe')

        atomd, vdms, protonmap= parse_structure_file(maefile, restrictive)

        if not read_probe:
            interactiond = parse_poseviewer_file(pv_file)
        else:
            interactiond = parse_probe_file(probe_file, protonmap)

        if not interactiond:
            print(f"\nError with {pdb_id}")
            continue

        vdmd = gen_vdms(vdms, atomd)

        populate_vdms(vdmd, atomd, interactiond)

        if not restrictive:
            vdm_file = d / (pdb_id +'_vdms.pkl')
        else:
            vdm_file = d / (pdb_id + '_restrictive_vdms.pkl')

        if read_probe:
            outfile = vdm_file.parent / (vdm_file.stem + '_probe_extracted_alternate_NH.pkl')
        else:
            outfile = vdm_file
        
        with open(outfile, 'wb') as fout:
            pickle.dump(vdmd, fout)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process Maestro Files To Extract VDMS")

    parser.add_argument("-i", "--input_dir", type=str, help="Input Directory")
    parser.add_argument("-restrictive", default = False, action= 'store_true', help= "Make iFG definitions more restrictive (strictly like they appear on side chains")
    parser.add_argument("-probe", default=False, action= 'store_true', help = "Parse interactions from probe")

    args = parser.parse_args()

    if args.restrictive:
        strict = True
    else:
        strict = False

    if args.probe:
        use_probe = True
    else:
        use_probe = False

    generate_all_vdms(args.input_dir, strict, use_probe)
    print(f"\n")