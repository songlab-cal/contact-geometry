import argparse
import glob
from os import PathLike
from pathlib import Path
from subprocess import Popen, PIPE
import sys
from typing import List, Union

import numpy as np
import pandas as pd

from schrodinger import structure
from schrodinger.structutils import analyze, measure

def lig_avg_bfact(lig:analyze.Ligand) -> float:
    '''
    Calculates the average b-factor (temperature factor) across a ligand. 

    Args:
    lig - ligand object from schrodinger
    Returns:
    float value of average b-factor across the ligand
    '''
    #deal with protons on some ligands that have Nonetype temperature factors with None check
    return np.average(np.array([a.temperature_factor for a in lig.st.atom if a.temperature_factor is not None]))

def get_chains_around_lig(cplx: structure.Structure, lig: analyze.Ligand, distance: float = 5.0) -> List[structure._Chain]:
    '''
    Calls measure to get chain object within distance A of ligand

    Args:
    cplx - the structure object holding everything in the PDB file
    lig - the ligand object of interest
    distance - maximum distance of atoms to consider

    Returns:
    List of schrodinger.structure._Chain objects that have atoms within distance Angstrom of the ligand
    '''
    atms = measure.get_atoms_close_to_structure(cplx, lig.st, distance)
    chs = list(set([cplx.atom[i].chain for i in atms]))
    return chs

def download_pdbs(base_dir: Union[str, PathLike], pdb_id_list: List[str]) -> None:

    '''
    Downloads a set of PDB IDs into base_dir/raw/ID/ID.pdb

    Args:
    base_dir - data directory for this project
    pdb_id_list - list of RCSB PDB IDs to download

    Returns:
    None
    '''

    base = Path(base_dir)
    raws = base / 'raw'

    if not raws.exists():
        raws.mkdir(parents=True)

    for pdb in pdb_id_list:
        #create a directory for each file
        outpath = raws / (str(pdb))
        outpath.mkdir(parents=True, exist_ok=False)

        pdb_url = f"https://files.rcsb.org/download/{pdb}.pdb"
        command = f"wget {pdb_url} -P {outpath}"

        pc = Popen(command, stdout= PIPE, stderr= PIPE, shell=True)
        o, e = pc.communicate()

        if pc.returncode != 0:
            print(f"Error with {pdb} : {e.decode('utf-8')}")

def convert_pdb_to_maestro(pdb_file: Union[str, PathLike]) -> None:
    '''
    Converts a PDB file to a Maestro File using Schrodinger's structconvert script

    Args:
    pdb_file - Filename or pathlike object of the pdb file to convert

    Returns:
    None
    '''

    pdbf = Path(pdb_file)
    maestro_file = pdbf.with_suffix(".mae")
    #will throw a deprecated warning about -ipdb and -omae since newest versions infer filetype through extension. Keeping as-is
    #since it still works and should be compatible with previous versions
    command = f"$SCHRODINGER/utilities/structconvert -ipdb {pdbf} -omae {maestro_file}"

    pc = Popen(command, stdout= PIPE, stderr= PIPE, shell=True)
    o, e = pc.communicate()

    if pc.returncode != 0:
        print(f"Error with {pdb_file}: {e.decode('utf-8')}")    


def process_structure(filename: Union[str, PathLike], ligname: str) -> None:
    '''
    Pre-processes the raw Maestro File for downstream applications
    This script will identify and keep chains that are in close proximity to the ligand of interest, and will change the ligand's
    chain to 'L' for ease of downstream processing.
    The output is a processed file of the form filename_proc.mae

    Args:
    filename - maestro file to be processed
    ligname - ligand's 3-letter code

    Returns:
    None
    '''
    alpha = list('ABCDEFGHIJKMNOPQRSTUVWXYZ')

    fin = Path(filename)

    #read in maestro file
    with structure.StructureReader(filename) as st:
        cplx = next(st)

    real_ligs = [l for l in analyze.find_ligands(cplx) if l.pdbres.strip() == ligname]
    
    #future fix - some structures fail the find_ligands command, need to figure out why
    if not real_ligs:
        print(f"\nError with {fin.stem}")
    else:
        #sort by b factor if there are multiple copies of the ligand, take the lig with lowest b-factor
        real_ligs.sort(key=lambda x: lig_avg_bfact(x))
        lig_to_take = real_ligs[0]

        chs = get_chains_around_lig(cplx, lig_to_take)
        
        #remove the ligand from the structure, still have a copy in ligs_to_take
        ##TODO: Remove the indices of all ligands that are not the primary ligand. Just extend the atom indexes for all ligands in real_ligs

        l_atoms = lig_to_take.atom_indexes
        cplx.deleteAtoms(l_atoms)

        #make a new structure of the chains around the ligand, rename alphabetically, but don't use chain L. Ligand will be chain L
        chains = [cplx.chain[c] for c in chs]
        for i,c in enumerate(chains):
            c.name = alpha[i]
        
        #create a new structure object with these chains
        newstruc = chains[0].extractStructure()
        if len(chains) > 1:
            for c in chains[1:]:
                newstruc.extend(c.extractStructure())

        #rename ligand
        ligst = lig_to_take.st
        lig_chain_name = list(ligst.chain)[0].name
        ligst.chain[lig_chain_name].name = 'L'
        #add ligand to new structure
        newstruc.extend(ligst)

        #write to outfile
        outfile = fin.parent / (fin.stem + '_proc.mae')
        newstruc.write(str(outfile))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process Maestro Files")
    parser.add_argument("-o", "--output_dir", help="Directory to save files to")
    parser.add_argument("-s", help="csv of structures to download")

    args = parser.parse_args()

    pdb_info = pd.read_csv(str(args.s), header=0, index_col=0)
    pdbs = list(pdb_info.rscb_id)
    #download pdb files
    download_pdbs(args.output_dir, pdbs)
    #convert them to Maestro
    for pdb_file in glob.glob(f"{args.output_dir}/raw/*/*.pdb"):
        convert_pdb_to_maestro(pdb_file)
    #process Maestro to keep chains around ligand, and rename ligand to chain 'L'
    for i,f in enumerate(glob.glob(f"{args.output_dir}/raw/*/*.mae")):
        fpath = Path(f)
        pdb_name = fpath.stem

        if '_proc' in pdb_name:
            #already processed file
            continue
        #because tqdm is for some reason disabled in schrodinger's python... use an old school progress reporter.
        sys.stdout.write(f"\rWorking on {i}/{len(pdb_info)} | {pdb_name}")
        sys.stdout.flush()
        
        info_df = pdb_info[pdb_info['rscb_id'] == str(pdb_name)]

        l3letter = np.array(info_df['lig_3letter'])[0]

        process_structure(f, l3letter)
