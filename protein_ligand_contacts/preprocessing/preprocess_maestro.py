import argparse
from os import PathLike
import glob
from pathlib import Path
from subprocess import Popen, PIPE
from typing import Union

def prep_file(infile: Union[str, PathLike]) -> Union[str, PathLike]:
    '''
    Rus the prepwizard on the input Maestro file.

    Args:
    infile - File to run prep wizard on
    Returns:
    name of the saved outfile
    '''

    in_f = Path(infile)
    out_f = in_f.parent / (in_f.stem + '_out.mae')
    outfile = str(out_f.relative_to(in_f.parent))

    command = f"$SCHRODINGER/utilities/prepwizard -WAIT -rehtreat -watdist 0 {in_f} {outfile}"
    pc = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    o, e, = pc.communicate()
    
    return out_f

def run_poseviewer(infile: Union[str, PathLike], hbond_max:float = 3.0) -> None:
    
    '''
    Runs poseviewer_interactions to identify all interactins involving the ligand. Ligand is identified as chain 'L' from the previous
    structure import steps.

    Args:
    infile - output of prepwizard (so that H-bonds are correctly identified)
    hbond_max - maximum distance to still count as a hydrogen bond

    Returns:
    None
    '''

    contacts = 'all' #find all contacts
    lig_asl = 'c. L' #ligand ASL - chain L

    command = f'$SCHRODINGER/run poseviewer_interactions.py -csv -hbond_max_dist {hbond_max} -hbond-heavy -contacts {contacts} -lig_asl "{lig_asl}" {infile}'
    pc = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    o, e = pc.communicate()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess Subsets of Mae Files")
    parser.add_argument("-fp", help="File Prefix")
    
    args = parser.parse_args()
    #ran in batches based on the first digit of the PDB ID. This was to get around the multiprocessing lock from prepwizard
    #this does, however, require multiple SCHRODINGER license files.
    for f in glob.glob(f"./{args.fp}[A-Z0-9][A-Z0-9][A-Z0-9]_proc.mae"):
        prepped_f = prep_file(f)
        run_poseviewer(prepped_f)

