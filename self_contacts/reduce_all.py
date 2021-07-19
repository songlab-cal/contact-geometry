from tqdm import tqdm
import probe_extraction as extract

pdb_codes = open('pdb_lists/unique_pdbs.txt').read().split(',')[:3]
pdb_files = ['data/raw_structs/' + pdb_code + '.pdb' for pdb_code in pdb_codes]
extract._reduce_pdb(pdb_files, outdir='data/reduce_out')
