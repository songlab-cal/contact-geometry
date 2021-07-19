from tqdm import tqdm
import probe_extraction

pdb_codes = open('pdb_lists/unique_pdbs.txt').read().split(',')[:3]
outdir = 'data/extract_out'
pdbdir = 'data/reduce_out'
probedir = 'data/probe_out'

for pdb_code in tqdm(pdb_codes):
    probe_file = probedir + '/' + pdb_code + '.csv'
    probe_extraction.extract_coords(probe_file, pdbdir, outdir=outdir)
