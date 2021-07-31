from tqdm import tqdm
import probe_extraction

pdb_codes = open('pdb_lists/unique_pdbs.txt').read().split(',')[:3]
pdb_chains = open('pdb_lists/abb8330_Data_S1.txt').read().splitlines()
pdb_chains = [line for line in pdb_chains if line[0] != '#']

for pdb_code in tqdm(pdb_codes):
    chains = [line[4] for line in pdb_chains if line[:4] == pdb_code]
    pdb_file = 'data/reduce_out/reduce_' + pdb_code + '.pdb'
    probe_extraction.parse_probe(pdb_file, chains, outdir='data/probe_out')
