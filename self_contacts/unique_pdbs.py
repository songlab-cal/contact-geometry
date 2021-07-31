chain_list = 'pdb_lists/abb8330_Data_S1.txt'
chains = open(chain_list, 'r').read().splitlines()
chains = [chain for chain in chains if chain[0] != '#']
unique_pdbs = set([chain[:4] for chain in chains])
open('pdb_lists/unique_pdbs.txt', 'w').write(",".join(unique_pdbs))
