import os
import pandas as pd
from collections import defaultdict
import numpy as np
import pickle as pkl

# CHANGE TO MATCH YOUR ENVIRONMENT
probe_path = '~/probe/probe.2.16.130520.linuxi386*'
reduce_path = '~/probe/reduce.3.23.130521.linuxi386*'

def _reduce_pdb(pdb_files, outdir):
	for pdb_file in pdb_files:
		suffix = pdb_file.split('/')[-1]
		cmd = [reduce_path, '-FLIP', pdb_file, '>', outdir+'/reduce_'+suffix]
		cmd = ' '.join(cmd)
		os.system(cmd)

def _make_cmd(chainname1, pdb_file, chainname2='',
	probe_sel_criteria='blt40 ogt99 not metal',
	outdir=None):
	chain1 = 'CHAIN' + chainname1
	chain2 = ''
	if chainname2 != '':
		probe_action = '-B'
		chain2 = 'CHAIN' + chainname2
	else:
		probe_action = '-SE'

	cmd = [probe_path, '-U -CON -DOCHO -MC -DE32',
		probe_action, '"', chain1, probe_sel_criteria, '"']

	if chainname2 != '':
		chain2_cmd = ['"', chain2, probe_sel_criteria, '"']
		cmd.extend(chain2_cmd)

	cmd.append(pdb_file)

	if outdir is not None:
		if outdir[-1] != '/': outdir += '/'
		outfile = ''.join([outdir, pdb_file.split('/')[-1][7:-4],
			'_', chainname1, chainname2, '.probe'])
		cmd.extend(['>', outfile])

	return ' '.join(cmd)

def _parse_probe_line(line, chainname1, chainname2):
	if chainname2 == '': chainname2 = chainname1

	spl = line.split(':')[1:]
	INTERACTION = spl[1]
	RESNUM1 = int(spl[2][2:6])
	RESNAME1 = spl[2][6:10].strip()
	NAME1 = spl[2][10:15].strip()
	ATOMTYPE1 = spl[12]
	RESNUM2 = int(spl[3][2:6])
	RESNAME2 = spl[3][6:10].strip()
	NAME2 = spl[3][10:15].strip()
	ATOMTYPE2 = spl[13]
	if RESNAME1 == 'HOH': NAME1 = 'O'
	if RESNAME2 == 'HOH': NAME2 = 'O'

	return (INTERACTION, chainname1,
		RESNUM1, RESNAME1, NAME1, ATOMTYPE1,
		chainname2, RESNUM2, RESNAME2,
		NAME2, ATOMTYPE2)

def parse_probe(pdb_file, chains,
	probe_sel_criteria='blt40 ogt99 not metal',
	outdir=None):
	col_names = ['interaction', 'chain1', 'resnum1', 'resname1', 'name1',
		'atomtype1', 'chain2', 'resnum2', 'resname2', 'name2', 'atomtype2']
	probe_df = pd.DataFrame(columns=col_names)
	for i, chainname1 in enumerate(chains):
		chains2 = ['']
		if i < len(chains) - 1: chains2 += chains[i+1:]
		for chainname2 in chains2:
			cmd = _make_cmd(chainname1, pdb_file, chainname2,
				probe_sel_criteria)
			probe_dict = defaultdict(list)
			probe_data = []
			with os.popen(cmd) as probefile:
				for line in probefile:
					line_data = _parse_probe_line(line, chainname1, chainname2)
					probe_dict[line_data[1:]].append(line_data[0])
			for info, interactions in probe_dict.items():
				if 'bo' in interactions: continue
				elif 'so' in interactions: continue
				elif 'hb' in interactions: interaction = 'hb'
				elif 'cc' in interactions: interaction = 'cc'
				elif 'wc' in interactions: continue
				if info[0] == info[5] and abs(int(info[1]) - int(info[6])) < 8: continue
				data = [interaction]
				data.extend(info)
				probe_data.append(data)
			probe_df = probe_df.append(pd.DataFrame(probe_data, columns=col_names))
	if outdir is None: outdir == '.'

	probe_df.to_csv(outdir+'/'+pdb_file.split('/')[-1][7:-4]+'.csv')


# Probe_df_file should be named *.csv, and the corresponding 
# pdb file outputted by reduce should exist at pdb_dir/reduce_*.pdb
def extract_coords(probe_df_file, pdb_dir, outdir):
	file_tag = probe_df_file.split('/')[-1][:-4]
	if pdb_dir[-1] != '/': pdb_dir += '/'
	data = dict()

	interaction_list = open(probe_df_file).read().splitlines()[1:]
	with open(pdb_dir+'reduce_'+file_tag+'.pdb', encoding='utf8', errors='ignore') as f:
		pdb_file = f.read().splitlines()
	pdb_file = [line for line in pdb_file if len(line) > 4 and line[:4] == 'ATOM']

	interaction_list = [interaction.split(',') for interaction in interaction_list]
	chain_res_pairs = dict()
	for entry in interaction_list:
		chain1 = entry[2]; chain2 = entry[7]
		resnum1 = int(entry[3]); resnum2 = int(entry[8])
		resname1 = entry[4]; resname2 = entry[9]
		name1 = entry[5]; name2 = entry[10]
		if chain1 == chain2 and resnum2 < resnum1:
			tmp = resnum1
			resnum1 = resnum2
			resnum2 = tmp
			tmp = resname1
			resname1 = resname2
			resname2 = tmp
			tmp = name1
			name1 = name2
			name2 = tmp
		interaction_key = (chain1, chain2, resnum1, resnum2, resname1, resname2)
		if interaction_key not in chain_res_pairs:
			chain_res_pairs[interaction_key] = {'hb':set(),'cc':set()}
		if entry[1] == 'hb':
			chain_res_pairs[interaction_key]['hb'].add((name1, name2))
		elif entry[1] == 'cc':
			chain_res_pairs[interaction_key]['cc'].add((name1, name2))

	for interaction in chain_res_pairs:
		(chain1, chain2, resnum1, resnum2, resname1, resname2) = interaction

		coords_1 = [line for line in pdb_file if line[21] == chain1]
		coords_1 = [line for line in coords_1 if int(line[22:26]) == resnum1]
		if (len(coords_1) == 0): continue
		atom_names_1 = [line[12:16] for line in coords_1]
		coords_1 = [[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in coords_1]
		coords_1 = np.array(coords_1)
		all_info_1 = (resname1, atom_names_1, coords_1)

		coords_2 = [line for line in pdb_file if line[21] == chain2]
		coords_2 = [line for line in coords_2 if int(line[22:26]) == resnum2]
		if (len(coords_2) == 0): continue
		atom_names_2 = [line[12:16] for line in coords_2]
		coords_2 = [[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in coords_2]
		coords_2 = np.array(coords_2)
		all_info_2 = (resname2, atom_names_2, coords_2)

		bond_list = chain_res_pairs[interaction]

		if resname1 < resname2:
			identifier = (resname1, resname2)
			if identifier not in data: data[identifier] = []
			data[identifier].append([all_info_1, all_info_2, bond_list])
		else:
			identifier = (resname2, resname1)
			bond_list = {key1: {val2[::-1] for val2 in bond_list[key1]} for key1 in bond_list}
			if identifier not in data: data[identifier] = []
			data[identifier].append([all_info_2, all_info_1, bond_list])
	with open(outdir+'/'+file_tag+'.pkl', 'wb') as f:
		pkl.dump(data, f)
