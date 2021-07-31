#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle as pkl
import numpy as np
import subprocess
from tqdm import tqdm
from decomp.functional_groups import SMARTS_ATOM_MAP
from scoring.utils import BACKBONE_ATOMS, AMINO_ACIDS, VDM_TYPES, IDEAL_ALA_REF

ifg_list = list(SMARTS_ATOM_MAP.keys())
vdm_types = list(set([item[0] + '-' + item[1] for item in VDM_TYPES]))


# In[2]:


in_dir = 'data/dcoid_probe/'
out_dir = 'prot_lig_contacts/'


# In[3]:


# Mobile and target are both N x 3
def kabsch(mobile, target):
    H = mobile.T @ target
    u, s, vh = np.linalg.svd(H)
    d = np.sign(np.linalg.det(vh.T @ u.T))
    R = vh.T @ np.array([[1, 0, 0],[0, 1, 0],[0, 0, d]]) @ u.T
    return R


# In[4]:


result = subprocess.run(['ls ' + in_dir + '*.pkl'], stdout=subprocess.PIPE, shell=True)
dcoid_files = result.stdout.decode("utf-8").splitlines()


# In[5]:


vdms_dcoid_true_true = {vdm_type:[] for vdm_type in vdm_types}
vdms_dcoid_true_false = {vdm_type:[] for vdm_type in vdm_types}
vdms_dcoid_false_true = {vdm_type:[] for vdm_type in vdm_types}
vdms_dcoid_false_false = {vdm_type:[] for vdm_type in vdm_types}

for i in tqdm(range(len(dcoid_files))):
    complex_name = dcoid_files[i].split('/')[-1][:4]
    vdms = pkl.load(open(dcoid_files[i], 'rb'))
    data = vdms['Ligand_1']
    
    for j, vdm in enumerate(data):
        ifg = repr(vdm).rsplit('_', 1)[0]
        if ifg not in SMARTS_ATOM_MAP.keys(): continue
        ifg_coords = vdm.coords.astype('float32')
        
        for contact in vdm.contacts.keys():
            aa = repr(contact)[-5:-2]
            if not aa in AMINO_ACIDS: continue
            backbone_coords = np.array(vdm.contacts[contact]['residue_coords'][:3], dtype='float32')
            #sidechain_coords = np.array(vdm.contacts[contact]['residue_coords'][3:], dtype='float32')
            
            ifg_coord = ifg_coords - backbone_coords.mean(axis=0)
            backbone_coords = backbone_coords - backbone_coords.mean(axis=0)
            ref_bb = IDEAL_ALA_REF - IDEAL_ALA_REF.mean(axis=0)
            R = kabsch(backbone_coords, ref_bb)
            
            ifg_coord = ifg_coord @ R.T
            backbone_coords = backbone_coords @ R.T
            
            bonds = vdm.contacts[contact]['interactions']
            entry_1 = np.concatenate((backbone_coords, ifg_coord), axis=0)
            #entry_1_atoms = np.array(vdm.contacts[contact]['residue_atoms'][:3] + vdm.atoms)
            #sidechain_atoms = np.array(vdm.contacts[contact]['residue_atoms'][3:])
            
            vdm_type = aa + '-' + ifg
            hbonding = (len(bonds['hb']) > 0)
            sidechain_bonding = not all([atom_pair[1] in BACKBONE_ATOMS for atom_pair in bonds['hb'] | bonds['cc']])
            
            if hbonding:
                if sidechain_bonding:
                    vdms_dcoid_true_true[vdm_type].append(entry_1)
                else:
                    vdms_dcoid_true_false[vdm_type].append(entry_1)
            else:
                if sidechain_bonding:
                    vdms_dcoid_false_true[vdm_type].append(entry_1)
                else:
                    vdms_dcoid_false_false[vdm_type].append(entry_1)

for vdm_type in vdm_types:
    vdms_dcoid_true_true[vdm_type] = np.array(vdms_dcoid_true_true[vdm_type])
    vdms_dcoid_true_false[vdm_type] = np.array(vdms_dcoid_true_false[vdm_type])
    vdms_dcoid_false_true[vdm_type] = np.array(vdms_dcoid_false_true[vdm_type])
    vdms_dcoid_false_false[vdm_type] = np.array(vdms_dcoid_false_false[vdm_type])


# In[6]:


for vdm_type in vdm_types:
    np.save(out_dir+vdm_type+'-h_bonding_True-sidechain_bonding_True.npy', vdms_dcoid_true_true[vdm_type])
    np.save(out_dir+vdm_type+'-h_bonding_True-sidechain_bonding_False.npy', vdms_dcoid_true_false[vdm_type])
    np.save(out_dir+vdm_type+'-h_bonding_False-sidechain_bonding_True.npy', vdms_dcoid_false_true[vdm_type])
    np.save(out_dir+vdm_type+'-h_bonding_False-sidechain_bonding_False.npy', vdms_dcoid_false_false[vdm_type])


# In[ ]:




