#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
import pickle as pkl
import numpy as np
from tqdm import tqdm
from scoring.utils import IDEAL_ALA_REF, VDM_TYPES
from decomp.functional_groups import SMARTS_ATOM_MAP
ifg_list = list(SMARTS_ATOM_MAP.keys())
vdm_types = list(set([item[0] + '-' + item[1] for item in VDM_TYPES]))


# In[10]:


in_dir = '../pdb_contacts/'
out_dir = 'prot_self_contacts/'


# In[11]:


files = os.listdir(in_dir)


# In[12]:


# Mobile and target are both N x 3
def kabsch(mobile, target):
    H = mobile.T @ target
    u, s, vh = np.linalg.svd(H)
    d = np.sign(np.linalg.det(vh.T @ u.T))
    R = vh.T @ np.array([[1, 0, 0],[0, 1, 0],[0, 0, d]]) @ u.T
    return R

def align(data):
    if data.size == 0: return data.copy()
    contacts = data.copy()
    contacts -= contacts[:, :3, :].mean(axis=1, keepdims=True)
    ref_backbone = IDEAL_ALA_REF - IDEAL_ALA_REF.mean(axis=0)
    for i in range(contacts.shape[0]):
        R = kabsch(contacts[i, :3, :], ref_backbone)
        contacts[i, :, :] = contacts[i, :, :] @ R.T
    return contacts


# In[13]:


data_true_true = {}
data_true_false = {}
data_false_true = {}
data_false_false = {}

for vdm_type in tqdm(vdm_types):
    match_files = [file for file in files if file.split('-h_bonding')[0] == vdm_type]
    for file in match_files:
        data = pkl.load(open(in_dir+file, 'rb'))[0]
        data = align(data)
        
        h_bonding = ('-h_bonding_True' in file)
        sidechain = ('-sidechain_bonding_True' in file)
        
        if h_bonding and sidechain:
            np.save(out_dir+vdm_type+'-h_bonding_True-sidechain_bonding_True.npy', data)
        elif h_bonding:
            np.save(out_dir+vdm_type+'-h_bonding_True-sidechain_bonding_False.npy', data)
        elif sidechain:
            np.save(out_dir+vdm_type+'-h_bonding_False-sidechain_bonding_True.npy', data)
        else: 
            np.save(out_dir+vdm_type+'-h_bonding_False-sidechain_bonding_False.npy', data)


# In[ ]:




