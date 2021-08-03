Data upload for the paper "Transferability of Geometric Patterns from Protein Self-Interactions to Protein-Ligand Interactions" submitted to PSB

The folders *_contacts each contain 1600 .npy files. For each of the 400 contact types examined in our paper, there is a file for each selection of h-bonding True/False and sidechain-bonding True/False.

Each file contains a numpy array of size N x M x 3, where N is the number of observed contacts for that contact type in the respective dataset. 
Each observation consists of the 3 backbone atom coordinates concatenated with the ifg coordinates, so M = 3 + the size of the ifg.
The atoms are (in order) exactly the atoms appearing before "OXT" in self_contacts/data/raw_contacts/aa_seqs.pkl for the aa, and in the same order as those in self_contacts/functional_groups.SMARTS_ATOM_MAP for the ifg.
These files can be found in the GitHub repo for the paper, https://github.com/songlab-cal/contact-geometry

The folder pdb_lists contains the list of PDB structures used to generate each dataset. The list of pdb chains for self-contacts was originally published by Polizzi and DeGrado [1], while the list of pdb files for protein-ligand contacts is the actives portion of the DCOID dataset [2].

[1] Polizzi, N. F., & DeGrado, W. F. (2020). A defined structural unit enables de novo design of small-molecule–binding proteins. Science, 369(6508), 1227-1233.
[2] Adeshina, Y. O., Deeds, E. J., & Karanicolas, J. (2020). Machine learning classification can reduce false positives in structure-based virtual screening. Proceedings of the National Academy of Sciences, 117(31), 18477-18488.