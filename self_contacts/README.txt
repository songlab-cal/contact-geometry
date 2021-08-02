To use this code, you need Python (developed in 3.7.10) as well as numpy, pickle, pandas, and tqdm. You also need Reduce and Probe binaries in your system. Then:

1) Place a text file of pdb chains in pdb_lists/. The file abb8330_Data_S1.txt is the list we use and is the same as in the earlier work by Polizzi and DeGrado.
If you make a new file, follow the same format. The comment lines are not required, but any comment lines should start with #.

2) Place pdb files in data/raw_structs/

3) Run unique_pdbs.py after changing the first line to point to your list

4) Set the path to probe and reduce in probe_extraction.py

5) Run reduce_all.py

6) Run probe_all.py after changing line 5 to point to your list of chains.

7) Run extract_all.py

8) Run combine_dicts.py

9) Run contact_reader.py

This will output unaligned contacts to data/contacts for each contact type and each choice of h-bonding and sidechain-bonding True and False, as pickle files. 
The pickle files each are a tuple with two elements. The first element is NxMx3, where N is the number of observations. 
Each observation consists of the 3 backbone atom coordinates concatenated with the ifg coordinates, so M = 3 + the size of the ifg. 
The second element of the array contains the sidechain coordinates (everything except N, C, Ca) for the same N observations. 
The atoms are (in order) exactly the atoms appearing before "OXT" in data/raw_contacts/aa_seqs.pkl for the aa, and in the same order as those in functional_groups.SMARTS_ATOM_MAP for the ifg.

To convert these pkl files to the .npy arrays that are provided in the dataset upload, you can use the script at the top level of the repo after changing the input and output directories.
