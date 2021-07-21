# contact-geometry

This repository provides the code needed to replicate results from the paper XXXXX submitted to PSB.

Our preprocessing pipeline can be run on your own list of PDB files to extract vdMs from structures of interest. All you need is a CSV file of pdb ids and the ligand of interest's 3-letter code. This does, however, require access to SCHRODINGER, and, in particular, a license for EPIK and IMPREF to properly identify hydrogen bonds. If you can find another way to correctly protonate your ligand, you can skip this step.

For the results in the paper, we provide pre-processed pkl files for each PDB structure. 
