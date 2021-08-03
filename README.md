# contact-geometry

This repository provides the code needed to replicate results from the paper "Transferability of Geometric Patterns from Protein Self-Interactions to Protein-Ligand Interactions" submitted to PSB.

Fully pre-processed data for both protein self-contacts and protein-ligand contacts are also available in the preprocessed folder

Our preprocessing pipeline can be run on your own list of PDB files to extract vdMs from structures of interest. All you need is a CSV file of pdb ids and the ligand of interest's 3-letter code. This does, however, require access to SCHRODINGER, and, in particular, a license for EPIK and IMPREF to properly identify hydrogen bonds. If you can find another way to correctly protonate your ligand, you can skip this step.

Parts of this code are modified from the Combs repository under the MIT license: https://github.com/npolizzi/Combs/blob/master/LICENSE.txt.

This repository is also available under the MIT license.
