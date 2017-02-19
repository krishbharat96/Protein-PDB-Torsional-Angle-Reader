# Protein-PDB-Torsional-Angle-Reader
This tool reads a protein pdb file and looks for relevant phi-psi angles. Trajectory files from MD simulations are deconstructed using 
parse files and are parsed using the program (Tor.sh) provided. This program identifies relevant Phi-Psi angles by reading the respective
protein pdb file (which, in this case, is 5 Alanines in a peptide sequenvce and determining the locations of Carbon and Nitrogen atoms 
involved in the formation of a Phi-Psi angle. There are three files provided: the Tor.sh file reads the pdb file and determines the
relevant locations; the parse.c files deconstructs the trajectory file from an MD simulation using the software OpenMM and determines the
values of the phi-psi angles corresponding to each step in the simulation; the 5ALA.pdb contains a sample PDB file that the Tor.sh file
can read for the Phi-Psi angles. 
