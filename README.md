# Protein PDB Torsional Angle Reader
This tool reads a protein pdb file and looks for relevant phi-psi angles. Trajectory files from MD simulations are deconstructed using parse files and are parsed based on the output of the program (Tor.sh) provided. This program identifies relevant Phi-Psi angles by reading the respective protein pdb file and directs the parse.c program to parse the relevant parts of the peptide sequence resulting in several parsed trajectory files displaying different Phi-Psi angles corresponding to steps of the OpenMM simulation. In this case, as an example the protein pdb file contains 5 Alanines in a peptide sequence. The program functions to first create a pdb file by accepting the user's input of the protein's peptide sequence, simulate the protein using the OpenMM Molecular Dynamics program and finally determine the locations of Carbon and Nitrogen atoms involved in the formation of a Phi-Psi angle. There are four files provided: the Tor.sh file reads the pdb file and determines the relevant locations; the parse.c files deconstructs the trajectory file from an MD simulation using the software OpenMM and determines the values of the phi-psi angles corresponding to each step in the simulation; the TrajectorySoftware.bash file can accept the user's input of a peptide sequence and outputs a pdb file containing the peptide sequence and an optimized number of water molecules; 5ALA.pdb contains a sample PDB file that the Tor.sh file can read for the Phi-Psi angles. 
