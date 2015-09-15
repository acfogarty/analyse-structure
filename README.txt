***************

Usage:
./analyseStructure.exe < input-analyseStructure

***************
Sample input-analyseStructure:

/*** input file for analyseStructure.exe ***/
/*** trajectory file ***/
GRO !coordfiletype, GRO, PDB or DCD
4 !number of trajectory files
../trj-run1.gro 946 !name of traj file, nframes in traj file
../trj-run2.gro 939 
../trj-run3.gro 938
../trj-run4.gro 951
32854 !natoms
0.25 !timestep btwn frames in ps
4 !nskip, frequency at which to skip over frames in traj files during analysis
/*** distances, angles and dihedrals to calculate ***/
7 !nmeasurements
/*** output filename, type [BON, ANG or DIH], indices of atoms involved ***/
dist_ASN59CA_TRP63CA.dat    BON 256 298 # comment
dist_ALA107CA_TRP63CA.dat   BON 349 298    
dist_ALA107CA_ASN59CA.dat   BON 349 256  
dist_ASN59CA_TRP62CA.dat    BON 256 274 
dist_ALA107CA_TRP62CA.dat   BON 349 274 
dist_TRP63CA_TRP62CA.dat    BON 298 274 
dih_TRP62CACBCGCD2.dat      DIH 274 276 279 293 
