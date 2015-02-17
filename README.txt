***************

Usage:
./analyseStructure.exe < input-analyseStructure

***************
Sample input-analyseStructure:

/*** input file for analyseStructure.exe ***/
/*** trajectory file ***/
../trj.gro !grofilename
32648 !natoms
100 !nframes
5 !nskip (only analyse every 5th frame)
/*** distances,angles and dihedrals to calculate ***/
3 !nmeasurements
/* format: type (BON/ANG/DIH), indices of atoms involved */
bond_1-2.dat BON 1 2
angle_1-2-3.dat ANG 1 2 3
dihedral_1-2-3-4.dat DIH 1 2 3 4

