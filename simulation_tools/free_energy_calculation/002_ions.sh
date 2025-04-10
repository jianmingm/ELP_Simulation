#!/bin/bash
gmx_mpi editconf -f bilayer_solvated.pdb -o bilayer_solvated.gro -c 
gmx_mpi grompp -f mdp_files/ions.mdp -c bilayer_solvated.gro -p System.top -o ions.tpr
echo 13 | gmx_mpi genion -s ions.tpr -o bilayer_ions.gro -p System.top -pname NA+ -nname CL -neutral

