#!/bin/bash
#SBATCH --job-name=mem-equi
#SBATCH --partition=andrewferguson-gpu
#SBATCH --time=8:00:00
#SBATCH --output=equi.out
#SBATCH --error=equi.err
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=20G
#SBATCH --account=pi-andrewferguson

source /home/mjianming/gromacs-2023.1/install/bin/GMXRC
module load cuda/12.0 gcc/10.2.0  openmpi

gmx_mpi grompp -f mdp_files/em.mdp -c bilayer_ions.gro -p System.top -o em.tpr -maxwarn 2
time gmx_mpi mdrun -v -deffnm em -ntomp 8

gmx_mpi grompp -f mdp_files/equilibration_nvt.mdp -c em.gro -p System.top -o nvt.tpr -maxwarn 2
time gmx_mpi mdrun -v -deffnm nvt -ntomp 8

gmx_mpi grompp -f mdp_files/equilibration_npt.mdp -c nvt.gro -p System.top -t nvt.cpt -o npt.tpr -maxwarn 2
time gmx_mpi mdrun -v -deffnm npt -ntomp 8

gmx_mpi grompp -f mdp_files/dynamics.mdp -c npt.gro -t npt.cpt -p System.top -o md.tpr -maxwarn 2
time gmx_mpi mdrun -v -deffnm md -ntomp 8
