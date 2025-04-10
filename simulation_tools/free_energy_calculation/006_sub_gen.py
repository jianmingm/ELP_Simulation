import os

parent_dir = os.path.basename(os.getcwd())

for num in range(0, 18):
    job_name = f"{parent_dir}-r1-{num}"

    command = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=gm4-pmext
#SBATCH --account=pi-andrewferguson
#SBATCH --qos=gm4-cpu
#SBATCH --time=16:00:00
#SBATCH --output=r1.out
#SBATCH --error=r1.err
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=2500

module load cuda/11.5 cmake/3.19 gcc/10.2.0
source /project2/andrewferguson/mjianming/Softwares/gromacs-2023.1/install/bin/GMXRC

FREE_ENERGY=`pwd`
echo "Free energy home directory set to $FREE_ENERGY"
MDP=$FREE_ENERGY/MDP
echo ".mdp files are stored in $MDP"

ntmpi=10

if [ -f "equi.gro" ]; then
    echo "Equi has done"
else
gmx_mpi grompp -f $MDP/em.mdp -c $FREE_ENERGY/init.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -o min.tpr -maxwarn 3
gmx_mpi mdrun -v -deffnm min -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

gmx_mpi grompp -f $MDP/equi.mdp -c $FREE_ENERGY/min.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top  -o equi.tpr -maxwarn 3
gmx_mpi mdrun -v -deffnm equi -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

fi

LAMBDA={num}

mkdir Lambda_${{LAMBDA}}; cd Lambda_${{LAMBDA}}
mkdir Production_MD
cd Production_MD

gmx_mpi grompp -f $MDP/Production_MD/md_${{LAMBDA}}.mdp -c ../../equi.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -t ../../equi.cpt -o md${{LAMBDA}}.tpr -maxwarn 3

gmx_mpi mdrun -v -deffnm md${{LAMBDA}} -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

echo "Ending. Job completed for lambda = ${{LAMBDA}}"

cd $FREE_ENERGY

exit;
"""

    script_filename = f"./restraint-1/sub_r1_{num}.sbatch"
    with open(script_filename, 'w') as script_file:
        script_file.write(command)

for num in range(0, 18):
    job_name = f"{parent_dir}-r2-{num}"

    command = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=gm4-pmext
#SBATCH --account=pi-andrewferguson
#SBATCH --qos=gm4-cpu
#SBATCH --time=16:00:00
#SBATCH --output=r1.out
#SBATCH --error=r1.err
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=2500

module load cuda/11.5 cmake/3.19 gcc/10.2.0
source /project2/andrewferguson/mjianming/Softwares/gromacs-2023.1/install/bin/GMXRC

FREE_ENERGY=`pwd`
echo "Free energy home directory set to $FREE_ENERGY"
MDP=$FREE_ENERGY/MDP
echo ".mdp files are stored in $MDP"

ntmpi=10

if [ -f "equi.gro" ]; then
    echo "Equi has done"
else
gmx_mpi grompp -f $MDP/em.mdp -c $FREE_ENERGY/init.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -o min.tpr -maxwarn 3
gmx_mpi mdrun -v -deffnm min -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

gmx_mpi grompp -f $MDP/equi.mdp -c $FREE_ENERGY/min.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top  -o equi.tpr -maxwarn 3
gmx_mpi mdrun -v -deffnm equi -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

fi

LAMBDA={num}

mkdir Lambda_${{LAMBDA}}; cd Lambda_${{LAMBDA}}
mkdir Production_MD
cd Production_MD

gmx_mpi grompp -f $MDP/Production_MD/md_${{LAMBDA}}.mdp -c ../../equi.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -t ../../equi.cpt -o md${{LAMBDA}}.tpr -maxwarn 3

gmx_mpi mdrun -v -deffnm md${{LAMBDA}} -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

echo "Ending. Job completed for lambda = ${{LAMBDA}}"

cd $FREE_ENERGY

exit;
"""

    script_filename = f"./restraint-2/sub_r2_{num}.sbatch"
    with open(script_filename, 'w') as script_file:
        script_file.write(command)

for num in range(0, 41):
    job_name = f"{parent_dir}-al-{num}"

    command = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=gm4-pmext
#SBATCH --account=pi-andrewferguson
#SBATCH --qos=gm4-cpu
#SBATCH --time=16:00:00
#SBATCH --output=r1.out
#SBATCH --error=r1.err
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=2500

module load cuda/11.5 cmake/3.19 gcc/10.2.0
source /project2/andrewferguson/mjianming/Softwares/gromacs-2023.1/install/bin/GMXRC

FREE_ENERGY=`pwd`
echo "Free energy home directory set to $FREE_ENERGY"
MDP=$FREE_ENERGY/MDP
echo ".mdp files are stored in $MDP"

ntmpi=10

LAMBDA={num}

    mkdir Lambda_$LAMBDA
    cd Lambda_$LAMBDA

    ##############################
    # ENERGY MINIMIZATION STEEP  #
    ##############################
    echo "Starting minimization for lambda = $LAMBDA..." 

    mkdir EM
    cd EM

    # Iterative calls to grompp and mdrun to run the simulations

    gmx_mpi grompp -f $MDP/EM/em_steep_$LAMBDA.mdp -c $FREE_ENERGY/init.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -o min$LAMBDA.tpr -maxwarn 3

    gmx_mpi mdrun -v -deffnm min$LAMBDA -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

    sleep 5

    #####################
    # NVT EQUILIBRATION #
    #####################
    echo "Starting constant pressure equilibration..."

    cd ../
    mkdir NVT
    cd NVT

    gmx_mpi grompp -f $MDP/NVT/nvt_$LAMBDA.mdp -c ../EM/min$LAMBDA.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -o nvt$LAMBDA.tpr -maxwarn 3

    gmx_mpi mdrun -v -deffnm nvt$LAMBDA -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

    echo "Constant volume equilibration complete."

    sleep 5

    #####################
    # NPT EQUILIBRATION #
    #####################
    echo "Starting constant pressure equilibration..."

    cd ../
    mkdir NPT
    cd NPT

    gmx_mpi grompp -f $MDP/NPT/npt_$LAMBDA.mdp -c ../NVT/nvt$LAMBDA.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -t ../NVT/nvt$LAMBDA.cpt -o npt$LAMBDA.tpr -maxwarn 3

    gmx_mpi mdrun -v -deffnm npt$LAMBDA -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

    echo "Constant pressure equilibration complete."

    sleep 5

    #################
    # PRODUCTION MD #
    #################
    echo "Starting production MD simulation..."

    cd ../
    mkdir Production_MD
    cd Production_MD

    gmx_mpi grompp -f $MDP/Production_MD/md_$LAMBDA.mdp -c ../NPT/npt$LAMBDA.gro -n $FREE_ENERGY/index.ndx -p $FREE_ENERGY/topol.top -t ../NPT/npt$LAMBDA.cpt -o md$LAMBDA.tpr -maxwarn 2

    gmx_mpi mdrun -v -deffnm md$LAMBDA -ntomp 1 -nt $ntmpi -ntmpi $ntmpi

    echo "Production MD complete."

    # End
    echo "Ending. Job completed for lambda = $LAMBDA"

exit;
"""

    script_filename = f"./alchem/sub_alchem_{num}.sbatch"
    with open(script_filename, 'w') as script_file:
        script_file.write(command)

