import mdtraj as md
import numpy as np
import os

single_elp = md.load('elp_cg.pdb')
single_elp_atom_number = single_elp.topology.n_atoms
command = f"echo \"a 1-{single_elp_atom_number}\nq\n\" | gmx_mpi make_ndx -f md.gro"
os.system(command)

if os.path.exists('index.ndx'):
    ### extract one ELP out of md.gro and reduce it into zero box
    os.system("echo 15 | gmx_mpi trjconv -f md.gro -s md.tpr -n index.ndx -o elp_single.gro")
    os.system("gmx_mpi editconf -f elp_single.gro -box 0 0 0 -o equi_elp_single.gro; rm elp_single.gro")

    ### Get the length of one equilibrated ELP molecule 
    elp = md.load('equi_elp_single.gro')
    elp_atom_number = elp.topology.n_atoms
    
    head = elp.topology.select('resid 0 and name BB')[0]
    tail = elp.topology.select(f'resid {elp.n_residues - 1} and name BB')[0]
    atom_pairs = np.array([[head, tail]])
    length = md.compute_distances(elp, atom_pairs)[0][0]  
    length = round(length,2)

    ### Get the highest Z position of upper plane of ELP membrane
    bilayer = md.load('md.gro')
    bilayer_atoms = bilayer.topology.select('protein')
    Z_positions = bilayer.xyz[0, bilayer_atoms, 2]
    Z_high = max(Z_positions)  # Find the highest Z position
    
    insert_XY_pos = round(bilayer.unitcell_lengths[0][0]/2,2)
    insert_Z_pos = round(Z_high + 3 + length/2, 2)  # Calculate the insertion Z position

    
    ### Write the position.dat for inserting ELPs into water phase
    with open('position.dat', 'w') as f:
        f.write(f'{insert_XY_pos} {insert_XY_pos} {insert_Z_pos}\n')
    
    os.system("echo 14 | gmx_mpi insert-molecules -f md.gro -ci equi_elp_single.gro -ip position.dat -o init.gro -rot none -replace W > tmp 2>&1")
    
    ### Get the number of replaced W moelcules
    with open('tmp', 'r') as file:
        for line in file:
            if 'Replaced' in line and 'residues' in line:
                num_residues = int(line.split()[1])
                break
    #os.system('rm tmp')
    #### Modify the System.top file
    with open('System.top', 'r') as file:
        lines = file.readlines()  # Read all lines into a list

    # Find and update the lines for W, elp_1, and elp_2
    for i, line in enumerate(lines):
        if 'W ' in line:
            parts = line.split()
            ori_num = parts[-1]
            parts[-1] = str(int(parts[-1]) - num_residues)
            lines[i] = ' '.join(parts) + ' ; '+ori_num+' \n'
        if 'Protein ' in line:
            parts = line.split() 
            ori_num = parts[-1]
            parts[-1] = str(int(parts[-1]) - 1)
            lines[i] = ' '.join(parts) + ' ; '+ori_num+' \n'
        elif '[ molecules ]' in line:
            lines.insert(i + 2, f'elp_1         1 \n')  # Add line for elp_1
            lines.insert(i + 6, f'elp_2         1 \n')  # Add line for elp_2
    with open('topol.top', 'w') as file:
        file.writelines(lines)
        
    #### Print information for makeing a new index file
    new_sys = md.load('init.gro')
    sys_atom = new_sys.topology.n_atoms
    print("##### Remake an index file using the following information! ##### \n")
    print(f" elp_1 : a 1-{elp_atom_number}")
    print(f" elp_2 : a {sys_atom-elp_atom_number+1}-{sys_atom}")
    print(f" elp_1_hydrophilic : a 1-10")
    print(f" elp_1_hydrophobic : a {elp_atom_number-9}-{elp_atom_number}")
    print(f" elp_2_hydrophilic : a {sys_atom-elp_atom_number+1}-{sys_atom-elp_atom_number+10}")
    print(f" elp_2_hydrophobic : a {sys_atom-9}-{sys_atom}")
    print(f" bilayer : 1&!15&!16 \n")
    print("##### End of index file definition ##### \n")
    
#    print("##### Run the Following ##### \n")
#    print(f"gmx_mpi make_ndx -f init.gro")
#    print(f"a 1-{elp_atom_number}")
#    print(f"a {sys_atom-elp_atom_number+1}-{sys_atom}")
#    print(f"a 1-10")
#    print(f"a {elp_atom_number-9}-{elp_atom_number}")
#    print(f"a {sys_atom-elp_atom_number+1}-{sys_atom-elp_atom_number+10}")
#    print(f"a {sys_atom-9}-{sys_atom}")
#    print(f"1&!15&!16")
#    print("name 15 elp_1")
#    print("name 16 elp_2")
#    print("name 17 elp_1_hydrophilic")
#    print("name 18 elp_1_hydrophobic")
#    print("name 19 elp_2_hydrophilic")
#    print("name 20 elp_2_hydrophobic")
#    print("name 21 bilayer")
#    print("q")
    command2 = f" echo \"a 1-{elp_atom_number}\na {sys_atom-elp_atom_number+1}-{sys_atom}\na 1-10\na {elp_atom_number-9}-{elp_atom_number}\na {sys_atom-elp_atom_number+1}-{sys_atom-elp_atom_number+10}\na {sys_atom-9}-{sys_atom}\n1&!15&!16\nname 15 elp_1\nname 16 elp_2\nname 17 elp_1_hydrophilic\nname 18 elp_1_hydrophobic\nname 19 elp_2_hydrophilic\nname 20 elp_2_hydrophobic\nname 21 bilayer\nq\n\" | gmx_mpi make_ndx -f init.gro"
    os.system(command2)
    print("##### New index file generated ##### \n")    

else:
    print("Make an index file first to denote one ELP molecule!\n")
    
