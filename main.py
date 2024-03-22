from __future__ import print_function
import argparse
import os
import pymol
from build_seq import build_seq
from modify_itp import *

if __name__ == "__main__":
    
    ############################
    ####   ELP Generation   ####
    ############################
    parser = argparse.ArgumentParser(description='Arguments for ELP generation, X1 is hydrophilic and X2 is hydrophobic')
    parser.add_argument('-x1', type=str, help='X1 is hydrophilic, choose from the following: {G, T, S, W, Y, H, E, Q, D, N, K, R}')
    parser.add_argument('-m', type=int, help='Repeat of the fragment 1: (VPGX1G)m')
    parser.add_argument('-x2', type=str, help='X2 is hydrophobic, choose from the following: {I, V, L, F, C, M, A}')
    parser.add_argument('-n', type=int, help='Repeat of the fragment 2: (VPGX2G)n')
    parser.add_argument('-save', dest='save', action='store_true', help='Save pdb and CG models')
    parser.add_argument('--no-save', dest='save', action='store_false', help='Do not save pdb')
    parser.add_argument('-save_dir', type=str, help='The Directory to store the aa and cg pdb, Default: the working directory')
    parser.set_defaults(save=True)

    args = parser.parse_args()

    print('X1 is {} and X2 is {}'.format(args.x1, args.x2))

    fragment1 = 'VPG{}G'.format(args.x1)
    fragment2 = 'VPG{}G'.format(args.x2)
    seq = fragment1*args.m + fragment2*args.n

    ELP_name = '{}{}{}{}'.format(args.x1, args.m, args.x2, args.n)
    whole_ELP_name = '({}){}({}){}'.format(fragment1, args.m, fragment2, args.n)

    dir_path = args.save_dir+"/"+ELP_name
        
    if args.save:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        print('\n### The ELP to be generated is {}, and the sequence is {} ###\n'.format(ELP_name, whole_ELP_name))
        build_seq(ELP_name, seq, ss="polypro")
        pymol.cmd.save(filename="{}/{}_AA.pdb".format(dir_path, ELP_name))
        #print("Directory '{}' was created.".format(dir_path))
        #with open("{}/ELP.dat".format(args.save_dir), 'w') as f:
        #    f.write(seq)

        print('\n### All-Atom ELP Successfully Generated and Saved to {} ! ###\n'.format(dir_path))

    else:
        with open("{}/ELP.dat".format(args.save_dir), 'a') as f:
            f.write(ELP_name+" : "+seq+"\n")
        print('### Only String is Saved to {} ! ###\n'.format(args.save_dir))


    ###########################
    #### ELP Martinization ####
    ###########################
    if args.save:
        cg_command = 'python2 martinize.py -f {}/{}_AA.pdb -o {}/System.top -x {}/elp_cg.pdb'.format(dir_path, ELP_name, dir_path, dir_path)

        os.system(cg_command) 
        # the Protein.itp is always generated when the command is executed
        # we modify it to have a neutral C terminal and then move to the desired directory
        os.system('cp Protein.itp Protein_ori.itp') 
    
        # Only the atom types and charges are changed 
        # The other bonded terms are not perturbed and thus can be copied from state A to state B by gromacs during gompp
        modify_C_ter('Protein.itp')
        add_mass('Protein.itp')
        gen_elp_1('Protein.itp', 'elp_1.itp')
        gen_elp_2('Protein.itp', 'elp_2.itp')
        os.system('mv *itp '+dir_path)
    
