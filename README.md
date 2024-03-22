### ELP Generation built on PyMOL

#### Prerequisite #### 
The code is built upon python2 and the build_sequence code by Robert L. Campbell.

The following packages are required

-- python2 : it is because PyMOL and its third party codes are majorly compatible with python2

-- PyMOL (version 2): can be installed by `conda install -c schrodinger pymol` (version 3 is not compatible with build_sequence code)

-- martinize.py : the python2 version (http://cgmartini.nl/index.php/tools2/proteins-and-bilayers/204-martinize)

#### Usage

The help message can be displayed by `python main_genelp.py -h`

example: `python main.py -x1 G -m 2 -x2 L -n 3 -save_dir . -save`

example: `python main.py -x1 G -m 2 -x2 L -n 3 -save_dir . --no-save`


usage: python main_genelp.py [-h] [-x1 X1] [-m M] [-x2 X2] [-n N] [-save or --no-save] [-save_dir SAVE_DIR]


Arguments for ELP generation, X1 is hydrophilic and X2 is hydrophobic

optional arguments:

  -h, --help          show this help message and exit

  -x1 X1              X1 is hydrophilic, choose from the following: {G, T, S, W, Y, H, E, Q, D, N, K, R}

  -m M                Repeat of the fragment 1: (VPGX1G)m

  -x2 X2              X2 is hydrophobic, choose from the following: {I, V, L, F, C, M, A}

  -n N                Repeat of the fragment 2: (VPGX2G)n

  -save               Save pdb and CG models

  --no-save           Do not save pdb

  -save_dir SAVE_DIR  The Directory to store the aa and cg pdb, Default: the
                      working directory


#### modify_itp.py

This is used to modify the C terminal to neutral, and to generate the itp files with state A and B for alchemical simulations.

#### batch_gen.sh

The gen.sh takes 4 arguments that include the number of the first residue, the number of the second residue, whether to save pdbs and the save directory

Then it will run the main.py and then output the sequence into ELP.dat in the specified directories

If set to `-save`, this will also generate the cg models and then create the elp_1.itp and elp_2.itp for subsequent alchemical transfer simulations. 

... will continue to generate codes to automate the whole alchemical simulation method
