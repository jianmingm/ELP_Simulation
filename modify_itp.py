def modify_C_ter(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Find the line number of the line above [ bonds ]
    bonds_line_index = lines.index('[ bonds ]\n')

    # Find the line above the [ bonds ] section
    line = lines[bonds_line_index-2].strip()
    if 'Qa' in line and '-1.0000' in line:
        # Modify the line
        modified_line = line.replace('Qa', 'P5').replace('-1.0000', ' 0.0000')
        lines[bonds_line_index-2] = '   ' + modified_line + '\n'  # Update the line in the list

        # Write the modified lines back to the file
        with open(filename, 'w') as file:
            file.writelines(lines)

def add_mass(filename):
    atom_mass_dict={'P5': 72, 'P4': 72, 'P3': 72, 'P2': 72, 'P1': 72, 'Nda': 72, 'Nd': 72, 'Na': 72, 'N0': 72, 'C5': 72,
    'C4': 72, 'C3': 72, 'C2': 72, 'C1': 72, 'Qda': 72, 'Qd': 72, 'Qa': 72, 'Q0': 72, 'SP5': 45, 'SP4': 45,
    'SP3': 45, 'SP2': 45, 'SP1': 45, 'SNda': 45, 'SNd': 45, 'SNa': 45, 'SN0': 45, 'SC5': 45, 'SC4': 45,
    'SC3': 45, 'SC2': 45, 'SC1': 45, 'SQda': 45, 'SQd': 45, 'SQa': 45, 'SQ0': 45, 'AC2': 72, 'AC1': 72,
    'BP4': 72}
    with open(filename, 'r') as file:
        lines = file.readlines()

    atom_section_start = lines.index('[ atoms ]\n') + 1
    atom_section_end = lines.index('[ bonds ]\n')

    for i in range(atom_section_start, atom_section_end):
        line = lines[i].strip().split()
        if len(line) >= 6 and line[1] in atom_mass_dict:
            mass_number = atom_mass_dict[line[1]]
            lines[i] = lines[i].rstrip().replace('; C', '') + ' ' + str(mass_number) + '\n'

    with open(filename, 'w') as file:
        file.writelines(lines)

          
def gen_elp_1(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    name_index = lines.index('[ moleculetype ]\n') + 2
    lines[name_index] = lines[name_index].rstrip().replace('Protein', 'elp_1')
    start_index = lines.index('[ atoms ]\n') + 1
    end_index = lines.index('[ bonds ]\n')

    for i in range(start_index, end_index):
        line = lines[i].strip().split()
        if len(line) == 8:
            lines[i] = lines[i].rstrip() + ' ' + line[1] + '_DUM 0.0000 '+ line[7] +'\n'

    with open(output_file, 'w') as f:
        f.writelines(lines)

def gen_elp_2(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    name_index = lines.index('[ moleculetype ]\n') + 2
    lines[name_index] = lines[name_index].rstrip().replace('Protein', 'elp_2')
    start_index = lines.index('[ atoms ]\n') + 1
    end_index = lines.index('[ bonds ]\n')

    for i in range(start_index, end_index):
        line = lines[i].strip().split()
        if len(line) == 8:
            lines[i] = lines[i].rstrip().replace(line[1], line[1]+'_DUM').replace(line[6], '0.0000')
            lines[i] = lines[i].rstrip()+ ' ' + line[1] + ' '+ line[6] +' '+ line[7] +'\n'

    with open(output_file, 'w') as f:
        f.writelines(lines)
        
