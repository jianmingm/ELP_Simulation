import os
def extract_numbers(index_file):
    with open(index_file, 'r') as file:
        lines = file.readlines()
    
    # Flags to indicate if we are in the desired section
    in_hydrophilic_section = False
    in_hydrophobic_section = False
    
    # Variables to hold the first and last numbers
    first_hydrophilic_number = None
    last_hydrophobic_number = None
    
    for line in lines:
        # Check if we've reached the hydrophilic section
        if '[ elp_2_hydrophilic ]' in line:
            in_hydrophilic_section = True
        elif '[ elp_2_hydrophobic ]' in line:
            in_hydrophilic_section = False
            in_hydrophobic_section = True
        elif in_hydrophilic_section:
            # Assuming the first line after the section header contains the numbers
            numbers = line.split()
            first_hydrophilic_number = numbers[0]  # Get the first number
            in_hydrophilic_section = False  # Reset flag
        elif in_hydrophobic_section:
            numbers = line.split()
            last_hydrophobic_number = numbers[-1]
            in_hydrophobic_section = False
    
    return first_hydrophilic_number, last_hydrophobic_number

def replace_values_in_file(template_file, replacements):
    with open(template_file, 'r') as file:
        content = file.read()
    
    # Replace the placeholders with actual values
    for key, value in replacements.items():
        content = content.replace(key, value)
    
    with open(template_file, 'w') as file:
        file.write(content)

index_file_path = 'index.ndx'

first_hydrophilic, last_hydrophobic = extract_numbers(index_file_path)

replacements = {
    'XXX': first_hydrophilic,
    'YYY': last_hydrophobic
}

for root, dirs, files in os.walk('MDP'):
    for file in files:
        if file.endswith('.mdp'):
            full_path = os.path.join(root, file)
            replace_values_in_file(full_path, replacements)

for root, dirs, files in os.walk('MDP-r1'):
    for file in files:
        if file.endswith('.mdp'):
            full_path = os.path.join(root, file)
            replace_values_in_file(full_path, replacements)

for root, dirs, files in os.walk('MDP-r2'):
    for file in files:
        if file.endswith('.mdp'):
            full_path = os.path.join(root, file)
            replace_values_in_file(full_path, replacements)


print(f"Replaced values: XXX with {first_hydrophilic} and YYY with {last_hydrophobic}")

os.system('mkdir alchem; mv MDP ./alchem; cp init.gro topol.top index.ndx *itp ./alchem')
os.system('mkdir restraint-1; mv MDP-r1 ./restraint-1/MDP; cp init.gro topol.top index.ndx *itp ./restraint-1')
os.system('mkdir restraint-2; mv MDP-r2 ./restraint-2/MDP; cp init.gro topol.top index.ndx *itp ./restraint-2')

os.system('cd alchem/MDP ; bash gen_mdp.sh')
os.system('cd restraint-1/MDP ; bash gen_mdp.sh')
os.system('cd restraint-2/MDP ; bash gen_mdp.sh')





