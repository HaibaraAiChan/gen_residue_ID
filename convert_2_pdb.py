import os


for filename in os.listdir('./ATP/Ligands/'):
    file = os.path.splitext(filename)
    new_name = file[0] + '.pdb'
    os.rename('./ATP/Ligands/' + filename, './ATP/Ligands/' + new_name)

