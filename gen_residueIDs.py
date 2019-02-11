import os
import numpy as np
import re
from scipy.spatial.distance import cdist
from biopandas.pdb import PandasPdb


def gen_residue_list(ligand_coords, protein_coords, protein_res_ID):
    # distance matrix for the removal of the ligand atoms that are close (<= 5) to any protein atoms
    dist = cdist(ligand_coords, protein_coords, 'euclidean')
    dist = np.array(dist)
    dist_refine = np.argwhere(dist < 5.1)
    residue_list = []
    for i in range(len(dist_refine)):
        idx = dist_refine[i][1]
        residue_list.append(protein_res_ID[idx])
    #
    #
    # for i in range(len(dist)):
    #     if np.all(dist[i, :] > 5.1):
    #         continue
    #     else:
    #         for j in range(len(dist[0])):
    #             if dist[i, j] < 5.1:
    #                 residue_list.append(protein_res_ID[j])

    residue_list, count = np.unique(residue_list, return_counts=True)

    return residue_list


def get_protein_coords(pdb_path):
    ppdb = PandasPdb().read_pdb(pdb_path)
    protein_df = ppdb.df['ATOM']
    protein_coords = np.array([protein_df['x_coord'], protein_df['y_coord'], protein_df['z_coord']]).T
    protein_res_ID = protein_df['residue_number']
    return protein_coords, protein_res_ID


def get_ligand_coords(pdb_path):
    ppdb = PandasPdb().read_pdb(pdb_path)
    ligand_df = ppdb.df['HETATM']
    ligand_coords = np.array([ligand_df['x_coord'], ligand_df['y_coord'], ligand_df['z_coord']]).T
    return ligand_coords


ligand_path = './ATP/Ligands/'
ligand_name = '1a27A00.pdb'
protein_path = './ATP/Proteins/'
protein_name = '1a27A.pdb'

output_path = './ATP/aux_files/'

ligand_coords = get_ligand_coords(ligand_path + ligand_name)
protein_coords, protein_res_ID = get_protein_coords(protein_path + protein_name)
residue_list = gen_residue_list(ligand_coords, protein_coords, protein_res_ID)
# tmp = residue_list.tostring()
# resi_ID_list = np.fromstring(tmp, dtype=int)
# resi_ID_list = np.array_str(resi_ID_list)
# print(type(resi_ID_list))
# # resi_ID_list.replace('\n', ' ')
# output = resi_ID_list[1:-1]


base = os.path.splitext(protein_name)
if not os.path.exists(output_path):
    os.makedirs(output_path)

f = open(output_path + base[0] + '_aux.txt', 'w')
f.write('BindingResidueIDs:')
for i in range(len(residue_list)):
    print(residue_list[i])
    if i < len(residue_list)-1:
        f.write(str(residue_list[i])+" ")
    else:
        f.write(str(residue_list[i])+ "\n")
f.write('BindingSiteCenter:')
f.close()
