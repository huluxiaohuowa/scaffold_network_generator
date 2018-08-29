from data import graph
from rdkit import Chem
from os import path


def get_sng_from_smiles(smiles):
    mol_graph = graph.get_mol_graph(smiles)
    ls_atom, ls_bond = mol_graph.graph_list_to_list()
    ls_scaffold = mol_graph.ls_mol_from_sng_u()
    ls_mol_atom_idx = []
    for i in range(len(ls_scaffold)):
        ls_atom_idx = []
        for j in ls_atom[i]:
            ls_atom_idx.append(j[1])
        ls_mol_atom_idx.append((Chem.MolToSmiles(ls_scaffold[i]),
                                ls_atom_idx))
    return ls_mol_atom_idx


def output_to_file(smiles, output_file='output.txt'):
    output_dir = path.join(path.dirname(__file__),
                           'datasets',
                           output_file)
    ls_mol_atom_idx = get_sng_from_smiles(smiles)
    with open(output_dir, 'a+') as f:
        for i in ls_mol_atom_idx:
            f.write(smiles + '\t' + i[0] + '\t' + ','.join([str(x) for x in i[1]]) + '\n')



