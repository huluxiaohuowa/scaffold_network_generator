from rdkit import Chem
import re
from copy import deepcopy
from . import data_struct
from os import path

ATOM_SYMBOLS = ['C', 'F', 'I', 'Cl', 'N', 'O', 'P', 'Br', 'S']
BOND_ORDERS = [Chem.BondType.AROMATIC,
               Chem.BondType.SINGLE,
               Chem.BondType.DOUBLE,
               Chem.BondType.TRIPLE]


def get_atom_type(atom):
    atom_symbol = ATOM_SYMBOLS.index(atom.GetSymbol())
    atom_charge = atom.GetFormalCharge()
    atom_hs = atom.GetNumExplicitHs()
    return [atom_symbol, atom_charge, atom_hs]


def get_bond_type(bond):
    return BOND_ORDERS.index(bond.GetBondType())


def tokenize(smiles):
    """Takes a SMILES and return a list of characters/tokens"""
    regex = '(\[[^\[\]]{1,6}\])'
    smiles = smiles.replace('Cl', 'L').replace('Br', 'R')
    char_list = re.split(regex, smiles)
    tokenized = []
    for char in char_list:
        if char.startswith('['):
            tokenized.append(char)
        else:
            chars = list(char)
            tokenized.extend(chars)
    tokenized.append('EOS')
    return tokenized


def _build_atom_types(file_atom_types=path.join(path.dirname(__file__),
                                                'datasets',
                                                'atom_types.txt'),
                      file_chembl=path.join(path.dirname(__file__),
                                            'datasets',
                                            'ChEMBL_cleaned_indexed.smi'
                                            )
                      ):
    atom_types = []
    with open(file_atom_types, 'w') as f:
        with open(file_chembl) as h:
            for line in h:
                if line == '':
                    break
                smiles = line.strip('\n').strip('\r')
                m = Chem.MolFromSmiles(smiles)
                for a in m.GetAtoms():
                    atom_type = get_atom_type(a)
                    if atom_type not in atom_types:
                        atom_types.append(atom_type)
        for atom_symbol in atom_types:
            f.write(str(atom_symbol) + '\n')


def _load_atom_types(file_atom_types=path.join(path.dirname(__file__),
                                               'datasets', 'atom_types.txt'
                                               )):
    atom_types = []
    with open(file_atom_types) as f:
        for line in f:
            atom_types.append(int(x) for x in line.strip('\n'))
    return atom_types


ATOM_TYPES = _load_atom_types()
BOND_TYPES = range(len(BOND_ORDERS))

NUM_ATOM_TYPES = len(ATOM_TYPES)
NUM_BOND_TYPES = len(BOND_TYPES)


def atom_to_index(atom):
    return get_atom_type(atom)


def bond_to_index(bond):
    return get_bond_type(bond)


def index_to_atom(idx):
    return ATOM_SYMBOLS[idx]


def index_to_bond(idx):
    return BOND_ORDERS(idx)

# def index_to_bond(mol, begin_id, end_id, idx):
#     mol.AddBond(begin_id, end_id, BOND_ORDERS[idx])


def label_gen(sssr):
    sssr_copy = deepcopy(sssr)
    sssr_label = deepcopy(sssr)
    sssr_all = []
    for i in sssr_copy:
        sssr_all += i
    for i in range(len(sssr_copy)):
        for j in range(len(sssr_copy[i])):
            if sssr_all.count(sssr_copy[i][j]) == 1:
                sssr_label[i][j] = 1
            else:
                sssr_label[i][j] = 0
    return sssr_copy, sssr_label


def d_sssr_single(sssr):
    d_sssr = []
    sssr_copy, sssr_label = label_gen(sssr)
    for i in range(len(sssr_copy)):
        if sum(sssr_label[i]) == 0:
            continue
        else:
            dd_sssr = []
            for j in range(len(sssr_copy[i])):
                if sssr_label[i][j] == 1:
                    dd_sssr.append(sssr_copy[i][j])
            if len(dd_sssr) > 0:
                d_sssr.append(dd_sssr)
    return d_sssr


def next_sssr(sssr, d_ring):
    sssr_copy = deepcopy(sssr)
    for i in sssr_copy:
        for j in i:
            if j in d_ring:
                sssr_copy.remove(i)
                break
    return sssr_copy


def graph_eq(graph1, graph2):
    if sorted(tuple(graph1.nodes)) == sorted(tuple(graph2.nodes)) and sorted(tuple(graph1.edges)) == sorted(tuple(graph2.edges)):
        return True
    else:
        return False


def get_mol_from_graph(symbol_charge_hs,
                       bond_start_end,
                       sanitize=True
                       ):

    chem = data_struct.get_mol_spec()
    mol = Chem.RWMol(Chem.Mol())
    for atom in symbol_charge_hs:
        mol.AddAtom(chem.index_to_atom(chem.atom_types.index(atom)))
    for bond in bond_start_end:
        chem.index_to_bond(mol, bond[0], bond[1], bond[2])


    if sanitize:
        try:
            mol = mol.GetMol()
            Chem.SanitizeMol(mol)
            return mol
        except:
            return None
    else:
        return None


# def get_mol_from_graph(
#         idx_atom_idx_symbol_charge_hs,
#         bond_start_end,
#         sanitize=True
# ):
#     try:
#         chem = data_struct.get_mol_spec()
#         mol = Chem.RWMol(Chem.Mol())
#         i = 0
#         dic = {}
#         for atom in idx_atom_idx_symbol_charge_hs:
#             mol.AddAtom(chem.index_to_atom(chem.atom_types.index((atom[2],
#                                                                   atom[3],
#                                                                   atom[4]
#                                                                   ))))
#             dic[atom[1]] = atom[0]
#             i += 1
#         for bond in bond_start_end:
#             chem.index_to_bond(mol, dic[bond[0]], dic[bond[1]], bond[2])
#
#         if sanitize:
#             mol = mol.GetMol()
#             Chem.SanitizeMol(mol)
#             return mol
#     except:
#         return None

def get_mol_from_graph_list(ls_ls_atom,
                            ls_ls_bond,
                            sanitize=True
                            ):
    mol_list = [get_mol_from_graph(ls_atom, ls_bond, sanitize) for ls_atom, ls_bond in zip(ls_ls_atom, ls_ls_bond)]
    return mol_list


# noinspection PyArgumentList
# def get_mol_from_graph(X, A, sanitize=True):
#     try:
#         mol = Chem.RWMol(Chem.Mol())
#
#         X, A = X.tolist(), A.tolist()
#         for i, atom_type in enumerate(X):
#             mol.AddAtom(data_struct.get_mol_spec().index_to_atom(atom_type))
#
#         for atom_id1, atom_id2, bond_type in A:
#             data_struct.get_mol_spec().index_to_bond(mol, atom_id1, atom_id2, bond_type)
#     except:
#         return None
#
#     if sanitize:
#         try:
#             mol = mol.GetMol()
#             Chem.SanitizeMol(mol)
#             return mol
#         except:
#             return None
#     else:
#         return mol


# def get_mol_from_graph_list(graph_list, sanitize=True):
#     mol_list = [get_mol_from_graph(X, A, sanitize) for X, A in graph_list]
#     return mol_list
