from rdkit import Chem
import re
import torch
import numpy as np
import random
import networkx as nx
import multiprocessing as mp
import os
from copy import deepcopy
from data import data_struct
from os import path
from data import data_struct as mol_spec
import linecache

ms = mol_spec.get_default_mol_spec()
ATOM_SYMBOLS = ms.atom_symbols
# []
# with open(path.join(
#     path.dirname(__file__),
#     'datasets',
#     'atom_types.txt')
# ) as f:
#     for line in f.readlines():
#         line = line.strip().split(',')
#         ATOM_SYMBOLS.append(line[0])

BOND_ORDERS = ms.bond_orders
__all__ = [
    'mol_gen', 
    'to_tensor', 
    'get_mol_from_array',
    'str_from_line',
]

# for test only:
# __all__ = ['mol_gen', 'to_tensor', 'collate_fn', 'get_array_from_mol']


# def sample_ordering(mol, num_samples, p=0.9, ms=mol_spec.get_default_mol_spec()):
#     """Sampling decoding routes of a given molecule `mol`
#     Parameters
#     ----------
#         mol : Chem.Mol
#             the given molecule (type: Chem.Mol)
#         num_samples: int
#             number of routes to sample (type: int)
#         p : float
#             1 - probability of making a mistake at each step (type: float)
#         ms : mol_spec.MoleculeSpec

#     Returns
#     -------
#         - route_list[i][j] the index of the atom reached at step j in sample i
#         - step_ids_list[i][j]: the step at which atom j is reach at sample i
#         - logp_list[i]: the log-likelihood value of route i
#     """

#     # build graph
#     atom_types, atom_ranks, bonds= [], [], []
#     for atom in mol.GetAtoms():
#         atom_types.append(ms.get_atom_type(atom))
#     for r in Chem.CanonicalRankAtoms(mol):
#         atom_ranks.append(r)
#     for b in mol.GetBonds():
#         idx_1, idx_2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
#         bonds.append([idx_1, idx_2])

#     # build nx graph
#     graph = nx.Graph()
#     graph.add_nodes_from(range(len(atom_ranks)))
#     graph.add_edges_from(bonds)

#     route_list, step_ids_list, logp_list = [], [], []
#     for i in range(num_samples):
#         step_ids, log_p = traverse_graph(graph, atom_ranks, atom_types, p=p)
#         step_ids_list.append(step_ids)
#         step_ids = np.argsort(step_ids)
#         route_list.append(step_ids)
#         logp_list.append(log_p)

#     # cast to numpy array
#     route_list, step_ids_list = np.array(route_list, dtype=np.int32), np.array(step_ids_list, dtype=np.int32)
#     logp_list = np.array(logp_list, dtype=np.float32)

#     return route_list, step_ids_list, logp_list

# def traverse_graph(graph, atom_ranks, atom_types, current_node=None, step_ids=None, p=0.9, log_p=0.0):
#     """ An recursive function for stochastic traversal of graph `graph`
#     Parameters
#     ----------
#         graph: nx.Graph
#             graph to traverse
#         atom_ranks : list
#             list storing the rank of each atom
#         atom_types: list
#             list storing the type of each atom
#         current_node : int
#             an integer indicating the current node, eq to None if the traversal is not yet started
#         step_ids : int
#             storing the step where each atom is traversed
#         p : float
#             1 - probability of making a mistake at each step (type: float)
#         log_p : float
#             the log-likelihood value

#     Returns
#     -------
#         tuple[list[int], float]
#             step_ids and log_p for the next traversal step
#     """

#     if current_node is None:
#         next_nodes = range(len(atom_ranks)) # the first step: include all atoms
#         next_nodes = sorted(next_nodes, key=lambda _x:atom_ranks[_x]) # sort by atom rank

#         step_ids = [-1, ] * len(next_nodes)
#     else:
#         next_nodes = graph.neighbors(current_node)  # get neighbor nodes
#         next_nodes = sorted(next_nodes, key=lambda _x:atom_ranks[_x])  # sort by atom rank

#         next_nodes = [n for n in next_nodes if step_ids[n] < 0] # filter visited nodes

#     # iterate through neighbors
#     while len(next_nodes) > 0:
#         if len(next_nodes)==1:
#             next_node = next_nodes[0]
#         elif random.random() >= (1 - p):
#             next_node = next_nodes[0]
#             log_p += np.log(p)
#         else:
#             next_node = next_nodes[random.randint(0, len(next_nodes) - 1)]
#             log_p += np.log((1.0 - p) / len(next_nodes))
#         step_ids[next_node] = max(step_ids) + 1
#         _, log_p = traverse_graph(graph, atom_ranks, atom_types, next_node, step_ids, p, log_p)
#         next_nodes = [n for n in next_nodes if step_ids[n] < 0] # filter visited nodes

#     return step_ids, log_p

# def reorder(atom_types, bond_info, route, step_ids):
#     """ Reorder atom and bonds according the decoding route
#     Parameters
#     ----------
#         atom_types : np.ndarray
#             storing the atom type of each atom, size: num_atoms
#         bond_info : np.ndarray
#             storing the bond information, size: num_bonds x 3
#         route : np.ndarray
#             route index
#         step_ids : np.ndarray
#             step index

#     Returns
#     -------
#         tuple[np.ndarray, np.ndarray, np.ndarray]
#             reordered atom_types and bond_info
#     """

#     atom_types, bond_info = np.copy(atom_types), np.copy(bond_info)

#     # sort by step_ids
#     atom_types = atom_types[route]
#     bond_info[:, 0], bond_info[:, 1] = step_ids[bond_info[:, 0]], step_ids[bond_info[:, 1]]
#     max_b, min_b = np.amax(bond_info[:, :2], axis=1), np.amin(bond_info[:, :2], axis=1)
#     bond_info = bond_info[np.lexsort([-min_b, max_b]), :]

#     # separate append and connect
#     max_b, min_b = np.amax(bond_info[:, :2], axis=1), np.amin(bond_info[:, :2], axis=1)
#     is_append = np.concatenate([np.array([True]), max_b[1:] > max_b[:-1]])
#     bond_info = np.concatenate([np.where(is_append[:, np.newaxis],
#                                          np.stack([min_b, max_b], axis=1),
#                                          np.stack([max_b, min_b], axis=1)),
#                                 bond_info[:, -1:]], axis=1)

#     return atom_types, bond_info, is_append

def str_from_line(file, idx):
    """
    Get string from a specific line
    Args:
        idx (int): index of line
        file (string): location of a file

    Returns:

    """
    return linecache.getline(file, idx + 1).strip()


def get_array_from_mol(mol, num_samples=1, p=0.9, ms=mol_spec.get_default_mol_spec()):

    atom_types, bond_info = [], []
    num_atoms, num_bonds = mol.GetNumAtoms(), mol.GetNumBonds()

    for atom_id, atom in enumerate(mol.GetAtoms()):
        atom_types.append(ms.get_atom_type(atom))

    for bond_id, bond in enumerate(mol.GetBonds()):
        bond_info.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), ms.get_bond_type(bond)])

    # shape:
    # atom_types: num_atoms
    # bond_info: num_bonds x 3
    atom_types, bond_info = np.array(atom_types, dtype=np.int32), \
                            np.array(bond_info, dtype=np.int32)

    # sample route
    route_list, step_ids_list, logp = sample_ordering(mol, num_samples, p, ms)
    # initialize paced molecule array data
    mol_array = []

    for sample_id in range(num_samples):
        # get the route and step_ids for the i-th sample
        route_i, step_ids_i = route_list[sample_id, :], step_ids_list[sample_id, :]
        # reorder atom types and bond info
        # note: bond_info [start_ids, end_ids, bond_type]
        atom_types_i, bond_info_i, is_append = reorder(atom_types, bond_info, route_i, step_ids_i)
        # atom type added at each step
        # -1 if the current step is connect
        atom_types_added = np.full([num_bonds,], -1, dtype=np.int32)
        atom_types_added[is_append] = atom_types_i[bond_info_i[:, 1]][is_append]
        # pack into mol_array_i
        # size: num_bonds x 4
        # note: [atom_types_added, start_ids, end_ids, bond_type]
        mol_array_i = np.concatenate([atom_types_added[:, np.newaxis], bond_info_i], axis=-1)
        # add initialization step
        init_step = np.array([[atom_types_i[0], -1, 0, -1]], dtype=np.int32)
        # concat into mol_array
        # size: (num_bonds + 1) x 4
        mol_array_i = np.concatenate([init_step, mol_array_i], axis=0)
        mol_array.append(mol_array_i)

    mol_array = np.stack(mol_array, axis=0)  # num_samples x (num_bonds + 1) x 4

    # Output size:
    # mol_array: num_samples x (num_bonds + 1) x 4
    # logp: num_samples

    return mol_array, logp


def get_mol_from_array(mol_array, sanitize=True, ms=mol_spec.get_default_mol_spec()):
    """
    Converting molecule array to Chem.Mol objects

    Parameters
    ----------
        mol_array : np.ndarray
            The array representation of molecules
            dtype: int, shape: [num_samples, num_steps, 4]
        sanitize : bool
            Whether to sanitize the output molecule, default to True
        ms : mol_spec.MoleculeSpec

    Returns
    -------
        list[Chem.Mol]
            The list of output molecules
    """

    # get shape information
    num_samples, max_num_steps, _ = mol_array.shape

    # initialize the list of output molecules
    mol_list = []

    # loop over molecules
    for mol_id in range(num_samples):
        try:
            mol = Chem.RWMol(Chem.Mol())  # initialize molecule
            for step_id in range(max_num_steps):
                atom_type, begin_ids, end_ids, bond_type = mol_array[mol_id, step_id, :].tolist()
                if end_ids == -1:
                    # if the actions is to terminate
                    break
                elif begin_ids == -1:
                    # if the action is to initialize
                    new_atom = ms.index_to_atom(atom_type)
                    mol.AddAtom(new_atom)
                elif atom_type == -1:
                    # if the action is to connect
                    ms.index_to_bond(mol, begin_ids, end_ids, bond_type)
                else:
                    # if the action is to append new atom
                    new_atom = ms.index_to_atom(atom_type)
                    mol.AddAtom(new_atom)
                    ms.index_to_bond(mol, begin_ids, end_ids, bond_type)
            if sanitize:
                mol = mol.GetMol()
                Chem.SanitizeMol(mol)
        except:
            mol = None
        mol_list.append(mol)

    return mol_list


def collate_fn(batch, k, p, ms=mol_spec.get_default_mol_spec()):  # batch = input [mol, mol, mol ...]
    # get the maximum number of atoms and bonds in this batch
    num_bonds_list = [mol.GetNumBonds() for mol in batch]
    max_num_steps = max(num_bonds_list) + 1

    mol_array, logp = [], []

    for mol_i in batch:
        # size:
        # mol_array_i : k x num_steps_i x 4
        # logp_i: k
        mol_array_i, logp_i = get_array_from_mol(mol_i, k, p, ms)

        # pad to the same length
        num_steps_i = mol_array_i.shape[1]
        mol_array_i = np.pad(mol_array_i,
                             pad_width=[[0, 0], [0, max_num_steps- num_steps_i], [0, 0]],
                             mode='constant', constant_values=-1)

        mol_array.append(mol_array_i)
        logp.append(logp_i)

    mol_array =np.stack(mol_array, axis=0)
    logp =np.stack(logp, axis=0)

    # Output size:
    # mol_array: batch_size x k x max_num_steps x 4
    # logp: batch_size x k

    return mol_array, logp


def mol_gen(directory='datasets', file_names=('ChEMBL_0.smi',),
            batch_size=128, num_epochs=5,
            k=5, p=0.5):
    """
    Iterate and pre-process the molecules in the given smi files
    Parameters
    ----------
        directory : str
            The location of input files, default to 'datasets'
        file_names : tuple[str]
            The dataset files, default to ('ChEMBL_0.smi')
        batch_size : int
            The batch size, default to 128
        num_epochs : int
            The number of epochs, default to 5
        k : int
            The number of route samples for each molecule, default to 5
        p : float
            1 -  the probability of making mistake, default to 0.5
    """

    # preparing data processing workers

    # create output queue
    queue = mp.Queue(len(file_names) * 3)

    # workers
    def target_w(_i):
        # create molecule reader
        _reader = Chem.SmilesMolSupplier(os.path.join(directory, file_names[_i]), nameColumn=-1) #, removeHs=False)
        print('indexing file {} ...'.format(file_names[_i]))
        print('number of conformers: {}'.format(len(_reader)))
        print('index complete !')

        if num_epochs is not None:
            for _ in range(num_epochs):
                _index = range(len(_reader))
                random.shuffle(_index)
                _batch = []
                for _idx in _index:
                    _batch.append(_reader[_idx])
                    if len(_batch) >= batch_size:
                        _outputs = collate_fn(_batch, k, p)
                        queue.put(_outputs)
                        _batch = []
            queue.put(None)
        else:
            while True:
                _index = range(len(_reader))
                random.shuffle(_index)
                _batch = []
                for _idx in _index:
                    _batch.append(_reader[_idx])
                    if len(_batch) >= batch_size:
                        _outputs = collate_fn(_batch, k, p)
                        queue.put(_outputs)
                        _batch = []

    # create and start process
    workers = [mp.Process(target=target_w, args=(idx,)) for idx in range(len(file_names))]
    for w in workers:
        w.start()

    end_counter = 0
    while end_counter < len(file_names):
        record = queue.get()
        if record is None:
            end_counter += 1
        else:
            yield record


def to_tensor(record, device=torch.device('cpu')):
    """Convert numpy record to pytorch tensors
    Parameters
    ----------
        record : tuple[np.ndarray]
            Input numpy records
        device : torch.device
            Device to place the pytorch tensors, default to cpu

    Returns
    -------
        tuple[torch.Tensor]
            The output pytorch tensors
    """
    mol_array, logp = record # type: np.ndarray

    mol_array = torch.tensor(mol_array, dtype=torch.long, device=device)  # type: torch.Tensor
    logp = torch.tensor(logp, dtype=torch.float32, device=device)  # type: torch.Tensor

    # Output size:
    # mol_array: batch_size x k x max_num_steps x 4
    # logp: batch_size x k

    return mol_array, logp

def get_num_lines(input_file):
    for num_lines, line in enumerate(open(input_file, 'r')):
        pass
    return num_lines + 1


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


# 返回sssr_copy = sssr_list （成环原子list的list）
# sssr_label 为相应原子的标记list的list，其中多次参与成环的原子为0，其余为1
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
    if graph1.nodes == graph2.nodes:
        return True
    else:
        return False


def get_mol_from_graph(symbol_charge_hs,
                       bond_start_end,
                       sanitize=True
                       ):

    chem = data_struct.get_default_mol_spec()
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


def get_mol_from_graph_list(ls_ls_atom,
                            ls_ls_bond,
                            sanitize=True
                            ):
    mol_list = [get_mol_from_graph(ls_atom, ls_bond, sanitize) for ls_atom, ls_bond in zip(ls_ls_atom, ls_ls_bond)]
    return mol_list


def id_mol(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol


def ls_inter(ls1, ls2):
    """

    Args:
        ls1 (list):
        ls2 (list):

    Returns:
        ls (list):
    """
    return list(set(ls1).intersection(set(ls2)))