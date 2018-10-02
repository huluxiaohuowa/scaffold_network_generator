from rdkit import Chem

from data.utils import get_num_lines
from . import graph
from .proto import *
import linecache
from multiprocessing import Pool
from os import path

# input_file = 'input_10.txt'
input_file = 'input.smi'
input_dir = path.join(path.dirname(__file__),
                      'datasets',
                      input_file
                      )


def get_sng_from_smiles(smiles):
    """

    Parameters
    ----------
        smiles: str
            A molecule SMILES

    Returns
        ls_mol_atom_idx: list
            list of tuples
    -------

    """
    mol_graph = graph.get_mol_graph(smiles)
    if any(mol_graph.sssr):
        ls_nh = mol_graph.hydro_nitro
        ls_np = mol_graph.ar_n_plus
        ls_atom, ls_bond = mol_graph.graph_list_to_list()
        ls_scaffold = mol_graph.ls_mol_from_sng_u()
        ls_mol_atom_idx = []
        for i in range(len(ls_scaffold)):
            ls_atom_idx = []
            for j in ls_atom[i]:
                ls_atom_idx.append(j[1])
            ls_mol_atom_idx.append((Chem.MolToSmiles(ls_scaffold[i]),
                                    ls_atom_idx))
        return ls_mol_atom_idx, ls_nh, ls_np
    else:
        return None


def smiles_from_line(idx=0, file=input_dir):
    return linecache.getline(file, idx+1).strip()


def sng_from_line(idx=0, file=input_dir):
    return get_sng_from_smiles(smiles_from_line(idx, file=file))


def sng_from_line_2_queue(idx, q, file=input_dir):
    try:
        sng = sng_from_line(idx, file=file)
        if sng is not None:
            q.put((idx, sng))
    except:
        print(smiles_from_line(idx=idx))


def sng_to_queue(q, processes=30, file=input_dir):
    """

    Parameters
    ----------
    q : multiprocessing.Manager().Queue
        queue
    processes : int
        num of processes
    file : str
        input file dir

    Returns
    -------

    """
    q.put(file)
    p = Pool(processes=processes)
    for i in range(get_num_lines(file)):
        p.apply_async(sng_from_line_2_queue,
                      (i, q,))
    p.close()
    p.join()


def scaffold_smiles_idx(idx, file=path.join(path.dirname(__file__),
                                            'datasets',
                                            'scaffolds.smi')):

    return smiles_from_line(idx, file=file)


def scaffold_mol_idx(idx, file=path.join(path.dirname(__file__),
                                         'datasets',
                                         'scaffolds.smi')):
    return Chem.MolFromSmiles(scaffold_smiles_idx(idx=idx, file=file))


# def protobuf_from_queue(q):
#     """
#     get protobuf message from a queue
#     :param: queue
#     :return: dict{scaffold idx: scaffold smiles}, dict{scaffold idx: (mol idx, [atom idx])}
#     """
#     scaffold_dict = DicIdxScaffolds()
#     dataset = DicScaffoldLs()
#     sng_dict_index = -1
#     file = q.get()
#     for i in range(get_num_lines(file)):
#         if i % 5000 == 0:
#             print(i)
#         mol_index, sng = q.get()
#         for mol_sng_index in range(len(sng)):
#             sng_now, sng_atoms = sng[mol_sng_index]  # 获取当前的scaffold骨架smiles及骨架原子的list
#             if sng_now not in scaffold_dict.dic_scaffold.values():
#                 sng_dict_index += 1
#                 scaffold_dict.dic_scaffold[sng_dict_index] = sng_now
#             sng_pb = TupMolLsatom()
#             sng_pb.idx_mol = mol_index
#             sng_pb.ls_atom.idx_atom.extend(sng_atoms)
#             dataset.idx_scaffold[sng_dict_index].dic_mol_atoms.extend([sng_pb])
#     return scaffold_dict, dataset


def data_from_queue(q, print_step=5000):
    """

    Parameters
    ----------
    q: queue
        queue
    print_step : int

    Returns
    -------
    dic_scaffold
        dict{scaffold smiles: (mol idx, [atom idx])}
    """
    dic_scaffold = DicSmScaffoldLs()
    file = q.get()
    print("Extracting scaffolds from" + file)
    i = 0
    num_lines = get_num_lines(file)
    while True:
        i += 1
        if i % print_step == 0:
            print(i)
        if i > num_lines:
            break
        try:
            mol_index, sng = q.get_nowait()
            for sm_sng, idx_atoms in sng[0]:
                sng_pb = TupMolLsatom()
                sng_pb.idx_mol = mol_index
                sng_pb.ls_atom.idx_atom.extend(idx_atoms)
                sng_pb.ls_nh.idx_atom.extend(sng[1])
                sng_pb.ls_np.idx_atom.extend(sng[2])
                dic_scaffold.smiles_scaffold[sm_sng].dic_mol_atoms.extend([sng_pb])
        except:
            continue

    #
    # for i in range(get_num_lines(file)):
    #     if i % print_step == 0:
    #         print(i)
    #     mol_index, sng = q.get_nowait()
    #     for sm_sng, idx_atoms in sng:
    #         sng_pb = TupMolLsatom()
    #         sng_pb.idx_mol = mol_index
    #         sng_pb.ls_atom.idx_atom.extend(idx_atoms)
    #         dic_scaffold.smiles_scaffold[sm_sng].dic_mol_atoms.extend([sng_pb])
    return dic_scaffold

# def protobuf_from_queue(q):
#     """
#     get protobuf message from a q
#     :param: queue
#     :return: dict{scaffold idx: scaffold smiles} dict{scaffold idx: (mol idx, [atom idx])}
#     """
#     scaffold_dict = DicIdxScaffolds()
#     dataset = DicScaffoldLs()
#     file = q.get()
#     for i in range(get_num_lines(file)):
#         mol_index, sng = q.get()
#         for mol_sng_index in range(len(sng)):
#             sng_now, sng_atoms = sng[mol_sng_index]  # 获取当前的scaffold骨架smiles及骨架原子的list
#             sng_dict_index = -1
#             sng_dict_len = len(scaffold_dict.dic_scaffold)  # scaffold dict 长度
#             for sng_dict_index_tmp in range(sng_dict_len):  # 在scaffold dict里查重
#                 if sng_now == scaffold_dict.dic_scaffold[sng_dict_index_tmp]:
#                     sng_dict_index = sng_dict_index_tmp
#                     break
#             if sng_dict_index == -1:  # 表明该scaffold未出现过
#                 scaffold_dict.dic_scaffold[sng_dict_len] = sng_now
#                 sng_dict_index = sng_dict_len
#             # 初始化Dicmollsatom
#             sng_pb = TupMolLsatom()
#             sng_pb.idx_mol = mol_index
#             sng_pb.ls_atom.idx_atom.extend(sng_atoms)
#             # 将该scaffold对应的 当前正在处理的mol_id和atom_list添加到dataset当中
#             dataset.idx_scaffold[sng_dict_index].dic_mol_atoms.extend([sng_pb])
#     return scaffold_dict, dataset

#
# def protobuf_from_file(file=input_dir):
#     # 初始化scaffold字典和 protobuf数据集
#     scaffold_dict = DicIdxScaffolds()
#     dataset = DicScaffoldLs()
#     # 为方便后面再套循环，先随便初始化一个分子index
#     for mol_index in range(get_num_lines(file)):
#         smiles = linecache.getline(file, mol_index + 1).strip()
#         sng_test = get_sng_from_smiles(smiles)
#         # 遍历正在处理的分子的scaffold
#         for mol_sng_index in range(len(sng_test)):
#             sng_now, sng_atoms = sng_test[mol_sng_index]  # 获取当前的scaffold骨架 smiles 及 骨架原子的list
#             sng_dict_index = -1
#             sng_dict_len = len(scaffold_dict.dic_scaffold)  # scaffold dict 长度
#             for sng_dict_index_tmp in range(sng_dict_len):  # 在scaffold dict里查重
#                 if sng_now == scaffold_dict.dic_scaffold[sng_dict_index_tmp]:
#                     sng_dict_index = sng_dict_index_tmp
#                     break
#             if sng_dict_index == -1:  # 表明该scaffold未出现过
#                 scaffold_dict.dic_scaffold[sng_dict_len] = sng_now
#                 sng_dict_index = sng_dict_len
#             # 初始化Dicmollsatom
#             sng_pb = TupMolLsatom()
#             sng_pb.idx_mol = mol_index
#             sng_pb.ls_atom.idx_atom.extend(sng_atoms)
#             # 将该scaffold对应的 当前正在处理的mol_id和atom_list添加到dataset当中
#             dataset.idx_scaffold[sng_dict_index].dic_mol_atoms.extend([sng_pb])
#     return scaffold_dict, dataset












