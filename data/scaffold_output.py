from rdkit import Chem

from data.utils import get_num_lines, ls_inter
from data import graph
from data.proto import *
import linecache
from multiprocessing import Pool
from os import path
from time import sleep
import json
import psycopg2

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
        # ls_nh = mol_graph.hydro_nitro
        # ls_np = mol_graph.ar_n_plus
        if mol_graph.graph_list_to_list() is not None:
            ls_atom, ls_bond = mol_graph.graph_list_to_list()
        else:
            return None
        ls_scaffold = mol_graph.ls_mol_from_sng_u()
        ls_mol_atom_idx = []
        for i in range(len(ls_scaffold)):
            ls_atom_idx = []
            for j in ls_atom[i]:
                ls_atom_idx.append(j[1])
            ls_nh = ls_inter(mol_graph.hydro_nitro, ls_atom_idx)
            ls_nh_cp = []
            ls_np = ls_inter(mol_graph.ar_n_plus, ls_atom_idx)
            ls_np_cp = []

            cp_a = [[a[0], a[1]] for a in ls_bond[i]]
            ls_ba = [aa for aaa in cp_a for aa in aaa]
            for a_nh in ls_nh:
                if ls_ba.count(a_nh) < 3:
                    ls_nh_cp.append(a_nh)
            for a_np in ls_np:
                if ls_ba.count(a_np) < 3:
                    ls_np_cp.append(a_np)

            ls_mol_atom_idx.append((Chem.MolToSmiles(ls_scaffold[i]),
                                    ls_atom_idx,
                                    ls_nh_cp,
                                    ls_np_cp))
        return ls_mol_atom_idx
    else:
        return None


def smiles_from_line(idx=0, file=input_dir):
    return linecache.getline(file, idx + 1).strip()


def sng_from_line(idx=0, file=input_dir):
    return get_sng_from_smiles(smiles_from_line(idx, file=file))


def sng_from_line_2_queue(idx, q, file=input_dir):
    while True:
        if q.full():
            print('full')
            sleep(0.3)
            continue
        else:
            try:
                sng = sng_from_line(idx, file=file)
                if sng is not None:
                    q.put((idx, sng))
                else:
                    q.put((idx, None))
                break
            except:
                q.put((idx, None))
                print(smiles_from_line(idx=idx, file=file))
                break


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
    """

    Args:
        idx (int):
        file (str):

    Returns:

    """

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


def sql_from_queue(
    q,
    dic_path,
    db_name,
    print_step=5000
):
    dic = json.load(open(dic_path))
    i = 0
    map_id = 0
    sc_key = 0
    file = q.get()
    num_lines = get_num_lines(file)

    # create a table
    conn = psycopg2.connect(**dic)
    cur = conn.cursor()
    cur.execute(
        f'''
        drop table if exists maps;
        drop table if exists scaffolds;
        drop table if exists mols;

        CREATE TABLE scaffolds (id integer PRIMARY KEY, smiles varchar);
        
        CREATE TABLE mols (id integer PRIMARY KEY, smiles varchar);
        
        CREATE TABLE maps (
            id integer PRIMARY KEY, 
            sc_id integer references scaffolds(id), 
            mol_id integer references mols(id),
            ls_atom_idx varchar,
            ls_nh varchar,
            ls_np varchar
        );
        '''
    )
    conn.commit()
    cur.close()
    conn.close()

    while True:
        if i % print_step == 0:
            print(i)
        if i >= num_lines:
            break
        mol_index, sng = q.get()
        i += 1
        mol_smiles = smiles_from_line(mol_index, file).strip()

        conn = psycopg2.connect(**dic)
        cur = conn.cursor()
        cur.execute(
            f'''
            insert into mols(id, smiles) values({mol_index}, '{mol_smiles}');

            '''
        )
        conn.commit()
        cur.close()
        conn.close()

        if sng is not None:
            for i_sng in sng:
                conn = psycopg2.connect(**dic)
                cur = conn.cursor()
                cur.execute(
                    f'''
                    select count(*) from scaffolds where smiles='{i_sng[0]}';
                    '''
                )
                num_ex = cur.fetchone()[0]
                if not num_ex:
                    sc_id = sc_key
                    cur.execute(
                        f'''
                        insert into scaffolds(id, smiles)
                            values ({sc_id}, '{i_sng[0]}');
                        '''
                    )
                    conn.commit()
                    sc_key += 1
                else:
                    cur.execute(
                        f'''
                        select id from scaffolds where smiles='{i_sng[0]}';

                        '''
                    )
                    sc_id = cur.fetchone()[0]
                cur.execute(
                    f'''
                    insert into maps(id, sc_id, mol_id, ls_atom_idx, ls_nh, ls_np)
                        values(
                            {map_id},
                            {sc_id},
                            {mol_index},
                            '{str(i_sng[1])}',
                            '{str(i_sng[2])}',
                            '{str(i_sng[3])}'
                        );
                    '''
                )
                conn.commit()
                map_id += 1
                cur.close()
                conn.close()


def data_from_queue(
    q,
    map_file,
    idx_file,
    print_step=5000,
):
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
    dic_scaffold = DicIdxLs()
    dic_sm_idx = DicSmIdx()
    file = q.get()
    print("Extracting scaffolds from" + file)
    i = 0
    num_lines = get_num_lines(file)
    idx_sc = 0
    while True:        
        if i % print_step == 0:
            print(i)
        if i >= num_lines:
            break
        mol_index, sng = q.get()
        i += 1
        try:
            if sng is not None:
                for sng_i in sng:
                    sng_pb = TupMolLsatom()
                    sng_pb.idx_mol = mol_index
                    sng_pb.ls_atom.idx_atom.extend(sng_i[1])
                    sng_pb.ls_nh.idx_atom.extend(sng_i[2])
                    sng_pb.ls_np.idx_atom.extend(sng_i[3])

                    if sng_i[0] not in dic_sm_idx.sm_sc.keys():
                        dic_sm_idx.sm_sc[sng_i[0]] = idx_sc
                        idx_sc += 1
                    dic_scaffold.smiles_scaffold[dic_sm_idx.sm_sc[sng_i[0]]].dic_mol_atoms.extend([sng_pb])
            else:
                continue
        except:
            continue
    dic_idx_sm = DicIdxSm()
    for k, v in dic_sm_idx.sm_sc.items():
        dic_idx_sm.sm_sc[v] = k

    # return dic_scaffold, dic_idx_sm
    with open(map_file, 'wb') as f:
        f.write(dic_scaffold.SerializeToString())
    with open(idx_file, 'wb') as f:
        f.write(dic_idx_sm.SerializeToString())
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


