import pymysql
import json
import abc
import pandas as pd
from os import path
from collections import Iterable
import random
import itertools
import math
from data import data_struct as mol_spec
import numpy as np
from rdkit import Chem
import networkx as nx
# from torch.utils.data import DataLoader
# from copy import deepcopy


# from utils import Block
# from enum import Enum

dic_mysql = json.load(open(path.join(path.dirname(__file__),'datasets', 'mysql.json'), 'r'))

__all__ = [
    "SQLDataset",
    "SQLSampler",
    "Block",
    "ExcludeBlock",
    "DFTransformer",
]


class Block(object):
    def __init__(self, block_id):
        super().__init__()
        self.num = block_id

    def __repr__(self):
        return self.num.__repr__()


class ExcludeBlock(object):
    def __init__(self, block_id, num_block=5):
        super().__init__()
        if isinstance(block_id, int):
            self.block_id = [block_id]
        elif isinstance(block_id, Iterable):
            self.block_id = block_id
        self.ls_all = list(range(num_block))
        for idx in self.block_id:
            assert idx in self.ls_all, "Block id out of range"
            self.ls_all.remove(idx)
        self.i = -1
    def __repr__(self):
        return self.ls_all.__repr__()
    def __iter__(self):
        return self
    def __next__(self):
        if self.i >= len(self.ls_all) - 1:
            raise StopIteration
        else:
            self.i += 1
            return self.ls_all[self.i]



class Dataset(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        super().__init__()

    @abc.abstractmethod
    def __getitem__(self, idx):
        raise NotImplementedError




class Sampler(object):
    """
    A generator
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        super().__init__()

    @abc.abstractmethod
    def __iter__(self):
        raise NotImplementedError


class Transformer(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        super().__init__()

    @abc.abstractmethod
    def mol_to_array(self):
        raise NotImplementedError

    @abc.abstractmethod
    def batch_to_array(self):
        raise NotImplementedError


class ProtoDataset(Dataset):

    def __init__(self):
        super().__init__()

    def __getitem__(self):
        raise NotImplementedError

    def __len__(self):
        raise NotImplementedError


class SQLDataset(Dataset):

    def __init__(self, mysql_param, scaffold_table, mol_table, map_table, num_fold=5):
        """
        Construction function
        Args:
            param_mysql: dic, dic_mysql = json.load(open('data/mysql.json','r'))
            table: name of table
        """
        super().__init__()
        self.db = pymysql.connect(**mysql_param)
        self.scaffold_table = scaffold_table
        self.mol_table = mol_table
        self.map_table = map_table
        self.num_mol = pd.read_sql(f'SELECT count(*) from {self.mol_table}',con=self.db).iloc[0,0]
        self.num_scaffold = pd.read_sql(f'SELECT count(*) from {self.scaffold_table}',con=self.db).iloc[0,0]
        self.num_map = pd.read_sql(f'SELECT count(*) from {self.map_table}', con=self.db).iloc[0,0]
        self.sc_fold = {}
        self.mol_fold = {}
        self.num_fold = num_fold
        for i in range(self.num_fold):
            self.sc_fold[i] = pd.read_sql(f"select id from {self.scaffold_table} where fold={i};", con=self.db)['id'].tolist()
            self.mol_fold[i] = pd.read_sql(f"select id from {self.mol_table} where fold={i};", con=self.db)['id'].tolist()
    @staticmethod
    def ls_to_string(ls):
        """
        Transform a list to a string which can be recognized by mysql
        eg. ls_to_string([2,3,4]) = "(2,3,4)"
        Args:
            ls: a list of int idx

        Returns:
            string:
        """
        return f"({','.join([str(i) for i in ls])})"

    def __getitem__(self, idx):
        """

        Args:
            idx: idx1,idx2 (int, list, Block, ExluceBlock or slice)
                 eg. ds = SQLDataset(dic_mysql,'scaffolds','molecules','scaffold_molecule')
                     ds[2,212315]
                     ds[2,Block(2)]
                     ds[2,[212315]]
                     ds[Block(2),Block(3)]
        Returns:
            pd.dataframe
        """
        if len(idx) == 2:
            scaffold_id, mol_id = idx
            k = None
        elif len(idx) == 3:
            scaffold_id, mol_id, k = idx
        else:
            raise TypeError
        if isinstance(scaffold_id, Block) and isinstance(mol_id, Block):
            return self._getitem_block_block(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, Block) and isinstance(mol_id, int):
            return self._getitem_block_block(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, Block) and isinstance(mol_id, list):
            return self._getitem_block_list(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, list) and isinstance(mol_id, Block):
            return self._getitem_list_block(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, list) and isinstance(mol_id, list):
            return self._getitem_list_list(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, list) and isinstance(mol_id, int):
            return self._getitem_list_int(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, int) and isinstance(mol_id, Block):
            return self._getitem_int_block(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, int) and isinstance(mol_id, list):
           return self._getitem_int_list(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, int) and isinstance(mol_id, int):
            return self._getitem_int_int(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, slice) and isinstance(mol_id, Block):
            if not any((scaffold_id.start, scaffold_id.stop)):
                return self._getitem_mol_block(mol_id, k)
            elif scaffold_id.start is None:
                scaffold_query = [i for i in range(0,scaffold_id.stop)]
                return self._getitem_list_block(scaffold_query, mol_id, k)
            elif scaffold_id.stop is None:
                scaffold_query = [i for i in range(scaffold_id.start, self.num_scaffold)]
                return self._getitem_list_block(scaffold_query, mol_id, k)
            else:
                scaffold_query = [i for i in range(scaffold_id.start, scaffold_id.stop)]
                return self._getitem_list_block(scaffold_query, mol_id, k)
        elif isinstance(scaffold_id, slice) and isinstance(mol_id, list):
            if not any((scaffold_id.start, scaffold_id.stop)):
                return self._getitem_mol_list(mol_id, k)
            elif scaffold_id.start is None:
                scaffold_query = [i for i in range(0,scaffold_id.stop)]
                return self._getitem_list_list(scaffold_query, mol_id, k)
            elif scaffold_id.stop is None:
                scaffold_query = [i for i in range(scaffold_id.start, self.num_scaffold)]
                return self._getitem_list_list(scaffold_query, mol_id, k)
            else:
                scaffold_query = [i for i in range(scaffold_id.start, scaffold_id.stop)]
                return self._getitem_list_list(scaffold_query, mol_id, k)
        elif isinstance(scaffold_id, slice) and isinstance(mol_id, int):
            if not any((scaffold_id.start, scaffold_id.stop)):
                return self._getitem_mol_int(mol_id, k)
            elif scaffold_id.start is None:
                scaffold_query = [i for i in range(0,scaffold_id.stop)]
                return self._getitem_list_int(scaffold_query, mol_id, k)
            elif scaffold_id.stop is None:
                scaffold_query = [i for i in range(scaffold_id.start, self.num_scaffold)]
                return self._getitem_list_int(scaffold_query, mol_id, k)
            else:
                scaffold_query = [i for i in range(scaffold_id.start, scaffold_id.stop)]
                return self._getitem_list_int(scaffold_query, mol_id, k)
        elif isinstance(scaffold_id, Block) and isinstance(mol_id, slice):
            if not any((mol_id.start, mol_id.stop)):
                return self._getitem_scaffold_block(scaffold_id, k)
            elif mol_id.start is None:
                mol_query = [i for i in range(0,mol_id.stop)]
                return self._getitem_block_list(scaffold_id, mol_query, k)
            elif mol_id.stop is None:
                mol_query = [i for i in range(mol_id.start, self.num_mol)]
                return self._getitem_block_list(scaffold_id, mol_query, k)
            else:
                mol_query = [i for i in range(mol_id.start, mol_id.stop)]
                return self._getitem_block_list(scaffold_id, mol_query, k)
        elif isinstance(scaffold_id, list) and isinstance(mol_id, slice):
            if not any((mol_id.start, mol_id.stop)):
                return self._getitem_scaffold_list(scaffold_id, k)
            elif mol_id.start is None:
                mol_query = [i for i in range(0,mol_id.stop)]
                return self._getitem_list_list(scaffold_id, mol_query, k)
            elif mol_id.stop is None:
                mol_query = [i for i in range(mol_id.start, self.num_mol)]
                return self._getitem_list_list(scaffold_id, mol_query, k)
            else:
                mol_query = [i for i in range(mol_id.start, mol_id.stop)]
                return self._getitem_list_list(scaffold_id, mol_query, k)
        elif isinstance(scaffold_id, int) and isinstance(mol_id, slice):
            if not any((mol_id.start, mol_id.stop)):
                return self._getitem_scaffold_int(scaffold_id, k)
            elif mol_id.start is None:
                mol_query = [i for i in range(0, mol_id.stop)]
                return self._getitem_int_list(scaffold_id, mol_query, k)
            elif mol_id.stop is None:
                mol_query = [i for i in range(mol_id.start, self.num_mol)]
                return self._getitem_int_list(scaffold_id, mol_query, k)
            else:
                mol_query = [i for i in range(mol_id.start, mol_id.stop)]
                return self._getitem_int_list(scaffold_id, mol_query, k)
        elif isinstance(scaffold_id, slice) and isinstance(mol_id, slice):
            if not any((scaffold_id.start, scaffold_id.stop, mol_id.start, mol_id.stop)):
                return self._getitem_all(k)
            else:
                if not any((scaffold_id.start, scaffold_id.stop)):
                    if mol_id.start is None:
                        mol_query = [i for i in range(0, mol_id.stop)]
                    elif mol_id.stop is None:
                        mol_query = [i for i in range(mol_id.start, self.num_mol)]
                    else:
                        mol_query = [i for i in range(mol_id.start, mol_id.stop)]
                    return self._getitem_mol_list(mol_query, k)
                if not any((mol_id.start, mol_id.stop)):
                    if scaffold_id.start is None:
                        scaffold_query = [i for i in range(0, scaffold_id.stop)]
                    elif scaffold_id.stop is None:
                        scaffold_query = [i for i in range(scaffold_id.start, self.num_scaffold)]
                    else:
                        scaffold_query = [i for i in range(scaffold_id.start, scaffold_id.stop)]
                    return self._getitem_scaffold_list(scaffold_query, k)
                else:
                    if mol_id.start is None:
                        mol_query = [i for i in range(0, mol_id.stop)]
                    elif mol_id.stop is None:
                        mol_query = [i for i in range(mol_id.start, self.num_mol)]
                    else:
                        mol_query = [i for i in range(mol_id.start, mol_id.stop)]
                    if scaffold_id.start is None:
                        scaffold_query = [i for i in range(0, scaffold_id.stop)]
                    elif scaffold_id.stop is None:
                        scaffold_query = [i for i in range(scaffold_id.start, self.num_scaffold)]
                    else:
                        scaffold_query = [i for i in range(scaffold_id.start, scaffold_id.stop)]
                    return self._getitem_list_list(scaffold_query, mol_query, k)
        elif isinstance(scaffold_id, ExcludeBlock) and isinstance(mol_id, ExcludeBlock):
            return self._getitem_exblock_exblock(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, ExcludeBlock) and isinstance(mol_id, Block):
            return self._getitem_exblock_block(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, ExcludeBlock) and isinstance(mol_id, list):
            return self._getitem_exblock_list(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, ExcludeBlock) and isinstance(mol_id, int):
            return self._getitem_exblock_int(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, Block) and isinstance(mol_id, ExcludeBlock):
            return self._getitem_block_exblock(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, list) and isinstance(mol_id, ExcludeBlock):
            return self._getitem_list_exblock(scaffold_id, mol_id, k)
        elif isinstance(scaffold_id, int) and isinstance(mol_id, ExcludeBlock):
            return self._getitem_int_exblock(scaffold_id, mol_id, k)

        elif isinstance(scaffold_id, slice) and isinstance(mol_id, ExcludeBlock):
            if not any((scaffold_id.start, scaffold_id.stop)):
                return self._getitem_mol_exblock(mol_id, k)
            elif scaffold_id.start is None:
                scaffold_query = [i for i in range(0,scaffold_id.stop)]
                return self._getitem_list_exblock(scaffold_query, mol_id, k)
            elif scaffold_id.stop is None:
                scaffold_query = [i for i in range(scaffold_id.start, self.num_scaffold)]
                return self._getitem_list_exblock(scaffold_query, mol_id, k)
            else:
                scaffold_query = [i for i in range(scaffold_id.start, scaffold_id.stop)]
                return self._getitem_list_exblock(scaffold_query, mol_id, k)
        elif isinstance(scaffold_id, ExcludeBlock) and isinstance(mol_id, slice):
            if not any((mol_id.start, mol_id.stop)):
                return self._getitem_scaffold_exblock(scaffold_id, k)
            elif mol_id.start is None:
                mol_query = [i for i in range(0,mol_id.stop)]
                return self._getitem_exblock_list(scaffold_id, mol_query, k)
            elif mol_id.stop is None:
                mol_query = [i for i in range(mol_id.start, self.num_mol)]
                return self._getitem_exblock_list(scaffold_id, mol_query, k)
            else:
                mol_query = [i for i in range(mol_id.start, mol_id.stop)]
                return self._getitem_exblock_list(scaffold_id, mol_query, k)
        else:
            raise TypeError


    def _getitem_block_block(self, scaffold_id, mol_id, k):
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold = {mol_id}
                              
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_block_int(self, scaffold_id, mol_id, k):
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id = {mol_id}
                            
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_block_list(self, scaffold_id, mol_id, k):
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id in {mol_id}
                            
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_list_block(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold = {mol_id}
                         
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_list_list(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id in {mol_id}
                           
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_list_int(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        mol_id = str(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id = {mol_id}
                    
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_int_block(self, scaffold_id, mol_id, k):
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold = {mol_id}
                           
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_int_list(self, scaffold_id, mol_id, k):
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id in {mol_id}
                           
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_int_int(self, scaffold_id, mol_id, k):
        scaffold_id, mol_id = str(scaffold_id), str(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id = {mol_id}
                         
                """
        return pd.read_sql(sql, con=self.db)
    def _getitem_exblock_exblock(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold in {mol_id}
                           
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_exblock_block(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        mol_id = str(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold = {mol_id}
                          
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_exblock_int(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id = {mol_id}
                              
                """
        return pd.read_sql(sql, con=self.db)
    def _getitem_exblock_list(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id in {mol_id}
                            
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_int_exblock(self, scaffold_id, mol_id, k):
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold in {mol_id}
                              
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_list_exblock(self, scaffold_id, mol_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold in {mol_id}
                             
                """
        return pd.read_sql(sql, con=self.db)
    def _getitem_block_exblock(self, scaffold_id, mol_id, k):
        scaffold_id = str(scaffold_id)
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold in {mol_id}
                          
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_scaffold_int(self, scaffold_id, k):
        scaffold_id= str(scaffold_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id
                    
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_scaffold_list(self, scaffold_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.id in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id
                         
                """
        if k is not None:
            sql = f"select * from ({sql}) as a order by rand() limit {k}"
            # sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_scaffold_block(self, scaffold_id, k):
        scaffold_id = str(scaffold_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold = {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id
                    
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_scaffold_exblock(self, scaffold_id, k):
        scaffold_id = self.ls_to_string(scaffold_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        where sc.fold in {scaffold_id} 
                        ) AS sc_query ON mol.id = sc_query.molecule_id
                        
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_mol_int(self, mol_id, k):
        mol_id = str(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id) 
                        AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.id = {mol_id}
                             
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_mol_list(self, mol_id, k):
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
               SELECT
                    sc_query.id,
                    sc_query.scaffold_id,
                    sc_query.scaffold_smiles,
                    sc_query.scaffold_fold,
                    mol.id AS molecule_id,
                    mol.smiles AS molecule_smiles,
                    mol.fold AS molecule_fold,
                    sc_query.atom_id_list,
                    sc_query.ls_np,
                    sc_query.ls_nh 
                FROM
                    {self.mol_table} AS mol
                    INNER JOIN (
                SELECT
                    sc.id AS scaffold_id,
                    sc.fold AS scaffold_fold,
                    sc.smiles AS scaffold_smiles,
                    sm.id AS id,
                    sm.molecule_id AS molecule_id,
                    sm.atom_id_list AS atom_id_list,
                    sm.ls_np AS ls_np,
                    sm.ls_nh AS ls_nh 
                FROM
                    {self.scaffold_table} AS sc
                    JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id) 
                    AS sc_query ON mol.id = sc_query.molecule_id 
                WHERE
                    mol.id in {mol_id}

            """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_mol_block(self, mol_id, k):
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold = {mol_id}
                              
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_mol_exblock(self,mol_id, k):
        mol_id = self.ls_to_string(mol_id)
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        ) AS sc_query ON mol.id = sc_query.molecule_id 
                    WHERE
                        mol.fold in {mol_id}
                             
                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)
    def _getitem_all(self, k):
        sql = f"""
                   SELECT
                        sc_query.id,
                        sc_query.scaffold_id,
                        sc_query.scaffold_smiles,
                        sc_query.scaffold_fold,
                        mol.id AS molecule_id,
                        mol.smiles AS molecule_smiles,
                        mol.fold AS molecule_fold,
                        sc_query.atom_id_list,
                        sc_query.ls_np,
                        sc_query.ls_nh 
                    FROM
                        {self.mol_table} AS mol
                        INNER JOIN (
                    SELECT
                        sc.id AS scaffold_id,
                        sc.fold AS scaffold_fold,
                        sc.smiles AS scaffold_smiles,
                        sm.id AS id,
                        sm.molecule_id AS molecule_id,
                        sm.atom_id_list AS atom_id_list,
                        sm.ls_np AS ls_np,
                        sm.ls_nh AS ls_nh 
                    FROM
                        {self.scaffold_table} AS sc
                        JOIN {self.map_table} AS sm ON sc.id = sm.scaffold_id 
                        ) AS sc_query ON mol.id = sc_query.molecule_id

                """
        if k is not None:
            sql = sql + f"ORDER BY RAND() LIMIT {k}"
        return pd.read_sql(sql, con=self.db)



class SQLSampler(Sampler):

    def __init__(self, dataset, block_id, mode, batch_size):
        super().__init__()
        self.dataset = dataset
        self.mode = mode
        self.batch_size = batch_size
        self.num_mol = self.dataset.num_mol
        self.num_scaffold = self.dataset.num_scaffold
        self.block_id = block_id
        # self.num_epoch = num_epoch
        ex_block = list(range(self.dataset.num_fold))
        ex_block.remove(self.block_id)
        self.ex = ex_block
        self.ids_sc_in = self.dataset.sc_fold[self.block_id]
        self.ids_mol_in = self.dataset.mol_fold[self.block_id]
        self.ids_sc_ex = list(itertools.chain.from_iterable([self.dataset.sc_fold[i] for i in self.ex]))
        self.ids_mol_ex = list(itertools.chain.from_iterable([self.dataset.mol_fold[i] for i in self.ex]))
        random.shuffle(self.ids_sc_in)
        random.shuffle(self.ids_sc_ex)
        random.shuffle(self.ids_mol_in)
        random.shuffle(self.ids_mol_ex)
        self.i = 0

    def __iter__(self):
        if self.mode == 'training':

            sc_df = []
            mol_df =[]
            while True:
                sc_id = self.ids_sc_ex[self.i % len(self.ids_sc_ex)]
                mol_id = self.ids_mol_ex[self.i % len(self.ids_mol_ex)]
                sc_df.append(self.dataset[sc_id, : , 1])
                mol_df.append(self.dataset[:, mol_id, 1])
                self.i += 1
                if self.i % int(self.batch_size/2) == 0:
                    sc_df = pd.concat(sc_df)
                    mol_df = pd.concat(mol_df)
                    batch_df = pd.concat([sc_df, mol_df])
                    # batch_df = deepcopy(sc_df)
                    sc_df = []
                    mol_df = []
                    yield batch_df #
        elif self.mode == 'test':
            sc_df = []
            mol_df = []
            while True:
                sc_id = self.ids_sc_in[self.i % len(self.ids_sc_in)]
                mol_id = self.ids_mol_in[self.i % len(self.ids_mol_in)]
                sc_df.append(self.dataset[sc_id, :, 1])
                mol_df.append(self.dataset[:, mol_id, 1])
                self.i += 1
                if self.i % int(self.batch_size / 2) == 0:
                    sc_df = pd.concat(sc_df)
                    mol_df = pd.concat(mol_df)
                    batch_df = pd.concat([sc_df, mol_df])
                    # batch_df = deepcopy(sc_df)
                    sc_df = []
                    mol_df = []
                    yield batch_df  #

        elif self.mode == 'evaluation':
            batch_df = []
            # for sc_ids_batch in self.divide(self.ids_sc_ex, self.batch_size):
            #     for idx, sc_id in enumerate(sc_ids_batch):
            #         if idx == 0:
            #             sc_df = self.dataset[sc_id, : ]
            #         else:
            #             sc_df = pd.concat([sc_df, self.dataset[sc_id, : ]])
            yield batch_df
        else:
            raise  ValueError

    @staticmethod
    def divide(ls, each):
        len_group = len(ls)
        num_group = math.ceil(len_group / each)
        return (ls[each * i:each * (i + 1)] for i in range(num_group))

    @staticmethod
    def get_group_dic(grouped_df):
        dic = {}
        for group in grouped_df:
            if group[0][0] not in dic:
                dic[group[0][0]] = [group[0][1]]
            else:
                dic[group[0][0]].append(group[0][1])
        return dic


class DFTransformer(Transformer):
    def __init__(self,
                 col_sc_smiles=2,
                 col_mol_smiles=5,
                 col_ls_atom=7,
                 col_ls_np=8,
                 col_ls_nh=9):
        super().__init__()
        self.col_sc_smiles = col_sc_smiles
        self.col_mol_smiles = col_mol_smiles
        self.col_ls_atom = col_ls_atom
        self.col_ls_np = col_ls_np
        self.col_ls_nh = col_ls_nh

    @staticmethod
    def smiles_to_mol(smiles):
        return Chem.MolFromSmiles(smiles)

    def mol_to_array(mol,
                     scaffold_nodes,
                     nh_nodes,
                     np_nodes,
                     k, p,
                     ms=mol_spec.get_default_mol_spec()):
        """
        Represent the molecule using `np.ndarray`
        Args:
            mol (Chem.Mol): The input molecule
            scaffold_nodes (Iterable): The location of scaffold represented as `list`/`np.ndarray`
            nh_nodes (Iterable): Nodes with modifications
            np_nodes (Iterable): Nodes with modifications
                dtype - np.int32, shape - [k, num_bonds + 1, 5]
            k (int): The number of importance samples
            p (float): Degree of uncertainty during route sampling, should be in (0, 1)
            ms (mol_spec.MoleculeSpec)
        Returns:
            mol_array (np.ndarray): The numpy representation of the molecule
            logp (np.ndarray): The log-likelihood of each route
                dtype - np.float32, shape - [k, ]
        """
        atom_types, bond_info = [], []
        num_atoms, num_bonds = mol.GetNumAtoms(), mol.GetNumBonds()

        # sample route
        scaffold_nodes = np.array(list(scaffold_nodes), dtype=np.int32)
        route_list, step_ids_list, logp = _sample_ordering(mol, scaffold_nodes, k, p)

        for atom_id, atom in enumerate(mol.GetAtoms()):
            if atom_id in nh_nodes:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
            if atom_id in np_nodes:
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
            atom_types.append(ms.get_atom_type(atom))

        for bond in mol.GetBonds():
            bond_info.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), ms.get_bond_type(bond)])

        # shape:
        # atom_types: num_atoms
        # bond_info: num_bonds x 3
        atom_types, bond_info = np.array(atom_types, dtype=np.int32), \
                                np.array(bond_info, dtype=np.int32)

        # initialize packed molecule array data
        mol_array = []

        for sample_id in range(k):
            # get the route and step_ids for the i-th sample
            route_i, step_ids_i = route_list[sample_id, :], step_ids_list[sample_id, :]

            # reorder atom types and bond info
            # note: bond_info [start_ids, end_ids, bond_type]
            atom_types_i, bond_info_i, is_append = _reorder(atom_types, bond_info, route_i, step_ids_i)

            # atom type added at each step
            # -1 if the current step is connect
            atom_types_added = np.full([num_bonds, ], -1, dtype=np.int32)
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

            # Mark up scaffold bonds
            is_scaffold = np.logical_and(mol_array_i[:, 1] < len(scaffold_nodes),
                                         mol_array_i[:, 2] < len(scaffold_nodes)).astype(np.int32)

            # Concatenate
            # shape: k x (num_bonds + 1) x 5
            mol_array_i = np.concatenate((mol_array_i, is_scaffold[:, np.newaxis]), axis=-1)

            mol_array.append(mol_array_i)

        mol_array = np.stack(mol_array, axis=0)  # num_samples x (num_bonds + 1) x 4

        # Output size:
        # mol_array: k x (num_bonds + 1) x 4
        # logp: k

        return mol_array, logp

    def batch_to_array(self):
        pass


def str_to_ls(str_ls):
    exec("return " + str_ls)

def _sample_ordering(mol, scaffold_nodes, k, p, ms=mol_spec.get_default_mol_spec()):
    """Sampling decoding routes of a given molecule `mol`
    Args:
        mol (Chem.Mol): the given molecule (type: Chem.Mol)
        scaffold_nodes (np.ndarray): the nodes marked as scaffold
        k (int): The number of importance samples
        p (float): Degree of uncertainty during route sampling, should be in (0, 1)
        ms (mol_spec.MoleculeSpec)
    Returns:
        route_list (np.ndarray):
            route_list[i][j] - the index of the atom reached at step j in sample i
        step_ids_list (np.ndarray):
            step_ids_list[i][j] - the step at which atom j is reach at sample i
        logp_list (np.ndarray):
            logp_list[i] - the log-likelihood value of route i
    """
    # build graph
    atom_types, atom_ranks, bonds = [], [], []
    for atom in mol.GetAtoms():
        atom_types.append(ms.get_atom_type(atom))
    for r in Chem.CanonicalRankAtoms(mol):
        atom_ranks.append(r)
    for b in mol.GetBonds():
        idx_1, idx_2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonds.append([idx_1, idx_2])
    atom_ranks = np.array(atom_ranks)

    # build nx graph
    graph = nx.Graph()
    graph.add_nodes_from(range(len(atom_ranks)))
    graph.add_edges_from(bonds)

    route_list, step_ids_list, logp_list = [], [], []
    for i in range(k):
        step_ids, log_p = _traverse(graph, atom_ranks, scaffold_nodes, p)
        step_ids_list.append(step_ids)
        step_ids = np.argsort(step_ids)
        route_list.append(step_ids)
        logp_list.append(log_p)

    # cast to numpy array
    route_list, step_ids_list = np.array(route_list, dtype=np.int32), np.array(step_ids_list, dtype=np.int32)
    logp_list = np.array(logp_list, dtype=np.float32)

    return route_list, step_ids_list, logp_list

# noinspection PyMethodMayBeStatic
def _reorder(atom_types, bond_info, route, step_ids):
    """ Reorder atom and bonds according the decoding route
    Args:
        atom_types (np.ndarray): storing the atom type of each atom, size: num_atoms
        bond_info (np.ndarray): storing the bond information, size: num_bonds x 3
        route (np.ndarray): route index
        step_ids (np.ndarray): step index
    Returns:
        atom_types, bond_info, is_append (np.ndarray): reordered atom_types and bond_info
    """

    atom_types, bond_info = np.copy(atom_types), np.copy(bond_info)

    # sort by step_ids
    atom_types = atom_types[route]
    bond_info[:, 0], bond_info[:, 1] = step_ids[bond_info[:, 0]], step_ids[bond_info[:, 1]]
    max_b, min_b = np.amax(bond_info[:, :2], axis=1), np.amin(bond_info[:, :2], axis=1)
    bond_info = bond_info[np.lexsort([-min_b, max_b]), :]

    # separate append and connect
    max_b, min_b = np.amax(bond_info[:, :2], axis=1), np.amin(bond_info[:, :2], axis=1)
    is_append = np.concatenate([np.array([True]), max_b[1:] > max_b[:-1]])
    bond_info = np.concatenate([np.where(is_append[:, np.newaxis],
                                         np.stack([min_b, max_b], axis=1),
                                         np.stack([max_b, min_b], axis=1)),
                                bond_info[:, -1:]], axis=1)

    return atom_types, bond_info, is_append

def _traverse(graph, atom_ranks, scaffold_nodes, p):
    """ An recursive function for stochastic traversal of the entire `graph`
    Args:
        graph (nx.Graph): The graph to be traversed
        atom_ranks (np.ndarray): A list storing the rank of each atom
        scaffold_nodes (np.ndarray): A list storing the index of atoms in the scaffold
        p (float): Degree of uncertainty during route sampling, should be in (0, 1)
    Returns:
        step_ids (np.ndarray): `step_ids` for the next traversal step
        log_p (np.ndarray): `log_p` for the next traversal step
    """
    step_ids = _traverse_scaffold(graph, atom_ranks, scaffold_nodes, p)
    if len(scaffold_nodes) < len(atom_ranks):
        step_ids, log_p = _traverse_chain(graph, atom_ranks, scaffold_nodes, step_ids, p)
    else:
        log_p = 0.0

    return step_ids, log_p

def _traverse_scaffold(graph, atom_ranks, scaffold_nodes, p,
                       current_node=None, step_ids=None):
    """ An recursive function for stochastic traversal of scaffold in `graph`"""
    # Initialize next_nodes and step_ids (if is None)
    if current_node is None:
        next_nodes = scaffold_nodes  # Initialize as the set of all scaffold nodes
        step_ids = np.full_like(atom_ranks, -1)  # Initialize step_ids as -1
    else:
        next_nodes = np.array(list(graph.neighbors(current_node)), dtype=np.int32)  # get neighbor nodes
        next_nodes = np.intersect1d(next_nodes, scaffold_nodes, assume_unique=True)  # Only scaffold nodes

    next_nodes = next_nodes[np.argsort(atom_ranks[next_nodes])]  # Sort by atom_ranks
    next_nodes = next_nodes[step_ids[next_nodes] == -1]  # Filter visited nodes

    # Iterate through neighbors
    while len(next_nodes) > 0:  # If there are unvisited neighbors
        if len(next_nodes) == 1:  # Only one neighbor is unvisited
            next_node = next_nodes[0]  # Visit this neighbor
        else:
            alpha = (p * np.array([1.0, ] + [0.0, ] * (len(next_nodes) - 1), dtype=np.float32) +
                 (1.0 - p) * np.array([1.0/len(next_nodes), ] * len(next_nodes), dtype=np.float32))
            next_node = np.random.choice(next_nodes, p=alpha)

        step_ids[next_node] = max(step_ids) + 1

        # Proceed to the next step
        step_ids = _traverse_scaffold(graph, atom_ranks, scaffold_nodes, p, next_node, step_ids)
        next_nodes = next_nodes[step_ids[next_nodes] == -1]  # Filter visited nodes

    return step_ids

def _traverse_chain(graph, atom_ranks, scaffold_nodes, step_ids, p,
                    current_node=None, log_p=0.0):
    """ An recursive function for stochastic traversal of side chains in `graph`
    Notes:
        The scaffold should be first traversed using `_traverse_scaffold`
    """
    # Initialize next_nodes
    if current_node is None:  # For the fist step
        next_nodes = set([])  # Initialize next_nodes as an empty set
        for scaffold_node_id in scaffold_nodes:  # Iterate through scaffold nodes
            # Add all nodes directly connected to scaffold nodes
            next_nodes = next_nodes.union(set(graph.neighbors(scaffold_node_id)))
        next_nodes = np.array(list(next_nodes), dtype=np.int32)  # Convert to ndarray
        next_nodes = np.setdiff1d(next_nodes, scaffold_nodes, assume_unique=True)  # Remove all scaffold nodes
    else:
        next_nodes = np.array(list(graph.neighbors(current_node)), dtype=np.int32)  # Get neighbor nodes
        next_nodes = np.setdiff1d(next_nodes, scaffold_nodes, assume_unique=True)  # Remove all scaffold nodes

    next_nodes = next_nodes[np.argsort(atom_ranks[next_nodes])]  # Sort by atom_ranks
    next_nodes = next_nodes[step_ids[next_nodes] == -1]  # Filter visited nodes

    # Iterate through neighbors
    while len(next_nodes) > 0:  # If there are unvisited neighbors
        if len(next_nodes) == 1:  # Only one neighbor is unvisited
            next_node = next_nodes[0]  # Visit this neighbor
            log_p_step = 0.0
        else:
            alpha = (p * np.array([1.0, ] + [0.0, ] * (len(next_nodes) - 1), dtype=np.float32) +
                 (1.0 - p) * np.array([1.0 / len(next_nodes), ] * len(next_nodes), dtype=np.float32))
            next_node_index = np.random.choice(np.arange(len(next_nodes)), p=alpha)
            next_node = next_nodes[next_node_index]
            log_p_step = np.log(alpha[next_node_index])

        # # If scaffold have been iterated
        # if not is_scaffold_iteration:
        #     log_p += log_p_step

        log_p += log_p_step
        step_ids[next_node] = max(step_ids) + 1

        # Proceed to the next step
        step_ids, log_p = _traverse_chain(graph, atom_ranks, scaffold_nodes, step_ids, p, next_node, log_p)
        next_nodes = next_nodes[step_ids[next_nodes] == -1]  # Filter visited nodes

    return step_ids, log_p