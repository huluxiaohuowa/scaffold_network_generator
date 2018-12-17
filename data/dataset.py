import pymysql
import json
import abc
import pandas as pd
from os import path
from collections import Iterable
import random
import itertools
import math

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
        random.shuffle(self.ids_mol_in)

    def __iter__(self):
        if self.mode == 'training':
            for sc_ids_batch in self.divide(self.ids_sc_ex, self.batch_size):
                for idx, sc_id in enumerate(sc_ids_batch):
                    if idx == 0:
                        sc_df = self.dataset[sc_id, : , 1]
                    else:
                        sc_df = pd.concat([sc_df, self.dataset[sc_id, : ,1]])
                yield sc_df
        elif self.mode == 'test':
            for sc_ids_batch in self.divide(self.ids_sc_in, self.batch_size):
                for idx, sc_id in enumerate(sc_ids_batch):
                    if idx == 0:
                        sc_df = self.dataset[sc_id, : , 1]
                    else:
                        sc_df = pd.concat([sc_df, self.dataset[sc_id, : ,1]])
                yield sc_df
        elif self.mode == 'evaluation':
            for sc_ids_batch in self.divide(self.ids_sc_ex, self.batch_size):
                for idx, sc_id in enumerate(sc_ids_batch):
                    if idx == 0:
                        sc_df = self.dataset[sc_id, : ]
                    else:
                        sc_df = pd.concat([sc_df, self.dataset[sc_id, : ]])
                yield sc_df
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

    def mol_to_array(self):
        pass

    def batch_to_array(self):
        pass

    @staticmethod
    def str_to_ls(str_ls):
        exec("return " + str_ls)
