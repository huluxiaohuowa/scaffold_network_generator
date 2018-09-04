
import data
from multiprocessing import Manager
from data import *

# from os import path


q = Manager().Queue()

data.sng_to_queue(q, processes=30, file='data/datasets/input.txt')

dic = data.data_from_queue(q)

print(dic.smiles_scaffold[scaffold_mol_idx(2)])

