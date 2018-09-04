import data
from multiprocessing import Manager, cpu_count
from data import *
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--file_input", help="The location of input file", type=str)
parser.add_argument("--file_output", help="The location of output file", type=str)
parser.add_argument("-n", help="Count of processes", type=int)
args = parser.parse_args()
# from os import path

if not args.file_input:
    file_input = path.join(path.dirname(__file__),
                           'data',
                           'datasets',
                           'input.txt')
else:
    file_input = args.file_input

if not args.file_output:
    file_output = path.join(path.dirname(__file__),
                           'data',
                           'datasets',
                           'scaffolds.bin')
else:
    file_output = args.file_output

if not args.n:
    n = cpu_count() - 1
else:
    n = args.n

q = Manager().Queue()

data.sng_to_queue(q, processes=n, file=file_input)

dic = data.data_from_queue(q)

print(dic.smiles_scaffold[scaffold_mol_idx(2)])

with open(file_output, 'wb') as f:
    f.write(dic.SerializeToString())

