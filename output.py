import data
from multiprocessing import Manager, cpu_count
from data import *
import argparse
from os import path
parser = argparse.ArgumentParser()
parser.add_argument("--file_input",
                    help="The location of input file",
                    type=str,
                    default=path.join(path.dirname(__file__),
                                      'data',
                                      'datasets',
                                      'input.txt'))
parser.add_argument("--file_output",
                    help="The location of output file",
                    type=str,
                    default=path.join(path.dirname(__file__),
                                      'data',
                                      'datasets',
                                      'scaffolds.bin'))
parser.add_argument("-n",
                    help="Count of processes",
                    type=int,
                    default=cpu_count()-1)
args = parser.parse_args()
q = Manager().Queue()
data.sng_to_queue(q, processes=n, file=file_input)
dic = data.data_from_queue(q)
with open(file_output, 'wb') as f:
    f.write(dic.SerializeToString())

