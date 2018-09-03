
import data
from multiprocessing import Manager
# from os import path


q = Manager().Queue()

data.sng_to_queue(q, processes=30, file='data/datasets/input.txt')

scaffold_dict, dataset = data.protobuf_from_queue(q)
